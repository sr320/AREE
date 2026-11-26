"""Offline tests for `aree fetch-fastq` (src/intake/fastq_stage.py), with and without the mirror.

A local HTTP server with byte-range support stands in for ENA and, under /raw/<ACCESSION>/, for the gannet mirror.
Mirror copies are different bytes from ENA's, as files converted from NCBI's
.sra are. No network
access beyond localhost.
"""
from __future__ import annotations

import csv
import hashlib
import http.server
import threading
import urllib.error

import pytest
from click.testing import CliRunner

import aree.cli as cli_mod
import intake.fastq_stage as fq
from aree.cli import main
from intake.fastq_mirror import LEDGER_COLUMNS
from intake.fastq_stage import (
    FastqError,
    ManifestFile,
    apply_mirror,
    fetch_all,
    fetch_one,
    load_ledger,
    plan,
    read_manifest,
    record_staged,
)

PAYLOAD = {f"SRR1_{m}.fastq.gz": bytes(range(256)) * (40 + m) for m in (1, 2)}


class _Handler(http.server.BaseHTTPRequestHandler):
    ignore_range = False
    requests: list[str] = []
    routes: dict[str, bytes] = {}   # full path -> body, checked before PAYLOAD by file name

    def do_GET(self):  # noqa: N802
        name = self.path.rsplit("/", 1)[-1]
        type(self).requests.append(f"{self.path if self.path in self.routes else name} "
                                   f"{self.headers.get('Range') or ''}".strip())
        data = self.routes.get(self.path, PAYLOAD.get(name) if not self.path.startswith("/raw/") else None)
        if data is None:
            self.send_error(404)
            return
        rng = self.headers.get("Range")
        if rng and not self.ignore_range:
            start = int(rng.split("=")[1].split("-")[0])
            body = data[start:]
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(data) - 1}/{len(data)}")
        else:
            body = data
            self.send_response(200)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def log_message(self, *args):
        pass


@pytest.fixture
def server(monkeypatch):
    # Clusters that route web traffic through a proxy (klone sets http_proxy,
    # with 127.0.0.1 missing from no_proxy) would send these requests to the
    # proxy, which cannot reach the test server.
    for var in ("no_proxy", "NO_PROXY"):
        monkeypatch.setenv(var, "127.0.0.1,localhost")
    _Handler.requests = []
    _Handler.ignore_range = False
    _Handler.routes = {}
    httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), _Handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_port}/vol1/fastq/SRR1", _Handler
    httpd.shutdown()


def _files(base, md5_override=None):
    return [ManifestFile("SRR1", name, f"{base}/{name}",
                         (md5_override or {}).get(name, hashlib.md5(data).hexdigest()), len(data))
            for name, data in PAYLOAD.items()]


def test_download_verify_and_skip(server, tmp_path):
    base, handler = server
    files = _files(base)
    res = fetch_all(files, tmp_path, jobs=2)
    assert res == {"downloaded": 2, "skipped": 0, "failed": {}}
    for name, data in PAYLOAD.items():
        assert (tmp_path / name).read_bytes() == data
        assert (tmp_path / f"{name}.md5").read_text().startswith(hashlib.md5(data).hexdigest())
        assert not (tmp_path / f"{name}.part").exists()
    handler.requests.clear()
    assert fetch_all(files, tmp_path)["skipped"] == 2
    assert handler.requests == []
    assert plan(files, tmp_path)["remaining_bytes"] == 0


def test_resume_from_partial(server, tmp_path):
    base, handler = server
    f = _files(base)[0]
    (tmp_path / f"{f.file}.part").write_bytes(PAYLOAD[f.file][:1000])
    (tmp_path / f"{f.file}.part.src").write_text(f.https_url + "\n")
    assert plan([f], tmp_path)["remaining_bytes"] == f.bytes - 1000
    assert fetch_one(f, tmp_path) == "downloaded"
    assert handler.requests == [f"{f.file} bytes=1000-"]
    assert (tmp_path / f.file).read_bytes() == PAYLOAD[f.file]


def test_server_ignoring_range_restarts_cleanly(server, tmp_path):
    base, handler = server
    handler.ignore_range = True
    f = _files(base)[0]
    (tmp_path / f"{f.file}.part").write_bytes(PAYLOAD[f.file][:1000])
    (tmp_path / f"{f.file}.part.src").write_text(f.https_url + "\n")
    assert fetch_one(f, tmp_path) == "downloaded"
    assert (tmp_path / f.file).read_bytes() == PAYLOAD[f.file]


def test_md5_mismatch_deletes(server, tmp_path):
    base, _ = server
    f = _files(base, {"SRR1_1.fastq.gz": "0" * 32})[0]
    with pytest.raises(FastqError, match="Deleted"):
        fetch_one(f, tmp_path)
    assert not any(tmp_path.iterdir())


def test_unverified_file_in_place_is_checked_not_refetched(server, tmp_path):
    base, handler = server
    f = _files(base)[0]
    (tmp_path / f.file).write_bytes(PAYLOAD[f.file])
    assert fetch_one(f, tmp_path) == "skipped"
    assert handler.requests == [] and (tmp_path / f"{f.file}.md5").exists()
    # A wrong file in place is replaced.
    (tmp_path / f"{f.file}.md5").unlink()
    (tmp_path / f.file).write_bytes(b"x" * f.bytes)
    assert fetch_one(f, tmp_path) == "downloaded"
    assert (tmp_path / f.file).read_bytes() == PAYLOAD[f.file]


def test_retries_then_succeeds_and_gives_up(tmp_path):
    data = b"ACGT" * 100
    f = ManifestFile("SRR9", "SRR9_1.fastq.gz", "ftp://host/SRR9_1.fastq.gz", hashlib.md5(data).hexdigest(), len(data))
    calls = []

    def flaky(url, part, offset):
        calls.append((url, offset))
        if len(calls) == 1:
            with open(part, "wb") as out:
                out.write(data[:150])
            raise urllib.error.URLError("reset")
        with open(part, "ab") as out:
            out.write(data[offset:])

    assert fetch_one(f, tmp_path, fetch=flaky, sleep=lambda s: None) == "downloaded"
    assert calls == [("https://host/SRR9_1.fastq.gz", 0), ("https://host/SRR9_1.fastq.gz", 150)]

    def dead(url, part, offset):
        raise ConnectionError("down")
    g = ManifestFile("SRR9", "SRR9_2.fastq.gz", "ftp://host/SRR9_2.fastq.gz", "a" * 32, 10)
    with pytest.raises(FastqError, match="after 4 attempts"):
        fetch_one(g, tmp_path, fetch=dead, sleep=lambda s: None)


def test_read_manifest(tmp_path):
    path = tmp_path / "m.tsv"
    path.write_text("run_accession\tmate\turl\tmd5\tbytes\n"
                    "SRR1\t1\tftp://h/SRR1_1.fastq.gz\tabc\t10\n"
                    "SRR2\t1\tftp://h/SRR2_1.fastq.gz\tdef\t20\n")
    files = read_manifest(path)
    assert [f.file for f in files] == ["SRR1_1.fastq.gz", "SRR2_1.fastq.gz"]
    assert [f.bytes for f in files] == [10, 20] and files[0].https_url == "https://h/SRR1_1.fastq.gz"
    assert [f.run_accession for f in read_manifest(path, ["SRR2"])] == ["SRR2"]
    with pytest.raises(FastqError, match="not in the manifest"):
        read_manifest(path, ["SRR3"])
    path.write_text(path.read_text().replace("\tabc\t", "\t\t"))
    with pytest.raises(FastqError, match="no MD5"):
        read_manifest(path)


def test_record_staged(server, tmp_path):
    base, _ = server
    files = _files(base)
    fetch_one(files[0], tmp_path / "raw")
    record = tmp_path / "study" / "fastq_staged.tsv"
    record_staged(files, tmp_path / "raw", record, machine="raven", today="2026-10-23")
    rows = list(csv.DictReader(record.open(), delimiter="\t"))
    assert [(r["file"], r["machine"]) for r in rows] == [("SRR1_1.fastq.gz", "raven")]
    fetch_one(files[1], tmp_path / "raw")
    record_staged(files, tmp_path / "raw", record, machine="mac", today="2026-10-24")
    rows = list(csv.DictReader(record.open(), delimiter="\t"))
    assert [(r["file"], r["machine"]) for r in rows] == [("SRR1_1.fastq.gz", "raven"), ("SRR1_2.fastq.gz", "mac")]


def test_second_downloader_into_same_dir_is_refused(server, tmp_path):
    import fcntl
    base, _ = server
    with open(tmp_path / ".aree-fetch.lock", "w") as held:
        fcntl.flock(held, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(FastqError, match="already downloading"):
            fetch_all(_files(base), tmp_path)
    assert fetch_all(_files(base), tmp_path)["downloaded"] == 2


# --- with the mirror -----------------------------------------------------------

ACC = "PRJNA1"
MIRROR = {name: b"converted-" + data[::-1] for name, data in PAYLOAD.items()}


def _md5(b: bytes) -> str:
    return hashlib.md5(b).hexdigest()


def _ledger(rows):
    lines = ["\t".join(LEDGER_COLUMNS)]
    for r in rows:
        lines.append("\t".join(str(r.get(c, "")) for c in LEDGER_COLUMNS))
    return ("\n".join(lines) + "\n").encode()


def _mirror_row(name, *, md5=None, ena_md5=None):
    return {"run_accession": "SRR1", "mate": name[-10], "file": name, "url": "https://s3/sra/SRR1",
            "md5": md5 or _md5(MIRROR[name]), "bytes": len(MIRROR[name]), "source": "ncbi_sra",
            "ena_md5": ena_md5 or _md5(PAYLOAD[name]), "mirrored_on": "2026-10-23"}


def _serve_mirror(handler, base, rows, files=MIRROR):
    root = base.split("/vol1")[0]
    handler.routes[f"/raw/{ACC}/mirror_manifest.tsv"] = _ledger(rows)
    for name, data in files.items():
        handler.routes[f"/raw/{ACC}/{name}"] = data
    return f"{root}/raw/{ACC}"


def test_apply_mirror_counts():
    files = _files("http://ena")
    ledger = {files[0].file: _mirror_row(files[0].file),
              files[1].file: _mirror_row(files[1].file, ena_md5="0" * 32)}
    out, counts = apply_mirror(files, ledger, "http://gannet/raw/X/")
    assert counts == {"mirror": 1, "not_mirrored": 0, "stale": 1}
    assert out[0].source == "mirror" and out[0].url == f"http://gannet/raw/X/{files[0].file}"
    assert out[0].fallback is files[0] and out[0].ena_md5 == files[0].md5
    assert out[1] is files[1]
    assert apply_mirror(files, None, "http://g")[1] == {"mirror": 0, "not_mirrored": 2, "stale": 0}


def test_copies_from_mirror_and_records_source(server, tmp_path):
    base, handler = server
    url = _serve_mirror(handler, base, [_mirror_row(n) for n in PAYLOAD])
    files, counts = apply_mirror(_files(base), load_ledger(f"{url}/mirror_manifest.tsv"), url)
    assert counts["mirror"] == 2
    assert fetch_all(files, tmp_path / "raw")["downloaded"] == 2
    for name in PAYLOAD:
        assert (tmp_path / "raw" / name).read_bytes() == MIRROR[name]
    assert all(r.startswith(f"/raw/{ACC}/") for r in handler.requests if "manifest" not in r)

    record = tmp_path / "fastq_staged.tsv"
    record_staged(files, tmp_path / "raw", record, machine="raven", today="2026-10-23")
    rows = list(csv.DictReader(record.open(), delimiter="\t"))
    assert {r["source"] for r in rows} == {"mirror"}
    assert rows[0]["md5"] == _md5(MIRROR[rows[0]["file"]]) and rows[0]["ena_md5"] == _md5(PAYLOAD[rows[0]["file"]])


def test_broken_mirror_copy_falls_back_to_ena(server, tmp_path):
    base, handler = server
    name = "SRR1_1.fastq.gz"
    url = _serve_mirror(handler, base, [_mirror_row(name, md5="f" * 32)], {name: MIRROR[name]})
    files, _ = apply_mirror(_files(base), load_ledger(f"{url}/mirror_manifest.tsv"), url)
    res = fetch_all(files, tmp_path, log=lambda m: None)
    assert res == {"downloaded": 2, "skipped": 0, "failed": {}}
    assert (tmp_path / name).read_bytes() == PAYLOAD[name]
    assert (tmp_path / f"{name}.md5").read_text().split()[2] == "ena"


def test_file_staged_from_ena_is_not_refetched_once_mirrored(server, tmp_path):
    base, handler = server
    ena = _files(base)
    fetch_all(ena, tmp_path)
    url = _serve_mirror(handler, base, [_mirror_row(n) for n in PAYLOAD])
    files, _ = apply_mirror(ena, load_ledger(f"{url}/mirror_manifest.tsv"), url)
    handler.requests.clear()
    assert fetch_all(files, tmp_path)["skipped"] == 2
    assert handler.requests == []
    record = tmp_path / "rec.tsv"
    record_staged(files, tmp_path, record)
    assert {r["source"] for r in csv.DictReader(record.open(), delimiter="\t")} == {"ena"}


def test_partial_from_ena_is_not_resumed_from_mirror(server, tmp_path):
    base, handler = server
    name = "SRR1_1.fastq.gz"
    ena = _files(base)[0]
    (tmp_path / f"{name}.part").write_bytes(PAYLOAD[name][:500])
    (tmp_path / f"{name}.part.src").write_text(ena.https_url + "\n")
    url = _serve_mirror(handler, base, [_mirror_row(name)], {name: MIRROR[name]})
    files, _ = apply_mirror([ena], load_ledger(f"{url}/mirror_manifest.tsv"), url)
    handler.requests.clear()
    assert fetch_all(files, tmp_path)["downloaded"] == 1
    assert handler.requests == [f"/raw/{ACC}/{name}"]   # from the start, no Range
    assert (tmp_path / name).read_bytes() == MIRROR[name]


def test_no_ledger_is_none(server):
    base, _ = server
    assert load_ledger(f"{base.split('/vol1')[0]}/raw/NOPE/mirror_manifest.tsv") is None


def _cli_repo(tmp_path, base, monkeypatch):
    repo = tmp_path / "repo"
    study = repo / "data" / "studies" / "X2020_Y"
    study.mkdir(parents=True)
    with open(study / "fastq_manifest.tsv", "w") as fh:
        fh.write("run_accession\tmate\turl\tmd5\tbytes\tsource\n")
        for f in _files(base):
            fh.write(f"SRR1\t1\t{f.url}\t{f.md5}\t{f.bytes}\tfastq_ftp\n")
    monkeypatch.setattr(cli_mod, "REPO_ROOT", repo)
    monkeypatch.setenv("AREE_RAW", str(tmp_path / "raw"))
    return repo, study


def test_cli_stages_from_ena_without_a_mirror(server, tmp_path, monkeypatch):
    base, _ = server
    _, study = _cli_repo(tmp_path, base, monkeypatch)
    args = ["fetch-fastq", "--study", "X2020_Y"]

    res = CliRunner().invoke(main, [*args, "--dry-run"])
    assert res.exit_code == 0, res.output
    assert "No mirror location known" in res.output
    assert "to download      : 2" in res.output and not (tmp_path / "raw").exists()

    with monkeypatch.context() as mp:
        mp.setattr(fq.shutil, "disk_usage", lambda p: type("U", (), {"free": 10})())
        res = CliRunner().invoke(main, args)
        assert res.exit_code != 0 and "not enough free space" in res.output

    res = CliRunner().invoke(main, args)
    assert res.exit_code == 0, res.output
    assert "2 downloaded" in res.output
    assert (tmp_path / "raw" / "X2020_Y" / "SRR1_2.fastq.gz").read_bytes() == PAYLOAD["SRR1_2.fastq.gz"]
    rows = list(csv.DictReader((study / "fastq_staged.tsv").open(), delimiter="\t"))
    assert {r["source"] for r in rows} == {"ena"} and len(rows) == 2


def test_cli_sources(server, tmp_path, monkeypatch):
    base, handler = server
    _cli_repo(tmp_path, base, monkeypatch)
    name = "SRR1_1.fastq.gz"
    url = _serve_mirror(handler, base, [_mirror_row(name)], {name: MIRROR[name]})
    args = ["fetch-fastq", "--study", "X2020_Y", "--mirror-url", url]

    res = CliRunner().invoke(main, [*args, "--source", "mirror", "--dest", str(tmp_path / "m")])
    assert res.exit_code != 0 and "not usable from the mirror" in res.output

    res = CliRunner().invoke(main, [*args, "--dest", str(tmp_path / "auto")])
    assert res.exit_code == 0, res.output
    assert "1 file(s) from the mirror, 1 not mirrored" in res.output
    assert (tmp_path / "auto" / name).read_bytes() == MIRROR[name]
    assert (tmp_path / "auto" / "SRR1_2.fastq.gz").read_bytes() == PAYLOAD["SRR1_2.fastq.gz"]

    res = CliRunner().invoke(main, [*args, "--source", "ena", "--dest", str(tmp_path / "ena")])
    assert res.exit_code == 0 and "Mirror" not in res.output
    assert (tmp_path / "ena" / name).read_bytes() == PAYLOAD[name]
