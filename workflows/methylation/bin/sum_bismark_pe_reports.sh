#!/usr/bin/env bash
# Sum Bismark paired-end alignment reports from the chunks of one library
# into one report in Bismark's own layout, so downstream QC that reads a
# per-library PE report (mapping efficiency, methylation by context) works
# unchanged on chunked runs.
#
#   sum_bismark_pe_reports.sh SAMPLE_ID chunk1_PE_report.txt chunk2_PE_report.txt ...
#
# Count lines ("label:<TAB>integer[<TAB>note]") are summed across chunks.
# Percentages are recomputed from the summed counts, not averaged.
set -euo pipefail
sample=$1; shift
[ "$#" -ge 1 ] || { echo "no chunk reports given" >&2; exit 1; }

awk -v sample="$sample" -v nfiles="$#" -F'\t' '
function pct(a, b) { return (b > 0) ? sprintf("%.1f%%", 100 * a / b) : "N/A" }
FNR == 1 { file++ }
# Layout and labels come from the first chunk; counts from every chunk.
file == 1 { line[++n] = $0 }
$2 ~ /^[0-9]+$/ { sum[$1] += $2 }
END {
    for (i = 1; i <= n; i++) {
        split(line[i], f, "\t")
        if (i == 1) {
            print "Bismark report for: " sample " (sum of " nfiles " chunk reports; " substr(line[i], index(line[i], "(version:") + 1)
        } else if (f[2] ~ /^[0-9]+$/) {
            out = f[1] "\t" sum[f[1]]
            if (f[3] != "") out = out "\t" f[3]
            print out
        } else if (f[1] == "Mapping efficiency:") {
            print f[1] "\t" pct(sum["Number of paired-end alignments with a unique best hit:"], sum["Sequence pairs analysed in total:"]) " "
        } else if (f[1] ~ /^C methylated in /) {
            ctx = f[1]; sub(/^C methylated in /, "", ctx); sub(/ context.*/, "", ctx)
            m = sum["Total methylated C\047s in " ctx " context:"]
            u = sum["Total unmethylated C\047s in " ctx " context:"]
            print f[1] "\t" pct(m, m + u)
        } else if (f[1] ~ /^Bismark completed in/) {
            continue
        } else {
            print line[i]
        }
    }
}' "$@"
