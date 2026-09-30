# AREE task runner. Run `make help` for the list.
.PHONY: help install test lint demo demo-clean register harmonize meta cards crosswalk intake intake-check real-study real-pool dashboard docs app clean

DEMO_STUDIES := GIGAS_HEAT01 GIGAS_OA02 GIGAS_PATH03 GIGAS_SAL04 GIGAS_LARV05 GIGAS_GROW06

help:
	@echo "AREE make targets:"
	@echo "  install      pip install the package with dev+app extras"
	@echo "  test         run the pytest suite"
	@echo "  lint         run ruff over src and tests"
	@echo "  demo         run the full demo pipeline (register -> harmonize -> meta -> cards)"
	@echo "  demo-clean   remove generated reports/ and the registry index, then run demo"
	@echo "  crosswalk    rebuild the real NCBI/UniProt identifier crosswalk (~380 MB download)"
	@echo "  intake       regenerate the real study's result tables from the published source"
	@echo "  intake-check verify the committed result tables still reproduce (no writes)"
	@echo "  real-study   register + harmonize the real study (HESSER2024_VCOR) against the real crosswalk"
	@echo "  real-pool    register + harmonize all three real studies against the real crosswalk"
	@echo "  dashboard    build docs/dashboard/data.json from current outputs, then render the site"
	@echo "  docs         render the Quarto documentation site"
	@echo "  app          launch the Streamlit interface"
	@echo "  clean        remove generated reports/ and caches"

install:
	pip install -e ".[dev,app,intake]"

test:
	pytest -q

lint:
	ruff check src tests

register:
	@for f in registry/studies/GIGAS_*.yaml; do aree register-study "$$f" --update; done
	@aree list-studies

harmonize:
	@for sid in $(DEMO_STUDIES); do aree harmonize --study "$$sid"; done

meta:
	aree meta-analyze --feature-type gene

cards:
	aree build-evidence-cards
	aree top-candidates

demo: register harmonize meta cards
	@echo ""
	@echo "Demo complete. Outputs:"
	@echo "  reports/evidence/evidence_table.tsv"
	@echo "  reports/meta_analysis/"
	@echo "  reports/evidence_cards/"
	@echo "  reports/top_candidates_summary.md"
	@echo "  reports/manifests/"

demo-clean: clean demo

crosswalk:
	aree build-crosswalk

REAL_CROSSWALK := data/reference/crosswalk/mgigas_gene_id_crosswalk.tsv
REAL_INTAKE := data/studies/HESSER2024_VCOR/intake.yaml

intake:
	aree intake-supplementary $(REAL_INTAKE)

intake-check:
	aree intake-supplementary $(REAL_INTAKE) --check

real-study: intake-check
	@test -f $(REAL_CROSSWALK) || (echo "Real crosswalk missing. Run: make crosswalk" && exit 1)
	AREE_CROSSWALK=$(REAL_CROSSWALK) aree register-study registry/studies/HESSER2024_VCOR.yaml --update
	AREE_CROSSWALK=$(REAL_CROSSWALK) aree harmonize --study HESSER2024_VCOR

# Everything the Pages workflow does to populate the dashboard: the real
# studies whose standardized result tables are committed under data/studies/.
# CALLA2026_OSHV (raw_reanalysis, comparisons still pending) is harmonized one
# committed comparison at a time, as docs/first_raw_reanalysis.md describes.
REAL_STUDIES := HESSER2024_VCOR DELISLE2020_OSHV_TEMP
RAW_STUDY := CALLA2026_OSHV

real-pool:
	@test -f $(REAL_CROSSWALK) || (echo "Real crosswalk missing. Run: make crosswalk" && exit 1)
	@for sid in $(REAL_STUDIES) $(RAW_STUDY); do \
		AREE_CROSSWALK=$(REAL_CROSSWALK) aree register-study registry/studies/$$sid.yaml --update; \
	done
	@for sid in $(REAL_STUDIES); do \
		AREE_CROSSWALK=$(REAL_CROSSWALK) aree harmonize --study $$sid; \
	done
	@for f in data/studies/$(RAW_STUDY)/$(RAW_STUDY)_*_dge_standardized.tsv; do \
		cid=$$(basename $$f _dge_standardized.tsv); cid=$${cid#$(RAW_STUDY)_}; \
		AREE_CROSSWALK=$(REAL_CROSSWALK) aree harmonize --study $(RAW_STUDY) --comparison $$cid --input $$f; \
	done

dashboard:
	aree build-dashboard
	quarto render docs

docs:
	quarto render docs

app:
	streamlit run app/main.py

clean:
	rm -rf reports/evidence reports/meta_analysis reports/evidence_cards reports/manifests
	rm -f reports/top_candidates_summary.md
	rm -rf .pytest_cache .ruff_cache
	find . -type d -name __pycache__ -prune -exec rm -rf {} +
