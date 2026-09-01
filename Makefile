.PHONY: test replot-relabund

test:
	pytest -q

# Regenerate relative-abundance stackbar plots from persisted pipeline
# outputs (counts + taxonomy CSVs) without rerunning the pipeline.
# Usage: make replot-relabund INPUT_DIR=path/to/pipeline/output
replot-relabund:
ifndef INPUT_DIR
	$(error INPUT_DIR is required, e.g. make replot-relabund INPUT_DIR=path/to/output)
endif
	python -m tag_analysis.replot $(INPUT_DIR)
