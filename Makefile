.PHONY: test replot-relabund redecontam

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

# Rerun decontamination (and its plots) on persisted pipeline outputs, without
# primer removal, DADA2 or taxonomy. Needs R + decontam.
# Usage: make redecontam INPUT_DIR=path/to/pipeline/output [LIBRARY_BATCHES=path/to/library_batches.csv]
redecontam:
ifndef INPUT_DIR
	$(error INPUT_DIR is required, e.g. make redecontam INPUT_DIR=path/to/output LIBRARY_BATCHES=batches.csv)
endif
	python -m tag_analysis.redecontam $(INPUT_DIR) $(if $(strip $(LIBRARY_BATCHES)),--library-batches $(LIBRARY_BATCHES))
