#!/bin/sh
export PROJECT_ROOT=$(PWD)
export env?=dev
INPUT_FILE ?= "default_batch"
OUTPUT_DIRECTORY_PATH ?= "./pdfs"

# Convert file to pdf
.PHONY: md-to-pdf
md-to-pdf:
	@echo "Input file: $(INPUT_FILE)"
	$(eval BASE_NAME := $(basename $(notdir $(INPUT_FILE))))
	$(eval OUTPUT_FILE := $(OUTPUT_DIRECTORY_PATH)/$(BASE_NAME).pdf)
	@echo "Output file path: $(OUTPUT_FILE)"
	pandoc -V geometry:margin=1in --from=gfm -t pdf "$(INPUT_FILE)" -o "$(OUTPUT_FILE)" --pdf-engine=xelatex

# E501 is line too long
.PHONY: fmt
fmt:
	black .
	flake8 --ignore=E501 .

# Non-modifying check used by CI. fmt above is for local use.
.PHONY: check-fmt
check-fmt:
	black --check .
	flake8 --ignore=E501 .

.PHONY: tests
tests:
	pytest tests/unit_tests/ tests/integration_tests/ -v

.PHONY: test-cov
test-cov:
	pytest tests/ --cov=software_module --cov-report=term-missing

.PHONY: ci
ci: check-fmt tests

.PHONY: clean
clean:
	echo "Clearing figures folder"
	rm -rf figures
	mkdir figures

## Setup (local)
.PHONY: install-reqs
install-reqs:
	pip install -r requirements.txt

.PHONY: init
init: clean
	python3 -m venv .env
	.env/bin/pip install --upgrade pip
	.env/bin/pip install -r requirements.txt
	@if [ "$(vscode)" = "True" ]; then\
	    .env/bin/python -m ipykernel install --user --name=$(project);\
	fi

.PHONY: template-update
template-update:
	git checkout main
	- git branch -d repro-repo-template-update
	git pull --rebase origin main
	git fetch template
	git checkout -b repro-repo-template-update
	git merge template/main --allow-unrelated-histories

# Create a new experiment from the template
.PHONY: new-experiment
new-experiment:
	@read -p "Experiment date (YYYY-MM-DD): " date; \
	read -p "Brief description (snake_case): " desc; \
	dir="experiments/processing/$${date}_$${desc}"; \
	mkdir -p "$${dir}/output"; \
	cp templates/experiment_A_processing.py "$${dir}/processing.py"; \
	cp templates/experiment_A_README.md "$${dir}/README.md"; \
	cp templates/STATUS.md "$${dir}/STATUS.md"; \
	sed -i'' -e "s/PROJECT_NAME/$${date} $${desc}/" "$${dir}/STATUS.md"; \
	echo "Created $${dir}/"
