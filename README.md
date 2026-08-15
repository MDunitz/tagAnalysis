# Reproducible Science

A template repository for reproducible scientific research. Clone this to start a new project with a standardized directory structure, common analysis utilities, and CI-ready tooling.

## Getting Started

**On a Mac** — if you are running a different OS you may need to update some of the commands.

Clone the repo (replace `$REPONAME` with your new project name):
```
git clone https://github.com/MDunitz/$REPONAME/
```

Initialize the virtual environment and install dependencies:
```
make init
```

Or manually:
```
python -m venv .env
source .env/bin/activate
pip install -r requirements.txt
```

If you want to easily merge updates from this template, set up the template remote once:
```
git remote add template https://github.com/MDunitz/repro.git
```
Then pull template updates anytime with:
```
make template-update
```

### Google Sheets Setup (optional)

The `software_module.google_sheets` module provides authenticated access to Google Sheets. To use it:

1. Create a Google Cloud project and enable the Sheets API
2. Create a service account and download the JSON key file
3. Share your Google Sheet with the service account email
4. Set the `GOOGLE_SHEETS_CREDENTIALS_PATH` environment variable to the key file path

## Make Targets

| Command | Description |
|---------|-------------|
| `make init` | Create venv and install requirements |
| `make new-experiment` | Scaffold a new experiment directory from templates |
| `make fmt` | Run `black` and `flake8` |
| `make tests` | Run `pytest` |
| `make ci` | Run fmt + tests |
| `make md-to-pdf INPUT_FILE=path.md` | Convert markdown to PDF (stored in `pdfs/`, not committed) |
| `make template-update` | Merge latest changes from the template repo |
| `make clean` | Clear the figures directory |

On a Mac you may need to install LaTeX for PDF conversion:
```
brew install --cask basictex
```

## Directory Structure

### **`experiments/`**
Where all experimental information lives, including any executed code — pipelines, scripts, and figure files.
  * **`processing/`**: Code used to *transform* data. This can include parsing text data, image segmentation/filtering, or simulations.
  * **`analysis/`**: Code to *draw conclusions* from an experiment or data set. This may include regression, dimensionality reduction, or calculation of various quantities.
  * **`exploratory/`**: A sandbox where you keep a record of your different approaches to transformation, interpretation, cleaning, or generation of data.
  * **`figures/`**: Code used to generate figures for finished work, presentations, or for any other use.
  * **`lab_notebook/`**: Per-experiment records with date-stamped protocol copies, observations, and deviations.
  * **`protocols/`**: Well-annotated, general descriptions of your experiments. These protocols should be descriptive enough for someone to follow your experiments independently.

### **`data/`**
All raw data collected from your experiments as well as copies of the transformed data from your processing code.
  * **`raw/`**: Raw dataframes imported (typically from Google Sheets) and cached for offline availability.
  * **`transformed/`**: Transformed dataframes cached for plotting and analysis.
  * **`local_data/`**: Larger data files for local analysis without git tracking.

### **`software_module/`**
Reusable code that is called (not executed directly). Functions here should be modular and used the same way from experiment to experiment, so that changes propagate by rerunning processing scripts rather than editing each experiment by hand.

Included utilities:
  * **`google_sheets.py`**: Import Google Sheets as DataFrames with local CSV caching.
  * **`calibration.py`**: Linear calibration curve fitting, LOD/LOQ calculation, dilution correction.
  * **`plotting.py`**: Bokeh figure helpers with project conventions (depth profiles with x-axis on top, auto-annotated calibration curves).

### **`miscellaneous/`**
Files that may not be code, but are important for reproducibility of your findings.
  * **`protocols/`**: General protocol descriptions (see also `experiments/protocols/`).
  * **`materials/`**: Manufacturer information, records of purity, lot and catalog numbers.
  * **`images/`**: Images embedded in miscellaneous documents.
  * **`instruments/`**: Instrument-specific data — links to manuals, explanatory diagrams, care details, calibration records.
  * **`lit_review/`**: Names and links for relevant literature, or concept summaries.

### **`tests/`**
All test suites for your code. *Any custom code you've written should be thoroughly and adequately tested to make sure you know how it is working.*

### **`templates/`**
Blank templates for experiments, processing scripts, and project status tracking. **Never directly copy one experiment's folder and rename the files** — always start from a template (or use `make new-experiment`).

Includes:
  * `experiment.md` — Full experiment description template
  * `experiment_A_README.md` — Lightweight experiment README
  * `experiment_A_processing.py` — Processing script scaffold with `software_module` imports
  * `STATUS.md` — Quick project orientation (current phase, next steps, key files)

## Required Files

1. **`LICENSE`**: A legal protection of your work. *It is important to think deeply about the licensing of your work, and is not a decision to be made lightly. See [this useful site](https://choosealicense.com/) for more information about licensing and choosing the correct license for your project.*

2. **`README.md`**: A descriptive yet succinct description of your research project and information regarding the structure outlined below.

## Acknowledgements

<p xmlns:dct="http://purl.org/dc/terms/" xmlns:vcard="http://www.w3.org/2001/vcard-rdf/3.0#">
  <a rel="license"
     href="http://creativecommons.org/publicdomain/zero/1.0/">
    <img src="http://i.creativecommons.org/p/zero/1.0/88x31.png" style="border-style: none;" alt="CC0" />
  </a>
  <br />
  To the extent possible under law,
  <a rel="dct:publisher"
     href="github.com/gchure/reproducible_research">
    <span property="dct:title">Griffin Chure</span></a>
  has waived all copyright and related or neighboring rights to
  <span property="dct:title">A template for using git as a platform for reproducible scientific research</span>.
This work is published from:
<span property="vcard:Country" datatype="dct:ISO3166"
      content="US" about="github.com/gchure/reproducible_research">
  United States</span>.
</p>

## License Information

This template was originally derived from Griffin Chure's CC0-licensed [`reproducible_research`](https://github.com/gchure/reproducible_research) template. Modifications and additions by Madison Dunitz are licensed under MIT (see `LICENSE`).
