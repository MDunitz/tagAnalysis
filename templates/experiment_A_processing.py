"""
Processing script for YEAR-MONTH-DATE DESCRIPTION experiment.

Usage:
    python processing.py
"""

import pandas as pd
from software_module import google_sheets, calibration, plotting
from bokeh.plotting import output_file, save
from bokeh.layouts import column

# Experimental parameters
DATE = ""
DESCRIPTION = ""

# Define a relative path to the data -- don't hardcode a path.
DATA_DIR = "../../data/raw"
OUTPUT_DIR = "output"

# --- Data Import ---
# Option A: From Google Sheet (first run caches locally)
# df = google_sheets.cache_sheet(
#     spreadsheet_id="YOUR_SHEET_ID",
#     output_path=f"{DATA_DIR}/{DATE}_{DESCRIPTION}_raw.csv",
#     worksheet_name="Sheet1",
# )

# Option B: From cached CSV
# df = pd.read_csv(f"{DATA_DIR}/{DATE}_{DESCRIPTION}_raw.csv")


# --- Processing ---
# processed_data = your_processing_here(df)


# --- Save outputs ---
# processed_data.to_csv(
#     f"{OUTPUT_DIR}/{DATE}_{DESCRIPTION}_processed.csv", index=False
# )
