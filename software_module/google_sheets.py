"""
Google Sheets data import utilities.

Provides authenticated access to Google Sheets for pulling experimental
data into pandas DataFrames. Supports both service account and OAuth
authentication flows.

Setup:
    1. Create a Google Cloud project and enable the Sheets API
    2. Create a service account and download the JSON key file
    3. Share your Google Sheet with the service account email
    4. Set GOOGLE_SHEETS_CREDENTIALS_PATH env var to the key file path,
       or pass it directly to get_client()
"""

import os
import pandas as pd
import gspread
from google.oauth2.service_account import Credentials

SCOPES = [
    "https://www.googleapis.com/auth/spreadsheets.readonly",
    "https://www.googleapis.com/auth/drive.readonly",
]


def get_client(credentials_path=None):
    """Return an authenticated gspread client.

    Parameters
    ----------
    credentials_path : str, optional
        Path to service account JSON key file.
        Falls back to GOOGLE_SHEETS_CREDENTIALS_PATH env var.
    """
    path = credentials_path or os.environ.get("GOOGLE_SHEETS_CREDENTIALS_PATH")
    creds = Credentials.from_service_account_file(path, scopes=SCOPES)
    return gspread.authorize(creds)


def sheet_to_dataframe(
    spreadsheet_id, worksheet_name=None, worksheet_index=0, credentials_path=None
):
    """Pull a Google Sheet worksheet into a pandas DataFrame.

    Parameters
    ----------
    spreadsheet_id : str
        The spreadsheet ID from the Google Sheets URL.
    worksheet_name : str, optional
        Name of the worksheet tab. If None, uses worksheet_index.
    worksheet_index : int
        Index of the worksheet tab (default 0 = first tab).
    credentials_path : str, optional
        Path to service account JSON key file.

    Returns
    -------
    pd.DataFrame
    """
    client = get_client(credentials_path)
    spreadsheet = client.open_by_key(spreadsheet_id)

    if worksheet_name:
        worksheet = spreadsheet.worksheet(worksheet_name)
    else:
        worksheet = spreadsheet.get_worksheet(worksheet_index)

    records = worksheet.get_all_records()
    return pd.DataFrame(records)


def cache_sheet(
    spreadsheet_id,
    output_path,
    worksheet_name=None,
    worksheet_index=0,
    credentials_path=None,
):
    """Pull a Google Sheet and save it as a CSV for offline use.

    Parameters
    ----------
    spreadsheet_id : str
        The spreadsheet ID from the Google Sheets URL.
    output_path : str
        Path to save the CSV file (typically in data/raw/).
    worksheet_name : str, optional
        Name of the worksheet tab.
    worksheet_index : int
        Index of the worksheet tab.
    credentials_path : str, optional
        Path to service account JSON key file.

    Returns
    -------
    pd.DataFrame
    """
    df = sheet_to_dataframe(
        spreadsheet_id, worksheet_name, worksheet_index, credentials_path
    )
    df.to_csv(output_path, index=False)
    print(f"Cached {len(df)} rows to {output_path}")
    return df
