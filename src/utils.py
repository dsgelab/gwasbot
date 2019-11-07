"""
Functions useful for both the FinnGen and UKBB bots.
"""

import logging
from datetime import date, datetime
from io import BytesIO
from os import getenv
from time import sleep

import pytz
from dateutil.relativedelta import relativedelta
from google.cloud import storage


logging.basicConfig(level=logging.INFO)


def check_done_pheno(save_file, failure_file):
    """Check already posted GWAS ids from save files"""
    logging.info("Checking already posted and failed phenotypes")

    try:
        with open(save_file) as f:
            posted = f.readlines()
            posted = list(map(lambda line: line.rstrip(), posted))
    except FileNotFoundError:
        logging.info(f"Nothing already posted in {save_file}.")
        posted = []

    try:
        with open(failure_file) as f:
            failed = f.readlines()
            failed = list(map(lambda line: line.rstrip(), failed))
    except FileNotFoundError:
        logging.info(f"No failed phenotypes to load in {failure_file}")
        failed = []

    return posted + failed


def gcloud_download(gs_path):
    """Download a file from Google Cloud Storage into an in-memory buffer."""
    logging.info(f"Downloading {gs_path} into memory")

    buffer = BytesIO()
    client = storage.Client()
    client.download_blob_to_file(gs_path, buffer)

    # Reset file position to be re-usable later by Pandas
    buffer.seek(0)

    return buffer


def wait(hour, timezone):
    """Wait until tomorrow 8am East coast time"""
    tomorrow = date.today() + relativedelta(days=1)
    tz = pytz.timezone(timezone)
    next_time = datetime(tomorrow.year, tomorrow.month, tomorrow.day, hour, 0, tzinfo=tz)

    now = datetime.now(tz=tz)
    seconds = (next_time - now).seconds

    logging.info(f"Sleeping for {seconds} s until next post")
    sleep(seconds)
