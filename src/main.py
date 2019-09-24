"""
Post GWAS info daily on twitter, alternate between UKBB and FinnGen data.
"""

import logging
from os import getenv
from pathlib import Path

import tweepy

from bot_ukbb import UKBBPoster
from utils import wait


def main():
    """Tweet a GWAS post daily"""
    ukbb = UKBBPoster(
        CORRELATION_FILE,
        H2_FILE,
        MANIFEST_FILE,
        SAVE_FILE_UKBB,
        GWAS_DIR_UKBB,
        GWAS_FILE_SUFFIX_UKBB,
        URI_PREFIX_UKBB
    )
    while True:
        ukbb.tweet()
        wait(hour=8, timezone='America/New_York')


if __name__ == '__main__':
    # Data directory
    DATA_PATH = getenv("DATA_PATH")
    assert DATA_PATH is not None, "DATA_PATH is not set"
    DATA_PATH = Path(DATA_PATH)

    # UKBB files
    CORRELATION_FILE = DATA_PATH / "geno_corr.csv"
    H2_FILE = DATA_PATH / "topline_h2.tsv"
    MANIFEST_FILE = DATA_PATH / "manifest.csv"
    SAVE_FILE_UKBB = DATA_PATH / "posted_ukbb.txt"

    # UKBB GWAS picture files
    GWAS_DIR_UKBB = DATA_PATH / "manhattan"
    GWAS_FILE_SUFFIX_UKBB = "_MF.png"

    # UKBB Google Storage API
    URI_PREFIX_UKBB = getenv("URI_PREFIX_UKBB")
    assert URI_PREFIX_UKBB is not None, "URI_PREFIX_UKBB is not set"
    if not URI_PREFIX_UKBB.endswith('/'):
        URI_PREFIX_UKBB += "/"
    assert URI_PREFIX_UKBB is not None, "URI_PREFIX_UKBB is not set"

    main()
