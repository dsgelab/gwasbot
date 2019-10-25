"""
Post GWAS info daily on twitter, alternate between UKBB and FinnGen data.
"""

import logging
from os import getenv
from pathlib import Path

from google.cloud import exceptions

from bot_ukbb import UKBBPoster
from bot_finngen import FGPoster
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

    fg = FGPoster(
        ENDPOINTS_FILE,
        SAVE_FILE_FG,
        GWAS_DIR_FG,
        GWAS_FILE_SUFFIX_FG,
        FINNGEN_URI_PREFIX,
    )

    turn = "UKBB"
    do_wait = True
    while True:
        if do_wait:
            wait(hour=8, timezone='America/New_York')

        try:
            if turn == "UKBB":
                ukbb.tweet()
            elif turn == "FG":
                fg.tweet()
        except exceptions.NotFound:
            # Skip the current phenotype and retry immediately with another one
            do_wait = False
        else:
            do_wait = True

        turn = "FG" if turn == "UKBB" else "UKBB"


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


    # FinnGen files
    ENDPOINTS_FILE = DATA_PATH / "endpoint_definitions.tsv"
    SAVE_FILE_FG = DATA_PATH / "posted_fg.txt"

    # FinnGen GWAS plot files
    GWAS_DIR_FG = DATA_PATH / "plots_fg"
    GWAS_FILE_SUFFIX_FG = "_MF.png"

    # FinnGen Google Storage
    FINNGEN_URI_PREFIX = getenv("FINNGEN_URI_PREFIX")
    assert FINNGEN_URI_PREFIX is not None, "FINNGEN_URI_PREFIX is not set"
    if not FINNGEN_URI_PREFIX.endswith("/"):
        FINNGEN_URI_PREFIX += "/"


    # Twitter API
    CONSUMER_KEY = getenv("CONSUMER_KEY")
    CONSUMER_SECRET = getenv("CONSUMER_SECRET")
    ACCESS_TOKEN = getenv("BOT_ACCESS_TOKEN")
    ACCESS_SECRET = getenv("BOT_ACCESS_SECRET")

    assert CONSUMER_KEY is not None, "CONSUMER_KEY is not set"
    assert CONSUMER_SECRET is not None, "CONSUMER_SECRET is not set"
    assert ACCESS_TOKEN is not None, "ACCESS_TOKEN is not set"
    assert ACCESS_SECRET is not None, "ACCESS_SECRET is not set"


    main()
