"""
Post GWAS info daily on twitter, alternate between UKBB and FinnGen data.
"""

import logging
from os import getenv
from pathlib import Path

from requests.exceptions import RequestException
from tweepy.error import TweepError

from bot_ukbb import UKBBPoster
from bot_finngen import FGPoster
from utils import wait


def main():
    """Tweet a GWAS post daily"""
    ukbb = UKBBPoster(
        H2_FILE,
        MANIFEST_FILE,
        SAVE_FILE_UKBB,
        FAILURE_FILE_UKBB,
        GWAS_DIR_UKBB,
        GWAS_FILE_SUFFIX_UKBB,
    )

    fg = FGPoster(
        SAVE_FILE_FG,
        FAILURE_FILE_FG,
        METADATA_FILE_FG,
        GWAS_DIR_FG,
        GWAS_GS_PREFIX,
    )

    turn = "UKBB"
    do_wait = True
    while True:
        if do_wait:
            wait(hour=8, timezone='America/New_York')
        if turn == "UKBB":
            poster = ukbb
        elif turn == "FG":
            poster = fg

        pheno = poster.get_pheno()
        try:
            poster.tweet(pheno)
        except (RequestException, TweepError) as exc:
            logging.error(f"Could not tweet: {exc}")
            poster.mark_failure(pheno)
            # Skip the current phenotype and retry immediately with another one
            do_wait = False
        except Exception as exc:  # Catch anything except KeyboardInterrupt
            # Don't retry posting so we can see that something has
            # gone wrong and investigate, whilst still allowing the
            # bot to post the next day.
            logging.error(f"Uncaught exception: {type(exc)} {exc}")
        else:
            poster.mark_posted(pheno)
            do_wait = True

        turn = "FG" if turn == "UKBB" else "UKBB"


if __name__ == '__main__':
    # Data directory
    DATA_PATH = getenv("DATA_PATH")
    assert DATA_PATH is not None, "DATA_PATH is not set"
    DATA_PATH = Path(DATA_PATH)

    # UKBB files
    H2_FILE = DATA_PATH / "topline_h2.tsv"
    MANIFEST_FILE = DATA_PATH / "Pan_UKBB_manifest_2020-09-07.csv"
    SAVE_FILE_UKBB = DATA_PATH / "posted_ukbb.txt"
    FAILURE_FILE_UKBB = DATA_PATH / "failure_ukbb.txt"

    # UKBB GWAS picture files
    GWAS_DIR_UKBB = DATA_PATH / "manhattan_UKBB_trans"
    GWAS_FILE_SUFFIX_UKBB = "_trans_MF.png"

    # FinnGen files
    SAVE_FILE_FG = DATA_PATH / "posted_fg.txt"
    FAILURE_FILE_FG = DATA_PATH / "failed_fg.txt"
    METADATA_FILE_FG = DATA_PATH / "images_finngen.js"
    GWAS_DIR_FG = DATA_PATH / "manhattan_FINNGEN"
    GWAS_GS_PREFIX = "gs://finngen-public-data-r4/summary_stats/finngen_R4_"

    main()
