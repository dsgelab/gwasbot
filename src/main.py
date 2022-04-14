"""
Post GWAS info daily on twitter, alternate between UKBB and FinnGen data.
"""

import logging
from os import getenv
from pathlib import Path

from requests.exceptions import RequestException
from tweepy.error import TweepError

from bot_bbj import BBJPoster
from bot_finngen import FGPoster
from bot_ukbb import UKBBPoster
from bot_metsim import METSIMPoster
from utils import wait


def poster_turn(iter):
    turns = [
        "METSIM", "UKBB", "FG", "BBJ"
    ]
    select = iter % len(turns)
    return turns[select]


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

    bbj = BBJPoster(
        SAVE_FILE_BBJ,
        FAILURE_FILE_BBJ,
        MANIFEST_FILE_BBJ,
        GWAS_DIR_BBJ,
    )

    metsim = METSIMPoster(
        SAVE_FILE_METSIM,
        FAILURE_FILE_METSIM,
        MANIFEST_FILE_METSIM,
        GWAS_DIR_METSIM,
    )

    do_wait = True
    iter = 0
    while True:
        if do_wait:
            wait(hour=8, timezone='America/New_York')

        turn = poster_turn(iter)
        if turn == "UKBB":
            poster = ukbb
        elif turn == "FG":
            poster = fg
        elif turn == "BBJ":
            poster = bbj
        elif turn == "METSIM":
            poster = metsim

        pheno = poster.get_pheno()
        try:
            poster.tweet(pheno)
        except Exception as exc:  # Catch anything except KeyboardInterrupt
            logging.error(f"Could not tweet: {exc}")
            poster.mark_failure(pheno)
            # Skip the current phenotype and retry immediately with another one
            do_wait = False
        else:
            poster.mark_posted(pheno)
            do_wait = True

        iter += 1


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
    GWAS_GS_PREFIX = "gs://finngen-public-data-r6/summary_stats/finngen_R6_"

    # BioBank Japan files
    SAVE_FILE_BBJ = DATA_PATH / "posted_bbj.txt"
    FAILURE_FILE_BBJ = DATA_PATH / "failure_bbj.txt"
    MANIFEST_FILE_BBJ = DATA_PATH / "images_bbj.js"
    GWAS_DIR_BBJ = DATA_PATH / "manhattan_BBJ"

    # METSIM files
    SAVE_FILE_METSIM = DATA_PATH / "posted_metsim.txt"
    FAILURE_FILE_METSIM = DATA_PATH / "failure_metsim.txt"
    MANIFEST_FILE_METSIM = DATA_PATH / "images_metsim.js"
    GWAS_DIR_METSIM = DATA_PATH / "manhattan_METSIM"

    main()
