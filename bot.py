"""
GWAS Bot: post one UKBB GWAS with info very 24h to twitter.
"""

import logging
from csv import excel_tab
from datetime import date, datetime
from io import BytesIO
from os import getenv
from pathlib import Path
from time import sleep

import numpy as np
import pandas as pd
import pytz
import tweepy
from dateutil.relativedelta import relativedelta
from google.cloud import storage


# Data files
CORRELATION_FILE = "geno_corr.csv"
H2_FILE = "topline_h2.tsv"
MANIFEST_FILE = "manifest.csv"
SAVE_FILE = "posted.txt"

# GWAS picture files
GWAS_DIR = Path("manhattan")
GWAS_FILE_SUFFIX = "_MF.png"

# Twitter API
CONSUMER_KEY = getenv("CONSUMER_KEY")
CONSUMER_SECRET = getenv("CONSUMER_SECRET")
ACCESS_TOKEN = getenv("BOT_ACCESS_TOKEN")
ACCESS_SECRET = getenv("BOT_ACCESS_SECRET")

# Google Storage API
UKBB_URI_PREFIX = getenv("UKBB_URI_PREFIX")
if not UKBB_URI_PREFIX.endswith('/'):
    UKBB_URI_PREFIX += "/"


# Setup logging
logging.basicConfig(level=logging.INFO)


def main():
    """Post a GWAS to Twitter once a day"""
    posted = check_posted(SAVE_FILE)
    _corr, phenos = load_correlations_phenos(CORRELATION_FILE)
    h2 = load_h2(H2_FILE, phenos)
    manifest = load_manifest(MANIFEST_FILE, phenos)

    # Pick a phenotype at random from the GWAS set, minus already posted ones
    to_post = list(set(phenos) - set(posted))
    np.random.shuffle(to_post)
    for pheno in to_post:
        post = build_post(pheno, manifest, h2)
        text = format_post(post)

        try:
            tweet(pheno, text)
        except tweepy.error.TweepError as exc:
            logging.error(f"Could not post tweet for pheno '{pheno}': {exc}")
        else:
            mark_posted(pheno, to_post)
            logging.info(f"Successfully posted pheno '{pheno}'")

        # Wait until tomorrow for next post
        wait()


def check_posted(filename):
    """Check already posted GWAS ids from a saved text file"""
    logging.info("Checking already posted phenotypes")
    try:
        with open(filename) as f:
            posted = f.readlines()
            posted = list(map(lambda line: line.rstrip(), posted))

    except FileNotFoundError:
        logging.info("No GWAS already posted.")
        posted = []

    return posted


def load_correlations_phenos(filename):
    """Extract the phenocodes of interest and their correlations"""
    logging.debug("Loading file: genetic correlation + phenos")
    corr = pd.read_csv(
        filename,
        usecols=[
            "p2",
            "p1",
            "rg",
            "se"
        ],
        dtype={
            "p2": "category",
            "p1": "category",
            "rg": "float64",
            "se": "float64"
        }
    )

    phenos = corr["p2"].unique()

    return (corr, phenos)


def load_h2(filename, phenos):
    """Load and clean heritability file"""
    logging.debug("Loading file: heritability")
    # Get h2 values from topline file
    topline = pd.read_csv(
        filename,
        dialect=excel_tab,
        usecols=[
            "phenotype",
            "h2_liability",
            "h2_liability_se",
        ],
    )

    # Rename topline pheno ending in _irnt
    irnt = topline["phenotype"].str.endswith("_irnt")
    topline.loc[irnt, "phenotype"] = topline.loc[irnt, "phenotype"].str.rstrip("_irnt")
    topline = topline.astype({
        "phenotype": "category"
    })

    # Filtering-out phenos for topline
    keep = topline["phenotype"].isin(phenos)
    topline = topline[keep]

    return topline


def load_manifest(filename, phenos):
    """Load phenos links and misc info"""
    logging.debug("Loading file: manifest")
    manifest = pd.read_csv(
        filename,
        usecols=[
            "Phenotype Code",
            "Phenotype Description",
            "UK Biobank Data Showcase Link",
            "File",
            "Dropbox File",
            "Sex",
        ]
    )

    manifest.rename(
        columns={
            "Phenotype Code": "phenotype",
            "Phenotype Description": "pheno_desc",
            "UK Biobank Data Showcase Link": "ukbb",
            "File": "filename",
            "Dropbox File": "dropbox",
            "Sex": "sex",
        },
        inplace=True
    )

    # Remove NaN phenotypes
    nan_phenos = manifest["phenotype"].isna()
    manifest = manifest[~ nan_phenos]

    # Keep only "both sexes"
    both_sexes = manifest["sex"] == "both_sexes"
    manifest = manifest[both_sexes]

    # Remove traits ending in _raw
    raw_phenos = manifest["phenotype"].str.endswith("_raw")
    manifest = manifest[~ raw_phenos]

    # Rename Manifest traits ending in _irnt,
    irnt_phenos = manifest["phenotype"].str.endswith("_irnt")
    manifest.loc[irnt_phenos, "phenotype"] = manifest.loc[irnt_phenos, "phenotype"].str.rstrip("_irnt")

    # Merge with the traits of interest
    keep = manifest["phenotype"].isin(phenos)
    manifest = manifest[keep]

    # Remove duplicated phenos
    dups = manifest["phenotype"].duplicated(keep="last")
    manifest = manifest[~ dups]

    return manifest


def build_post(pheno, manifest, h2):
    """Gather all the info we need for the twitter post"""
    logging.info("Building twitter post")
    # Get phenotype info from manifest
    pheno_info = manifest.loc[
        manifest.phenotype == pheno,
        ["pheno_desc", "dropbox", "ukbb"]
    ].iloc[0]  # go from a pd.Series to a scalar value

    # Genetic correlations
    ukbbrg = f"https://ukbb-rg.hail.is/rg_summary_{pheno}.html"

    # Heritability
    heritability = h2_ci(pheno, h2)

    # Top SNP
    snp = top_snp(pheno, manifest)
    gnomad_snp = snp.replace(":", "-")
    gnomad = f"https://gnomad.broadinstitute.org/variant/{gnomad_snp}"

    post = {
        "pheno": pheno_info["pheno_desc"],
        "ukbb_link": pheno_info["ukbb"],
        "download": pheno_info["dropbox"],
        "ukbbrg_link": ukbbrg,
        "heritability": heritability,
        "top_snp": snp,
        "gnomad_link": gnomad,
    }
    return post


def h2_ci(pheno, h2):
    """Get the heritability h2 95% confidence interval"""
    logging.info("Getting h2 values")
    h2data = h2[h2.phenotype == pheno].iloc[0]
    ci_min = h2data.h2_liability - 1.96 * h2data.h2_liability_se
    ci_max = h2data.h2_liability + 1.96 * h2data.h2_liability_se
    return {
        "h2": h2data.h2_liability,
        "ci_min": ci_min,
        "ci_max": ci_max,
    }


def top_snp(pheno, manifest):
    """Find the top SNP after downloading the phenotype data"""
    logging.info("Finding top SNP")
    data = download_gwas(pheno, manifest)

    # Load data
    df = pd.read_csv(
        data,
        usecols=["variant", "pval"],
        dialect=excel_tab,
        compression="gzip"
    )

    # Find top SNP
    snp = df.sort_values(by="pval").iloc[0].variant

    return snp


def download_gwas(pheno, manifest):
    """Download GWAS from google cloud storage to memory"""
    logging.info("Downloading GWAS")
    # Get the filename for this phenotype
    filename = manifest.loc[manifest.phenotype == pheno, "filename"].iloc[0]

    # Build the google storage API client
    client = storage.Client()
    blob = UKBB_URI_PREFIX + filename

    # Download the GWAS file in memory
    buffer = BytesIO()
    client.download_blob_to_file(blob, buffer)

    # Reset file position to be re-usable later by Pandas
    buffer.seek(0)

    return buffer


def format_post(post):
    """Textual representation of the post"""
    logging.debug("Formatting twitter post")
    # Shorten the phenocode name if too long for the tweet
    limit = 65  # based on a few tests
    if len(post['pheno']) > limit:
        pheno = post['pheno'][:limit - 2] + "…"
    else:
        pheno = post['pheno']

    # Don't show description link if it is "nan"
    if pd.isna(post['ukbb_link']):
        desc = ""
    else:
        desc = f"\n📚 Description {post['ukbb_link']}\n"

    text = f"""{pheno}
{desc}
⬇️ Download {post['download']}

🧬 Genetic correlations {post['ukbbrg_link']}

👪 Heritability
{post['heritability']['h2']:.2f} [{post['heritability']['ci_min']:.2f}, {post['heritability']['ci_max']:.2f}]

📊 GWAS top hit
{post['top_snp']} {post['gnomad_link']}"""

    return text


def tweet(pheno, text):
    """Make a post to twitter"""
    logging.info("Posting to twitter")
    auth = tweepy.OAuthHandler(CONSUMER_KEY, CONSUMER_SECRET)
    auth.set_access_token(ACCESS_TOKEN, ACCESS_SECRET)
    api = tweepy.API(auth)

    # Upload GWAS image
    filepath = GWAS_DIR / (pheno + GWAS_FILE_SUFFIX)
    gwas = api.media_upload(str(filepath))

    # Post to twitter
    api.update_status(text, media_ids=[gwas.media_id_string])


def mark_posted(pheno, to_post):
    """Record that the given pheno was posted to twitter"""
    logging.info("Marking phenotype as posted")
    with open(SAVE_FILE, "a") as f:
        f.write(f"{pheno}\n")
    to_post.remove(pheno)


def wait():
    """Wait until tomorrow 8am East coast time"""
    tomorrow = date.today() + relativedelta(days=1)
    tz = pytz.timezone('America/New_York')
    next_time = datetime(tomorrow.year, tomorrow.month, tomorrow.day, 8, 0, tzinfo=tz)

    now = datetime.now(tz=tz)
    seconds = (next_time - now).seconds

    logging.info(f"Sleeping for {seconds} s until next post")
    sleep(seconds)


if __name__ == '__main__':
    assert CONSUMER_KEY is not None, "CONSUMER_KEY is not set"
    assert CONSUMER_SECRET is not None, "CONSUMER_SECRET is not set"
    assert ACCESS_TOKEN is not None, "ACCESS_TOKEN is not set"
    assert ACCESS_SECRET is not None, "ACCESS_SECRET is not set"
    assert UKBB_URI_PREFIX is not None, "UKBB_URI_PREFIX is not set"
    main()
