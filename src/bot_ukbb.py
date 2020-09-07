import logging
from csv import excel_tab
from io import BytesIO

import numpy as np
import pandas as pd
import requests

from gwas_poster import GWASPoster
from utils import check_done_pheno


class UKBBPoster(GWASPoster):

    def __init__(
            self,
            h2_file,
            manifest_file,
            save_file,
            failure_file,
            gwas_dir,
            gwas_file_suffix):
        """Setup and load the data for UKBB"""
        super().__init__(save_file, failure_file)
        # Set data paths
        self.gwas_dir = gwas_dir
        self.gwas_file_suffix = gwas_file_suffix

        # Get the phenotypes to post
        posted = check_done_pheno(save_file, failure_file)
        phenos = load_phenos(gwas_dir)
        self.to_post = list(set(phenos) - set(posted))
        np.random.shuffle(self.to_post)

        self.h2 = load_h2(h2_file, phenos)
        self.manifest = load_manifest(manifest_file, phenos)

    def build_post(self, pheno):
        """Gather all the info we need for the twitter post"""
        logging.info("Building twitter post for UKBB")

        clean_pheno = pheno.replace("_irnt", "")  # for resources that don't use "_irnt" suffix

        # Get phenotype info from manifest
        pheno_info = self.manifest.loc[
            (self.manifest.phenocode == clean_pheno) | (self.manifest.phenocode_coding == clean_pheno),
            ["pheno_desc", "aws"]
        ].iloc[0]  # go from a pd.Series to a scalar value

        # Genetic correlations
        ukbbrg = f"https://ukbb-rg.hail.is/rg_summary_{clean_pheno}.html"

        # Heritability
        heritability = h2_ci(pheno, self.h2)

        # Heritability link
        ukbbh2 = f"https://nealelab.github.io/UKBB_ldsc/h2_summary_{pheno}.html"

        # Top SNP
        snp_cpra = top_snp(pheno_info.aws)

        # Look up SNP on h38, as UKBB uses h37 but OpenTargets uses h38
        gnomad_cpra_h37 = snp_cpra.replace(":", "-")
        gnomad_cpra_h38 = h37_to_h38(gnomad_cpra_h37)

        opentargets_cpra = gnomad_cpra_h38.replace("-", "_")
        opentargets = f"https://genetics.opentargets.org/variant/{opentargets_cpra}"

        post = {
            "pheno": pheno_info["pheno_desc"],
            "ukbb_link": get_ukbb_link(clean_pheno),
            "download": pheno_info["aws"],
            "ukbbrg_link": ukbbrg,
            "ukbbh2_link": ukbbh2,
            "heritability": heritability,
            "top_snp": snp_cpra,
            "opentargets_link": opentargets,
        }
        return post

    def format_post(self, post):
        """Textual representation of the post"""
        logging.debug("Formatting twitter post")
        # Shorten the phenocode name if too long for the tweet
        limit = 35  # based on trial and error
        if len(post['pheno']) > limit:
            pheno = post['pheno'][:limit - 2] + "…"
        else:
            pheno = post['pheno']

        # Don't show description link if it's not on UKBB portal
        if post['ukbb_link'] is None:
            desc = ""
        else:
            desc = f"\n🇬🇧 Description {post['ukbb_link']}\n"

        text = f"""{pheno}
{desc}
⬇️ Download {post['download']}

🧬 Genetic correlations {post['ukbbrg_link']}

👪 Heritability
{post['heritability']['h2']:.2f} [{post['heritability']['ci_min']:.2f}, {post['heritability']['ci_max']:.2f}] {post['ukbbh2_link']}

📊 GWAS top hit
{post['top_snp']} {post['opentargets_link']}"""

        return text

    def gwas_img(self, pheno):
        """Local path to the GWAS plot"""
        return self.gwas_dir / (pheno + self.gwas_file_suffix)


def load_phenos(gwas_dir):
    """Get all the phenocodes"""
    phenos = [plot.stem.replace("_trans_MF", "") for plot in gwas_dir.iterdir()]
    return phenos


def get_ukbb_link(pheno):
    """Try to guess the UKBB portal link for a phenotype.

    Not all phenotypes are in the UKBB portal, so this function checks
    the UKBB portal first.

    This is made harder by the UKBB portal not returning a 404 for
    phenotypes not found.
    """
    link = f"https://biobank.ctsu.ox.ac.uk/crystal/field.cgi?id={pheno}"
    error_msg = "Sorry, an internal error prevents field listing"
    resp = requests.get(link)

    try:
        resp.raise_for_status()
    except requests.exceptions.HTTPError:
        link = None
    else:
        if error_msg in resp.text:
            link = None

    return link


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
            "phenocode",
            "coding",
            "description",
            "aws_link",
            "pheno_sex",
        ],
        dtype={
            "phenocode": "object",
            "coding": "object",
        }
    )

    # Create phenotype column as phenocode_coding
    # Some GWAS filename are "phenocode", some are "phenocode_coding".
    # Need to look for both.
    has_coding = manifest.coding.notna()
    manifest["phenocode_coding"] = ""
    manifest.loc[has_coding, "phenocode_coding"] = manifest.loc[has_coding, :].phenocode.map(str) + "_" + manifest.loc[has_coding, :].coding.map(str)
    manifest.loc[~ has_coding, "phenocode_coding"] = manifest.loc[~ has_coding, "phenocode"]

    manifest.rename(
        columns={
            "description": "pheno_desc",
            "aws_link": "aws",
            "pheno_sex": "sex",
        },
        inplace=True
    )

    # Remove NaN phenotypes
    nan_phenos = manifest["phenocode"].isna()
    manifest = manifest[~ nan_phenos]

    # Keep only "both sexes"
    both_sexes = manifest["sex"] == "both_sexes"
    manifest = manifest[both_sexes]

    # Remove traits ending in _raw
    raw_phenos = manifest["phenocode"].str.endswith("_raw")
    manifest = manifest[~ raw_phenos]

    # Remove duplicated phenos
    dups = manifest["phenocode_coding"].duplicated(keep="last")
    manifest = manifest[~ dups]

    return manifest


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


def top_snp(aws_url):
    """Find the top SNP after downloading the phenotype data"""
    logging.info("Finding top SNP")

    logging.debug("Downloading GWAS from AWS")
    resp = requests.get(aws_url)
    resp.raise_for_status()
    buffer = BytesIO(resp.content)

    # Load data
    logging.debug("Loading GWAS in Pandas and finding top SNP")
    df = pd.read_csv(
        buffer,
        usecols=[
            "chr",
            "pos",
            "ref",
            "alt",
            "beta_meta",
            "se_meta",
            "low_confidence_EUR"
        ],
        dtype={
            # We specify column types to prevent "Column type mismatch".
            # We don't specify "low_confidence_EUR" here so Pandas can
            # magically assign True, False and NaN.
            "chr": "object",
            "pos": "int",
            "ref": "object",
            "alt": "object",
            "beta_meta": "float",
            "se_meta": "float",
        },
        dialect=excel_tab,
        compression="gzip"
    )

    # Combine chr, pos, ref, alt into a "variant" column
    df["variant"] = df.chr + ":" + df.pos.map(str) + ":" + df.ref + ":" + df.alt

    # Find top SNP
    df = df[df.low_confidence_EUR==False]
    df['zstat'] = abs(df.beta_meta/df.se_meta)
    snp = df.sort_values(by="zstat",ascending=False).iloc[0].variant

    return snp


def h37_to_h38(gnomad_cpra):
    """Translate a SNP ID (as CPRA) on h37 to h38"""
    gnomad_api = 'https://gnomad.broadinstitute.org/api/'

    # 1. h37 -> rsid
    payload = {
        # Set the query as raw GrahQL query because we don't need a library just for that
        'query': 'query GnomadVariant($variantId: String, $rsid: String, $datasetId: DatasetId!) { variant(variantId: $variantId, rsid: $rsid, dataset: $datasetId) {rsid} }',
        'variables': {
            'datasetId': 'gnomad_r2_1',
            'variantId': gnomad_cpra
        }
    }
    resp = requests.post(gnomad_api, json=payload)
    resp.raise_for_status()
    rsid = (
        resp.json()
        ["data"]
        ["variant"]
        ["rsid"]
    )

    # 2. rsid -> h38
    payload = {
        'query': 'query GnomadVariant($variantId: String, $rsid: String, $datasetId: DatasetId!) { variant(variantId: $variantId, rsid: $rsid, dataset: $datasetId) {variantId} }',
        'variables': {
            'datasetId': 'gnomad_r3',
            'rsid': rsid
        }
    }
    resp = requests.post(gnomad_api, json=payload)
    resp.raise_for_status()
    cpra_h38 = (
        resp.json()
        ["data"]
        ["variant"]
        ["variantId"]
    )

    return cpra_h38
