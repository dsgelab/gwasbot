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

        # Get phenotype info from manifest
        pheno_info = self.manifest.loc[
            self.manifest.phenotype == pheno,
            ["pheno_desc", "dropbox", "ukbb"]
        ].iloc[0]  # go from a pd.Series to a scalar value

        # Genetic correlations
        corr_pheno = pheno.rstrip("_irnt")  # this portal doesn't use _irnt names
        ukbbrg = f"https://ukbb-rg.hail.is/rg_summary_{corr_pheno}.html"

        # Heritability
        heritability = h2_ci(pheno, self.h2)

        # Heritability link
        ukbbh2 = f"https://nealelab.github.io/UKBB_ldsc/h2_summary_{pheno}.html"

        # Top SNP
        snp = top_snp(pheno_info.dropbox)
        gnomad_snp = snp.replace(":", "-")
        gnomad = f"https://gnomad.broadinstitute.org/variant/{gnomad_snp}"

        post = {
            "pheno": pheno_info["pheno_desc"],
            "ukbb_link": pheno_info["ukbb"],
            "download": pheno_info["dropbox"],
            "ukbbrg_link": ukbbrg,
            "ukbbh2_link": ukbbh2,
            "heritability": heritability,
            "top_snp": snp,
            "gnomad_link": gnomad,
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

        # Don't show description link if it is "nan"
        if pd.isna(post['ukbb_link']):
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
{post['top_snp']} {post['gnomad_link']}"""

        return text

    def gwas_img(self, pheno):
        """Local path to the GWAS plot"""
        return self.gwas_dir / (pheno + self.gwas_file_suffix)


def load_phenos(gwas_dir):
    """Get all the phenocodes"""
    phenos = [plot.stem.rstrip("_MF") for plot in gwas_dir.iterdir()]
    return phenos


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

    # Merge with the traits of interest
    keep = manifest["phenotype"].isin(phenos)
    manifest = manifest[keep]

    # Remove duplicated phenos
    dups = manifest["phenotype"].duplicated(keep="last")
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


def top_snp(dropbox_url):
    """Find the top SNP after downloading the phenotype data"""
    logging.info("Finding top SNP")

    # Dropbox will redirect to the actual file content for user agents
    # such as curl and wget. If we don't specicify a user agent then
    # we get an HTML response and would need to parse it to find the
    # relevant download link.
    headers = {"user-agent": "curl/7.58"}
    resp = requests.get(dropbox_url, headers=headers)
    resp.raise_for_status()
    buffer = BytesIO(resp.content)

    # Load data
    df = pd.read_csv(
        buffer,
        usecols=["variant", "pval", "beta", "se", "low_confidence_variant"],
        dialect=excel_tab,
        compression="gzip"
    )

    # Find top SNP
    df = df[df.low_confidence_variant==False]
    df['zstat'] = abs(df.beta/df.se)
    snp = df.sort_values(by="zstat",ascending=False).iloc[0].variant

    return snp
