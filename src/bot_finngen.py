import json
import logging
from csv import excel_tab
from io import BytesIO

import requests
import numpy as np
import pandas as pd

from gwas_poster import GWASPoster
from utils import check_done_pheno


class FGPoster(GWASPoster):

    def __init__(
            self,
            save_file,
            failure_file,
            metadata_file,
            gwas_dir):
        """Setup and load the data for FinnGen"""
        super().__init__(save_file, failure_file)

        # Load metadata file with phenos information
        with open(metadata_file) as f:
            self.metadata = json.load(f)

        # Set data path
        self.gwas_dir = gwas_dir

        # Get the phenotypes to post
        posted = check_done_pheno(save_file, failure_file)
        phenos = [record["phenocode"] for record in self.metadata]

        self.to_post = list(set(phenos) - set(posted))
        np.random.shuffle(self.to_post)

    def build_post(self, pheno):
        """Gather the data that will be displayed in the tweeter post"""
        logging.info("Building twitter post for FinnGen")
        meta = self.get_meta(pheno)

        top_variant = find_top_hit(meta["download"])
        variant_link = f"https://r2.finngen.fi/variant/{top_variant}"

        post_data = {
            'pheno_longname': meta["description"],
            'download_link': meta["download"],
            'pheweb_pheno_link': meta["pheweb_link"],
            'top_variant': top_variant,
            'pheweb_variant_link': variant_link,
        }

        return post_data

    def format_post(self, post):
        """Textual representation of the twitter post"""
        logging.info("Formatting twitter post")
        limit = 150
        if len(post['pheno_longname']) > limit:
            pheno = post['pheno'][:limit - 2] + "…"
        else:
            pheno = post['pheno_longname']

        text = f"""{pheno}

🇫🇮 PheWeb {post['pheweb_pheno_link']}

⬇️ Download {post['download_link']}

📊 GWAS top hit
{post['top_variant']} {post['pheweb_variant_link']}
""".strip()

        return text

    def gwas_img(self, pheno):
        """Local path to the GWAS plot"""
        meta = self.get_meta(pheno)
        return self.gwas_dir / meta["file"]

    def get_meta(self, pheno):
        """Get meta data on a specific pheno"""
        meta = list(filter(lambda r: r["phenocode"] == pheno, self.metadata))
        assert len(meta) == 1
        meta = meta[0]
        return meta


def find_top_hit(url):
    """Download a GWAS data file to an in-memory buffer"""
    logging.info("Finding GWAS top hit")

    # Download the GWAS data
    logging.info(f"Downloading GWAS data at {url}")
    resp = requests.get(url)
    buffer = BytesIO(resp.content)

    # Put data into a DataFrame to get variant ID and p-value
    logging.info("Finding top hit with Pandas")
    df = pd.read_csv(
        buffer,
        usecols=["#chrom", "pos", "ref", "alt", "pval"],
        dialect=excel_tab,
        compression="gzip")

    # Find top-variant using "pval"
    variant = df.iloc[df.pval.idxmin()]

    # Convert variant info to variant CPRA
    cpra = f"{variant['#chrom']}:{variant.pos}-{variant.ref}-{variant.alt}"

    return cpra
