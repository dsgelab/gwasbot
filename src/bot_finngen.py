import logging
from csv import excel_tab

import numpy as np
import pandas as pd
from google.api_core.exceptions import NotFound

from gwas_poster import GWASPoster
from utils import check_posted
from utils import gcloud_download


class FGPoster(GWASPoster):

    def __init__(
            self,
            endpoints_file,
            save_file,
            gwas_dir,
            gwas_file_suffix,
            gcloud_uri_prefix):
        """Setup and load the data for FinnGen"""
        super().__init__(save_file)

        # Set data paths
        self.gwas_dir = gwas_dir
        self.gwas_file_suffix = gwas_file_suffix
        self.gcloud_uri_prefix = gcloud_uri_prefix

        # Get the phenotypes to post
        posted = check_posted(save_file)
        phenos = get_phenos(gwas_dir, "_MF.png")

        self.to_post = list(set(phenos) - set(posted))
        np.random.shuffle(self.to_post)

        # Load the endpoint definitions
        self.endpoints = get_endpoint_definitions(endpoints_file)


    def build_post(self, pheno):
        """Gather the data that will be displayed in the tweeter post"""
        download_link = ""  # TODO

        pheweb_pheno_link = f"http://r2.finngen.fi/pheno/{pheno}"

        try:
            top_variant = find_variant_top_hit(pheno, self.gcloud_uri_prefix)
        except NotFound:
            logging.warning(f"GWAS data not found for pheno: {pheno}")
            top_variant = None
            pheweb_variant_link = None
        else:
            pheweb_variant_link = f"http://r2.finngen.fi/variant/{top_variant}"

        post_data = {
            'pheno_longname': get_pheno_longname(pheno, self.endpoints),
            'download_link': download_link,
            'pheweb_pheno_link': pheweb_pheno_link,
            'top_variant': top_variant,
            'pheweb_variant_link': pheweb_variant_link,
        }

        return post_data

    def format_post(self, post):
        """Textual representation of the twitter post"""
        if post['top_variant'] is not None:
            top_hit = f"\n📊 GWAS top hit\n{post['top_variant']} {post['pheweb_variant_link']}\n"
        else:
            top_hit = ""

        text = f"""{post['pheno_longname']}

⬇️ Download {post['download_link']}

🧬 PheWeb {post['pheweb_pheno_link']}
{top_hit}
""".strip()
        return text

    def gwas_img(self, pheno):
        """Local path to the GWAS plot"""
        return self.gwas_dir / (pheno + self.gwas_file_suffix)


def get_phenos(gwas_dir, file_suffix):
    """Get the list of phenotypes to post"""
    files = filter(lambda f: f.name.endswith(file_suffix), gwas_dir.iterdir())
    phenos = map(lambda f: f.name.rstrip(file_suffix), files)
    phenos = list(phenos)
    return phenos


def get_endpoint_definitions(endpoints_file):
    """Load all the endpoints from the FinnGen endpoints file"""
    df = pd.read_csv(
        endpoints_file,
        usecols=["NAME", "LONGNAME"],
        dialect=excel_tab,
    ).drop(0)  # remove the comment line

    return df


def get_pheno_longname(pheno, endpoints):
    """Return the phenotype longname / description"""
    longname = endpoints.loc[endpoints.NAME == pheno, "LONGNAME"].iloc[0]
    return longname


def find_variant_top_hit(pheno, uri_prefix):
    """Given a phenotype name, return the variant with the lowest p-value"""
    # Download GWAS data from Google Cloud
    gwas = download_gwas(pheno, uri_prefix)

    # Load file using pandas
    df = pd.read_csv(
        gwas,
        usecols=["#chrom", "pos", "ref", "alt", "pval"],
        dialect=excel_tab,
        compression="gzip"
    )

    # Find top-variant using "pval"
    variant = df.iloc[df.pval.idxmin()]

    # Convert variant info to variant CPRA
    cpra = f"{variant['#chrom']}:{variant.pos}-{variant.ref}-{variant.alt}"

    return cpra


def download_gwas(pheno, gcloud_uri_prefix):
    """Download GWAS from Google Cloud Storage into memory"""
    blob = gcloud_uri_prefix + pheno + ".gz"
    buffer = gcloud_download(blob)
    return buffer
