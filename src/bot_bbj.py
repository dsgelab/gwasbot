import json
import logging

import numpy as np

from gwas_poster import GWASPoster
from utils import check_done_pheno


class BBJPoster(GWASPoster):

    def __init__(self, save_file, failure_file, manifest_file, gwas_dir):
        super().__init__(save_file, failure_file)

        # Load manifest with pheno name, description, filename
        with open(manifest_file) as f:
            manifest = json.load(f)

            # Reformat manifest with a phenocode index
            self.manifest = {}
            for record in manifest:
                pheno = record.pop("phenocode")
                self.manifest[pheno] = record

        self.gwas_dir = gwas_dir

        # Build list of phenos to post
        posted = check_done_pheno(save_file, failure_file)
        phenos = [phenocode for phenocode in self.manifest.keys()]
        self.to_post = list(set(phenos) - set(posted))
        np.random.shuffle(self.to_post)

    def build_post(self, pheno):
        logging.info("Building twitter post for BioBank Japan")
        info = self.manifest[pheno]

        post_data = {
            "longname": info["description"],
            "pheweb_link": info["pheweb_link"],
        }

        return post_data

    def format_post(self, post):
        logging.info("Formatting twitter post")
        text = f"""{post['longname']}

🇯🇵 Download and PheWeb: {post['pheweb_link']}
"""
        return text.strip()

    def gwas_img(self, pheno):
        file = self.manifest[pheno]["file"]
        return self.gwas_dir / file
