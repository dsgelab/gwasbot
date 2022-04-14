import json
import logging

import numpy as np

from gwas_poster import GWASPoster
from utils import check_done_pheno


class METSIMPoster(GWASPoster):

    def __init__(self, save_file, failure_file, manifest_file, gwas_dir):
        super().__init__(save_file, failure_file)

        # Load manifest with pheno name, description, filename
        with open(manifest_file) as f:
            manifest = json.load(f)

            # Reformat manifest with a biochemical name index
            self.manifest = {}
            for record in manifest:
                bio_name = record.pop("BIOCHEMICAL_NAME")
                self.manifest[bio_name] = record

        self.gwas_dir = gwas_dir

        # Build list of biochemical names to post
        posted = check_done_pheno(save_file, failure_file)
        bio_names = [mm for mm in self.manifest.keys()]
        self.to_post = list(set(bio_names) - set(posted))
        np.random.shuffle(self.to_post)

    def build_post(self, bio_name):
        logging.info("Building twitter post for METSIM")
        info = self.manifest[bio_name]

        post_data = {
            "longname": bio_name,
            "pheweb_link": info["pheweb_link"],
        }

        return post_data

    def format_post(self, post):
        logging.info("Formatting twitter post")
        text = f"""{post['longname']}

🇫🇮 Download and PheWeb: {post['pheweb_link']}
"""
        return text.strip()

    def gwas_img(self, bio_name):
        file = self.manifest[bio_name]["file"]
        return self.gwas_dir / file
