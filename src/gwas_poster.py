import logging
from os import getenv

import tweepy


class GWASPoster(object):

    def __init__(self, save_file, failure_file):
        self.save_file = save_file
        self.failure_file = failure_file

        self.consumer_key = getenv("CONSUMER_KEY")
        self.consumer_secret = getenv("CONSUMER_SECRET")
        self.access_token = getenv("ACCESS_TOKEN")
        self.access_secret = getenv("ACCESS_SECRET")

        assert self.consumer_key is not None, "CONSUMER_KEY is not set"
        assert self.consumer_secret is not None, "CONSUMER_SECRET is not set"
        assert self.access_token is not None, "ACCESS_TOKEN is not set"
        assert self.access_secret is not None, "ACCESS_SECRET is not set"

    def get_pheno(self):
        """Get the next phenotype to post"""
        if len(self.to_post) > 0:
            # NOTE self.to_post needs to be a list defined in sub-classes
            return self.to_post[0]
        else:
            raise ValueError("no more phenotype to pick from")

    def build_post(self, pheno):
        raise NotImplementedError

    def format_post(self, post):
        raise NotImplementedError

    def gwas_img(self, pheno):
        raise NotImplementedError

    def tweet(self, pheno):
        logging.info(f"Posting {pheno} to twitter")

        # Gather post info
        post = self.build_post(pheno)
        text = self.format_post(post)
        img_filepath = self.gwas_img(pheno)

        # Connecting to Twitter API
        auth = tweepy.OAuthHandler(self.consumer_key, self.consumer_secret)
        auth.set_access_token(self.access_token, self.access_secret)
        api = tweepy.API(auth)

        # Upload GWAS image
        gwas = api.media_upload(str(img_filepath))

        # Post to twitter
        api.update_status(text, media_ids=[gwas.media_id_string])

    def mark_posted(self, pheno):
        """Record that the given pheno was posted to twitter"""
        with open(self.save_file, "a") as f:
            f.write(f"{pheno}\n")
        self.to_post.remove(pheno)

    def mark_failure(self, pheno):
        """Record that a given pheno has trouble being posted to twitter"""
        with open(self.failure_file, "a") as f:
            f.write(f"{pheno}\n")
        self.to_post.remove(pheno)
