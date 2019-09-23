"""
Functions useful for both the FinnGen and UKBB bots.
"""

import logging
from datetime import date, datetime
from io import BytesIO
from os import getenv
from time import sleep

import tweepy
import pytz
from dateutil.relativedelta import relativedelta
from google.cloud import storage


logging.basicConfig(level=logging.INFO)


# Twitter API
CONSUMER_KEY = getenv("CONSUMER_KEY")
CONSUMER_SECRET = getenv("CONSUMER_SECRET")
ACCESS_TOKEN = getenv("BOT_ACCESS_TOKEN")
ACCESS_SECRET = getenv("BOT_ACCESS_SECRET")

assert CONSUMER_KEY is not None, "CONSUMER_KEY is not set"
assert CONSUMER_SECRET is not None, "CONSUMER_SECRET is not set"
assert ACCESS_TOKEN is not None, "ACCESS_TOKEN is not set"
assert ACCESS_SECRET is not None, "ACCESS_SECRET is not set"



def gcloud_download(gs_path):
    """Download a file from Google Cloud Storage into an in-memory buffer."""
    logging.info(f"Downloading {gs_path} into memory")

    buffer = BytesIO()
    client = storage.Client()
    client.download_blob_to_file(gs_path, buffer)

    # Reset file position to be re-usable later by Pandas
    buffer.seek(0)

    return buffer


def tweet(text, img_filepath):
    """Make a post to twitter"""
    logging.info("Posting to twitter")
    auth = tweepy.OAuthHandler(CONSUMER_KEY, CONSUMER_SECRET)
    auth.set_access_token(ACCESS_TOKEN, ACCESS_SECRET)
    api = tweepy.API(auth)

    # Upload GWAS image
    gwas = api.media_upload(str(img_filepath))

    # Post to twitter
    api.update_status(text, media_ids=[gwas.media_id_string])


def wait(hour, timezone):
    """Wait until tomorrow 8am East coast time"""
    tomorrow = date.today() + relativedelta(days=1)
    tz = pytz.timezone(timezone)
    next_time = datetime(tomorrow.year, tomorrow.month, tomorrow.day, hour, 0, tzinfo=tz)

    now = datetime.now(tz=tz)
    seconds = (next_time - now).seconds

    logging.info(f"Sleeping for {seconds} s until next post")
    sleep(seconds)
