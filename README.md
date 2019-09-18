# GWASbot


## Setup

Create a python virtual environment and install the dependencies:
```bash
python3 -m venv venv
source venv/bin/activate
pip install -U pip
pip install -r requirements.txt
```

## Running

```bash
export GOOGLE_APPLICATION_CREDENTIALS="/path/to/keyfile.json"  # GCP service account credentials
export CONSUMER_KEY=... CONSUMER_SECRET=... BOT_ACCESS_TOKEN=... BOT_ACCESS_SECRET=...  # credentials from Twitter
source venv/bin/activate
export UKBB_URI_PREFIX=...  # path to GWAS files, excluding filename, starts with gs://...
python bot.py
```
