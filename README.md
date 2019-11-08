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
export CONSUMER_KEY=... CONSUMER_SECRET=... ACCESS_TOKEN=... ACCESS_SECRET=...  # credentials from Twitter
export DATA_PATH=...  # path to the data directory containing the input files
export URI_PREFIX_UKBB=...  # Google Storage path to GWAS files, excluding filename, starts with gs://
source venv/bin/activate
python main.py
```

For long-running scenario, one could do the previous commands in a `tmux` session:

```bash
tmux

# previous commands
# ...

# exit with Ctrl-b d
```

Then one can logout of the server running the GWASbot without it stopping.
To reconnect to the session: `tmux attach`.
