import os.path as op
import pathlib

NUM_PROCESSES = 8

# used for downloading, calucation of time series etc
YEARS = range(1995, 2019)
MONTHS = range(1, 13)

DISTANCE_FACTORS = 2, 3, 4, 6

LOG_FILE = op.join(op.dirname(__file__), '..', 'data', 'logfile.log')

INTERIM_DIR = pathlib.Path(__file__).parent.parent / 'data' / 'interim'

EXTERNAL_DIR = pathlib.Path(__file__).parent.parent / 'data' / 'external'

FIGURES_DIR = pathlib.Path(__file__).parent.parent / 'figures'

FIGSIZE = (12, 7.5)