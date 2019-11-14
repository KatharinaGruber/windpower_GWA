import argparse
import datetime
import glob
import math
import numpy as np
import os
import pandas as pd
import rasterio
import seaborn as sns
import statsmodels.api as sm
import time
import xarray as xr

import sys
sys.path.append('../')

from functools import reduce
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d

from utils import power_curve
from utils import windpower_simulation_era5
from utils import windpower_simulation_merra2
from utils import windpower_simulation_era5_large
from utils import windpower_simulation_merra2_large


from dask.diagnostics import ProgressBar
ProgressBar().register()

from paths_usa import *

os.chdir(script_path)
'''
os.system('python sim_oneyear.py -year 2018')
os.system('python sim_oneyear.py -year 2017')
os.system('python sim_oneyear.py -year 2011')
os.system('python sim_oneyear.py -year 2010')
'''
os.system('python cor_oneyear.py -year 2018')
os.system('python cor_oneyear.py -year 2017')
os.system('python cor_oneyear.py -year 2011')
os.system('python cor_oneyear.py -year 2010')


