# data i/o
import os

import numpy as np
import pandas as pd


import xarray as xr
import datetime
import matplotlib.pyplot as plt

cor_2y = xr.open_dataset('/data/users/kgruber/results/USA/correlations/USA_2y.nc').c

cor_1y1 = xr.open_dataset('/data/users/kgruber/results/USA/correlations/USA_1y1.nc').c
cor_1y2 = xr.open_dataset('/data/users/kgruber/results/USA/correlations/USA_1y2.nc').c

def plot1():
    plt.imshow(cor_2y,cmap='gist_rainbow')
    plt.colorbar()
    plt.savefig('/data/users/kgruber/results/USA/correlations/USA2y.png',dpi=300)
    

def plot2():
    plt.imshow(cor_1y1,cmap='gist_rainbow')
    plt.colorbar()
    plt.savefig('/data/users/kgruber/results/USA/correlations/USA1y1.png',dpi=300)

def plot3():
    plt.imshow(cor_1y2,cmap='gist_rainbow')
    plt.colorbar()
    plt.savefig('/data/users/kgruber/results/USA/correlations/USA1y2.png',dpi=300)

def plot4():
    plt.imshow(cor_2y-cor_1y1,cmap='nipy_spectral')
    plt.colorbar()
    plt.savefig('/data/users/kgruber/results/USA/correlations/USA2y-1y1.png',dpi=300)

def plot5():
    plt.imshow(cor_2y - cor_1y2,cmap='nipy_spectral')
    plt.colorbar()
    plt.savefig('/data/users/kgruber/results/USA/correlations/USA2y-1y2.png',dpi=300)
    return

def plot6():
    plt.imshow(cor_1y1 - cor_1y2,cmap='nipy_spectral')
    plt.colorbar()
    plt.savefig('/data/users/kgruber/results/USA/correlations/USA1y1-1y2.png',dpi=300)
    
if(not os.path.isfile('/data/users/kgruber/results/USA/correlations/USA2y.png')):
    print('plot1')
    plot1()
    exit()
    
if(not os.path.isfile('/data/users/kgruber/results/USA/correlations/USA1y1.png')):
    print('plot2')
    plot2()
    exit()
    
if(not os.path.isfile('/data/users/kgruber/results/USA/correlations/USA1y2.png')):
    print('plot3')
    plot3()
    exit()

if(not os.path.isfile('/data/users/kgruber/results/USA/correlations/USA2y-1y1.png')):
    print('plot4')
    plot4()
    exit()
    
if(not os.path.isfile('/data/users/kgruber/results/USA/correlations/USA2y-1y2.png')):
    print('plot5')
    plot5()
    
if(not os.path.isfile('/data/users/kgruber/results/USA/correlations/USA1y1-1y2.png')):
    print('plot6')
    plot6()
    exit()