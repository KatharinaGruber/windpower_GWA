# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 11:37:01 2019

@author: KatharinaG
"""

import datetime
import glob
import numpy as np
import os
from calendar import monthrange
import subprocess

max_number_retries = 10

def download_month(date, lon1, lat1, lon2, lat2, var, user, password, outpath):
    '''
    function for downloading MERRA-2 data per month

    input parameters:

    date	  yearmonth to download in format yyyymm as str
    lon1,lat1,lon2,lat2    geographical extent in degrees as float
    var    variables to download, either str or list of str
    user,password    credentials as str
    outpath     path to save final files as str
    '''

    try:
        os.chdir(outpath)
    except Exception:
        print('Directory does not exist')
        exit()

    y = date[:4]
    m = date[4:]
    if int(y) not in range(1980, datetime.datetime.now().year + 1):
        print('invalid year')
        exit()
    if (int(y) == datetime.datetime.now().year) & (int(m) > datetime.datetime.now().month):
        print('invalid date')
        exit()
    if (int(m) > 12) | (int(m) < 1):
        print('invalid month')
        exit()
    dates = np.array(range(1, monthrange(int(y), int(m))[1] + 1)) + int(date) * 100
    dates = dates.astype('<U8')

    x1 = str(round((lat1+90)/0.5))
    y1 = str(round((lon1+180)/0.625))
    x2 = str(round((lat2+90)/0.5))
    y2 = str(round((lon2+180)/0.625))

    xy = '[0:23]['+x1+':'+x2+']['+y1+':'+y2+']'

    if type(var) == list:
        vs = [v + xy for v in var]
        vstr = vs[0]
        for v in vs[1:]:
            vstr = vstr + ',' + v
    else:
        vstr = var + xy

    if int(y) in range(1980, 1992):
        ds = '100'
    elif int(y) in range(1992, 2001):
        ds = '200'
    elif int(y) in range(2001, 2011):
        ds = '300'
    elif int(y) in range(2011, 2020):
        ds = '400'
    else:
        print('I do not know the dataset number yet...')
        exit()
    
    #fl = 0
    #while fl < len(dates):
    for date in dates:
        if outpath + '/MERRA2_' + ds + '.tavg1_2d_slv_Nx.' + date + '.nc4.nc' not in glob.glob(outpath + '/*'):
            print('downloading ' + date)
            for i in range(max_number_retries):
                r = os.system('wget -q --content-disposition --user ' + user + ' --password ' +
                #subprocess.check_call('wget -q --content-disposition --user ' + user + ' --password ' +
                          password + ' https://goldsmr4.gesdisc.eosdis.nasa.gov/' +
                          'opendap/MERRA2/M2T1NXSLV.5.12.4/' + y + '/' + m + '/MERRA2_' + ds +
                          '.tavg1_2d_slv_Nx.' + date + '.nc4.nc4?' + vstr + ',time,lat[' + x1 +
                          ':' + x2 + '],lon[' + y1 + ':' + y2 + ']')
                if r==0:
                    break
                else:
                    print('retry ' + str(i) + ' downloading ' + date)
                if (r!=0)&(i==(max_number_retries-1)):
                    print('download ' + date + ' failed!')
        #fl = len(glob.glob(outpath + '/MERRA2_' + ds + '.tavg1_2d_slv_Nx.' + y + m + '??.nc4.nc'))
