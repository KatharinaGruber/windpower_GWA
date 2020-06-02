import glob
import urllib.request
from paths_nz import nz_path
"""
Download wind power generation data New Zealand Aug 1997 - Dez 2019
URL: https://www.emi.ea.govt.nz/Wholesale/Datasets/Generation/Generation_MD/
Source: Electricity Authority

"""

url_root = 'https://www.emi.ea.govt.nz/Wholesale/Datasets/Generation/Generation_MD/'
destination = nz_path + '/generation/'
file_extension = '_Generation_MD.csv'

for year in range(1997, 2020):
    for month in range(1, 13):
        if (year != 1997)|(month>7):
            url = url_root + str(year) + str(('%02d') % month) + file_extension

            dest = destination + str(year) + str(('%02d') % month) + file_extension
            if dest not in glob.glob(destination + '*.csv'):
                try:
                    urllib.request.urlretrieve(url, dest)
                except Exception:
                    print('File ' + str(year) + str(('%02d') % month) + file_extension + ' does not exist')
                    pass
