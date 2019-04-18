#!/usr/bin/env python
# coding: utf-8
""" Download ERA5 South Africa

adapted from
https://github.com/inwe-boku/wind_repower_usa/blob/master/scripts/download_wind_era5.py
script for downloading era5 single levels data with cds api
get account at:
https://cds.climate.copernicus.eu/user/register?destination=%2F%23!%2Fhome
check license agreement
install csdapi and key: https://cds.climate.copernicus.eu/api-how-to
(to create .cdsapirc file create text file and rename to .cdsapirc.)
if using conda, make sure to install it via conda prompt
see data you can download:
https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form
data are availabe from 2000 until now
for variables available see link above
"""

import os
import os.path as op
import logging
import cdsapi

from logging_config import setup_logging


# folder where data shall be stored
DOWNLOAD_DIR = os.getcwd()

# and also find the bounding box
lon1 = 16
lon2 = 33
lat1 = -35
lat2 = -22
bounding_box = [lat2, lon1, lat1, lon2]

# acronym of the country
country = 'ZAF'

YEARS = range(2014, 2018)
MONTHS = list(range(1, 13))

setup_logging(op.join(DOWNLOAD_DIR, 'download.log'))


def main():
    # API documentation for downloading a subset:
    # https://confluence.ecmwf.int/display/CKB/Global+data%3A+Download+data+from+ECMWF+for+a+particular+area+and+resolution
    # https://retostauffer.org/code/Download-ERA5/

    logging.info("Downloading bounding_box=%s for years=%s and months=%s",
                 bounding_box, YEARS, MONTHS)

    c = cdsapi.Client()

    for year in YEARS:
        for month in MONTHS:
            print(str(year)+str(month))
            filename = op.join(DOWNLOAD_DIR,
                               f'era5_wind_{country}_{year}{month:02d}.nc')

            if op.exists(filename):
                logging.info(f"Skipping {filename}, already exists!")
                continue

            logging.info(f"Starting download of {filename}...")

            for i in range(5):
                try:
                    c.retrieve(
                        'reanalysis-era5-single-levels',
                        {
                            'product_type': 'reanalysis',
                            'format': 'netcdf',
                            'variable': [
                                '100m_u_component_of_wind',
                                '100m_v_component_of_wind',
                                '10m_u_component_of_wind',
                                '10m_v_component_of_wind'
                            ],
                            'year': f'{year}',
                            'month': [
                                f'{month:02d}'
                            ],
                            'area': bounding_box,
                            'day': [f"{day:02d}" for day in range(1, 32)],
                            'time': [f"{hour:02d}:00" for hour in range(24)],
                        },
                        f"{filename}.part"
                    )
                except Exception as e:
                    logging.warning("Download failed: %s", e)
                else:
                    logging.info(f"Download of {filename} successful!")
                    os.rename(f"{filename}.part", filename)
                    break
            else:
                logging.warning("Download failed permanently!")


def _cdsapi_download_with_timeout(self, url, size, target):
    """Copied from cdsapi.api, see below."""
    from cdsapi.api import bytes_to_string, time, requests

    if target is None:
        target = url.split('/')[-1]

    self.info("Downloading %s to %s (%s)", url, target, bytes_to_string(size))
    start = time.time()

    r = self.robust(requests.get)(url, stream=True, verify=self.verify,
                                  timeout=20)
    try:
        r.raise_for_status()

        total = 0
        with open(target, 'wb') as f:
            for chunk in r.iter_content(chunk_size=1024):
                if chunk:
                    f.write(chunk)
                    total += len(chunk)
    finally:
        r.close()

    assert total == size

    elapsed = time.time() - start
    if elapsed:
        self.info("Download rate %s/s", bytes_to_string(size / elapsed))

    return target


def patch_cdsapi():
    """Monkey patch the cdsapi and add 20s timeout for hanging downloads.
    See also:
    https://jira.ecmwf.int/servicedesk/customer/portal/1/CUS-7104
    > I am using the cdsapi package to download ERA5 data. I am having
    > troubles with frozen download connections. I don't really know why,
    > maybe this is due to my internet connection.
    > I have not experienced that such a frozen connection recovers, so I
    > guess it might make sense to add a timeout. This is also what the
    > documentation suggests:
    > > Most requests to external servers should have a timeout attached, in
    > > case the server is not responding in a timely manner. By default,
    > > requests do not time out unless a timeout value is set explicitly.
    > > Without a timeout, your code may hang for minutes or more.
    > http://docs.python-requests.org/en/master/user/advanced/#timeouts
    > https://stackoverflow.com/questions/45267003/python-requests-hanging-freezing
    > This leads to an exception of the following form after 20 seconds of no
    > response:
    > HTTPConnectionPool(host='136.156.132.105', port=80): Read timed out.
    > I am not sure if it is possible to recover from this situation in a way
    > without re-downloading everything. But even a simple retry by catching
    > this exception makes the download easier.
    """

    # to make sure that cdsapi hasn't changed, let's patch only
    import inspect
    import hashlib
    cdsapi_src = inspect.getsource(cdsapi.api).encode('utf-8')
    expected_md5 = '1e93fc6bbd1cc825a10845f47c59835f'
    cdsapi_md5 = hashlib.md5(cdsapi_src).hexdigest()
    if cdsapi_md5 != expected_md5:
        raise RuntimeError("Could not patch cdsapi, md5sum of api.py"
                           f"has changed. {cdsapi_md5}")

    cdsapi.api.Result._download = _cdsapi_download_with_timeout


if __name__ == '__main__':
    patch_cdsapi()
main()
