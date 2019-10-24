# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 09:00:37 2019

@author: KatharinaG

adapted from source:
https://github.com/inwe-boku/wind_repower_usa/blob/master/scripts/download_wind_era5.py
"""

import os
import os.path as op
import logging

import cdsapi

from logging_config import setup_logging

from paths_usa import era_path

DOWNLOAD_DIR = era_path

YEARS = range(2000, 2019)
MONTHS = list(range(1, 13))

setup_logging(op.join(DOWNLOAD_DIR, 'download.log'))


def main():
    # API documentation for downloading a subset:
    # https://confluence.ecmwf.int/display/CKB/Global+data%3A+Download+data+from+ECMWF+for+a+particular+area+and+resolution
    # https://retostauffer.org/code/Download-ERA5/

    north = 68
    south = 12
    east = -173
    west = -64

    # Format for downloading ERA5: North/West/South/East
    bounding_box = "{}/{}/{}/{}".format(north, west, south, east)

    # old data downloaded with:  "68/-172/16/-61"

    logging.info("Downloading bounding_box=%s for years=%s and months=%s",
                 bounding_box, YEARS, MONTHS)

    c = cdsapi.Client()

    for year in YEARS:
        for month in MONTHS:
            filename = op.join(DOWNLOAD_DIR,
                               f'era5_wind_USA_{year}{month:02d}.nc')

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

    if total != size:
        raise Exception("Download failed: downloaded %s byte(s) out of %s" % (total, size))

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

    # add md5 hashes if you know that the version can be safely monkey patched:
    expected_md5s = ('1e93fc6bbd1cc825a10845f47c59835f',
                     '3a96092ca74d684a26615958e8297c4f')

    cdsapi_md5 = hashlib.md5(cdsapi_src).hexdigest()
    if cdsapi_md5 not in expected_md5s:
        logging.warning("Could not patch cdsapi, unexpected md5sum of "
                        f"api.py: {cdsapi_md5}")
        return

    cdsapi.api.Result._download = _cdsapi_download_with_timeout


if __name__ == '__main__':
    patch_cdsapi()
main()