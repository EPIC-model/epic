from tools.epic_config import *
from datetime import datetime
from dateutil.tz import tzlocal
import time

def write_nc_info(ncfile, file_type):
    ncfile.setncattr('EPIC_version', package_version)
    ncfile.setncattr('file_type', file_type)
    ncfile.setncattr('Conventions', cf_version)
    write_nc_timestamp(ncfile)

def write_nc_timestamp(ncfile):
    # 9 March 2022
    # https://stackoverflow.com/questions/32490629/getting-todays-date-in-yyyy-mm-dd-in-python
    # https://docs.python.org/3/library/datetime.html
    # https://stackoverflow.com/questions/415511/how-to-get-the-current-time-in-python
    # https://stackoverflow.com/questions/35057968/get-system-local-timezone-in-python
    ncfile.setncattr('creation_date', datetime.today().strftime('%Y/%m/%d'))
    ncfile.setncattr('creation_time', datetime.now().strftime('%H:%M:%S'))
    ncfile.setncattr('creation_zone', "UTC" + datetime.now(tzlocal()).strftime('%z'))


def write_nc_parameters(ncfile, params = {}):
    pgrp = ncfile.createGroup("parameters")
    for key, val in params.items():
        pgrp.setncattr(key, val)
