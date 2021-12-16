import xarray as xr
import argparse
import os.path

parser = argparse.ArgumentParser()

parser.add_argument('input_file_name',type=str)
args = parser.parse_args()
input_file_name=args.input_file_name

file_ext = os.path.splitext(input_file_name)[1]
file_root = os.path.splitext(input_file_name)[0]

if not (file_ext=='.hdf5'):
    raise argparse.ArgumentTypeError('Argument filename must end on ".hdf5"')

ds=xr.open_dataset(input_file_name)
ds.to_netcdf(file_root+'.nc')
