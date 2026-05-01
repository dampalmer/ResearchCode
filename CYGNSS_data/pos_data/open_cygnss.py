import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt

# Open the netCDF file
ds = nc.Dataset('cygnss_extracted_data.nc')

print(ds.variables)
