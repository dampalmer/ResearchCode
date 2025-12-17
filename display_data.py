#!/usr/bin/env python3
"""
Display the structure and sample data from the extracted CYGNSS NetCDF file.
"""

import netCDF4
import numpy as np

def display_cygnss_data(filename):
    """Display the structure and sample data from the CYGNSS extracted file."""
    
    # Open the extracted data file
    nc = netCDF4.Dataset(filename)
    
    print('=== CYGNSS EXTRACTED DATA STRUCTURE ===')
    print()
    print('Variables in the file:')
    for var_name in nc.variables:
        var = nc.variables[var_name]
        print(f'  {var_name}: {var.shape} - {var.long_name if hasattr(var, "long_name") else "No description"}')
    
    print()
    print('Global attributes:')
    for attr in nc.ncattrs():
        print(f'  {attr}: {nc.getncattr(attr)}')
    
    print()
    print('=== SAMPLE DATA (first 10 records) ===')
    print()
    
    # Get sample data
    timestamps = nc.variables['ddm_timestamp_utc'][:10]
    pos_x = nc.variables['sc_pos_x'][:10]
    pos_y = nc.variables['sc_pos_y'][:10]
    pos_z = nc.variables['sc_pos_z'][:10]
    sat_ids = nc.variables['satellite_id'][:10]
    
    print('Sample | Timestamp (UTC) | Sat ID | Position X (m) | Position Y (m) | Position Z (m)')
    print('-------|-----------------|--------|----------------|----------------|----------------')
    
    for i in range(10):
        print(f'{i+1:6d} | {timestamps[i]:15.3f} | {sat_ids[i]:6d} | {pos_x[i]:14.1f} | {pos_y[i]:14.1f} | {pos_z[i]:14.1f}')
    
    print()
    print('=== DATA SUMMARY BY SATELLITE ===')
    print()
    
    # Get all satellite IDs and count samples per satellite
    all_sat_ids = nc.variables['satellite_id'][:]
    unique_sats = np.unique(all_sat_ids)
    
    for sat_id in unique_sats:
        count = np.sum(all_sat_ids == sat_id)
        print(f'Satellite CYG{sat_id:02d}: {count:,} samples')
    
    print()
    print(f'Total samples: {len(all_sat_ids):,}')
    print(f'Total satellites: {len(unique_sats)}')
    
    print()
    print('=== TIMESTAMP RANGE ===')
    all_timestamps = nc.variables['ddm_timestamp_utc'][:]
    print(f'First timestamp: {all_timestamps[0]:.3f} UTC')
    print(f'Last timestamp:  {all_timestamps[-1]:.3f} UTC')
    print(f'Time span: {all_timestamps[-1] - all_timestamps[0]:.3f} seconds')
    
    print()
    print('=== POSITION STATISTICS ===')
    print(f'Position X range: {np.min(pos_x):.1f} to {np.max(pos_x):.1f} m')
    print(f'Position Y range: {np.min(pos_y):.1f} to {np.max(pos_y):.1f} m')
    print(f'Position Z range: {np.min(pos_z):.1f} to {np.max(pos_z):.1f} m')
    
    # Calculate orbital radius
    orbital_radius = np.sqrt(pos_x**2 + pos_y**2 + pos_z**2)
    print(f'Orbital radius range: {np.min(orbital_radius):.1f} to {np.max(orbital_radius):.1f} m')
    print(f'Average orbital radius: {np.mean(orbital_radius):.1f} m')
    
    nc.close()

if __name__ == '__main__':
    display_cygnss_data('cygnss_extracted_data.nc')
