#!/usr/bin/env python3
"""
Extract specific fields from CYGNSS NetCDF files and organize by satellite ID.

This script extracts the following fields from each .nc file:
- ddm_timestamp_utc
- sc_pos_x
- sc_pos_y  
- sc_pos_z

The satellite ID is extracted from the filename pattern cygXX and added as a field.
All data is combined into a single output NetCDF file organized by satellite.
"""

import os
import re
import numpy as np
import netCDF4
from pathlib import Path
import argparse
from datetime import datetime

def extract_satellite_id(filename):
    """Extract satellite ID from filename pattern cygXX."""
    match = re.match(r'cyg(\d+)', filename)
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"Could not extract satellite ID from filename: {filename}")

def process_netcdf_file(filepath):
    """Extract required fields from a single NetCDF file."""
    print(f"Processing: {filepath}")
    
    try:
        with netCDF4.Dataset(filepath, 'r') as nc:
            # Extract required fields
            timestamp = nc.variables['ddm_timestamp_utc'][:]
            pos_x = nc.variables['sc_pos_x'][:]
            pos_y = nc.variables['sc_pos_y'][:]
            pos_z = nc.variables['sc_pos_z'][:]
            
            # Get satellite ID from filename
            satellite_id = extract_satellite_id(os.path.basename(filepath))
            
            # Create satellite ID array
            sat_id_array = np.full(len(timestamp), satellite_id)
            
            return {
                'timestamp': timestamp,
                'pos_x': pos_x,
                'pos_y': pos_y,
                'pos_z': pos_z,
                'satellite_id': sat_id_array,
                'filename': os.path.basename(filepath)
            }
            
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return None

def create_output_netcdf(data_list, output_path):
    """Create output NetCDF file with all extracted data."""
    print(f"Creating output file: {output_path}")
    
    # Combine all data
    all_timestamps = []
    all_pos_x = []
    all_pos_y = []
    all_pos_z = []
    all_sat_ids = []
    all_filenames = []
    
    for data in data_list:
        if data is not None:
            all_timestamps.extend(data['timestamp'])
            all_pos_x.extend(data['pos_x'])
            all_pos_y.extend(data['pos_y'])
            all_pos_z.extend(data['pos_z'])
            all_sat_ids.extend(data['satellite_id'])
            all_filenames.extend([data['filename']] * len(data['timestamp']))
    
    # Convert to numpy arrays
    timestamps = np.array(all_timestamps)
    pos_x = np.array(all_pos_x)
    pos_y = np.array(all_pos_y)
    pos_z = np.array(all_pos_z)
    sat_ids = np.array(all_sat_ids)
    filenames = np.array(all_filenames, dtype='S100')  # String array
    
    # Create output NetCDF file
    with netCDF4.Dataset(output_path, 'w', format='NETCDF4') as nc:
        # Create dimensions
        nc.createDimension('sample', len(timestamps))
        
        # Create variables
        timestamp_var = nc.createVariable('ddm_timestamp_utc', 'f8', ('sample',))
        timestamp_var.units = 'seconds since 2000-01-01 00:00:00'
        timestamp_var.long_name = 'DDM timestamp in UTC'
        
        pos_x_var = nc.createVariable('sc_pos_x', 'f8', ('sample',))
        pos_x_var.units = 'm'
        pos_x_var.long_name = 'Spacecraft position X component'
        
        pos_y_var = nc.createVariable('sc_pos_y', 'f8', ('sample',))
        pos_y_var.units = 'm'
        pos_y_var.long_name = 'Spacecraft position Y component'
        
        pos_z_var = nc.createVariable('sc_pos_z', 'f8', ('sample',))
        pos_z_var.units = 'm'
        pos_z_var.long_name = 'Spacecraft position Z component'
        
        sat_id_var = nc.createVariable('satellite_id', 'i4', ('sample',))
        sat_id_var.long_name = 'CYGNSS satellite ID (cygXX)'
        sat_id_var.units = 'dimensionless'
        
        # Assign data
        timestamp_var[:] = timestamps
        pos_x_var[:] = pos_x
        pos_y_var[:] = pos_y
        pos_z_var[:] = pos_z
        sat_id_var[:] = sat_ids
        
        # Handle filename strings
        max_filename_len = max(len(f) for f in filenames)
        nc.createDimension('filename_len', max_filename_len)
        filename_var = nc.createVariable('source_filename', 'S1', ('sample', 'filename_len'))
        filename_var.long_name = 'Source NetCDF filename'
        
        # Fill filename array
        filename_array = np.zeros((len(filenames), max_filename_len), dtype='S1')
        for i, filename in enumerate(filenames):
            if isinstance(filename, bytes):
                filename_bytes = filename
            else:
                filename_bytes = filename.encode('utf-8')
            filename_array[i, :len(filename_bytes)] = [c for c in filename_bytes]
        filename_var[:] = filename_array
        
        # Add global attributes
        nc.title = 'CYGNSS Extracted Data - Timestamps and Spacecraft Positions'
        nc.description = 'Extracted fields from CYGNSS Level 1 NetCDF files'
        nc.history = f'Created on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
        nc.source = 'CYGNSS Level 1 Science Data Record Version 3.2'
        nc.institution = 'University of Michigan Space Physics Research Lab (SPRL)'
        nc.total_samples = len(timestamps)
        nc.satellites_included = sorted(list(set(sat_ids)))
        
        print(f"Output file created with {len(timestamps)} total samples")
        print(f"Satellites included: {sorted(list(set(sat_ids)))}")

def main():
    parser = argparse.ArgumentParser(description='Extract CYGNSS data fields and organize by satellite')
    parser.add_argument('input_dir', help='Directory containing .nc files')
    parser.add_argument('-o', '--output', default='cygnss_extracted_data.nc', 
                       help='Output NetCDF file name')
    parser.add_argument('--pattern', default='*.nc', 
                       help='File pattern to match (default: *.nc)')
    
    args = parser.parse_args()
    
    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        print(f"Error: Input directory {input_dir} does not exist")
        return 1
    
    # Find all .nc files
    nc_files = list(input_dir.glob(args.pattern))
    if not nc_files:
        print(f"No .nc files found in {input_dir}")
        return 1
    
    print(f"Found {len(nc_files)} NetCDF files to process")
    
    # Process all files
    data_list = []
    for nc_file in sorted(nc_files):
        data = process_netcdf_file(nc_file)
        data_list.append(data)
    
    # Filter out None results
    valid_data = [d for d in data_list if d is not None]
    
    if not valid_data:
        print("No valid data found")
        return 1
    
    print(f"Successfully processed {len(valid_data)} files")
    
    # Create output file
    output_path = input_dir.parent / args.output
    create_output_netcdf(valid_data, output_path)
    
    print(f"Extraction complete! Output saved to: {output_path}")
    return 0

if __name__ == '__main__':
    
    exit(main())
