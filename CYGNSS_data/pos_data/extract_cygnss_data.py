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
from datetime import datetime, timedelta
try:
    from dateutil import parser
except ImportError:
    parser = None, timedelta
try:
    from dateutil import parser
except ImportError:
    print("Warning: python-dateutil not installed. Install with: pip install python-dateutil")
    parser = None

def extract_satellite_id(filename):
    """Extract satellite ID from filename pattern cygXX."""
    match = re.match(r'cyg(\d+)', filename)
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"Could not extract satellite ID from filename: {filename}")

def parse_timestamp_units(units_str):
    """Parse timestamp units string to extract the reference datetime.
    Format: 'seconds since YYYY-MM-DD HH:MM:SS...'"""
    # Extract the datetime part after "seconds since "
    match = re.search(r'seconds since (.+)', units_str)
    if match:
        dt_str = match.group(1).strip()
        # Remove timezone info if present
        dt_str = dt_str.replace('Z', '').split('+')[0]
        if dt_str.count('-') > 2:  # Has timezone offset
            dt_str = '-'.join(dt_str.split('-')[:3]) + ' ' + dt_str.split()[-1] if ' ' in dt_str else dt_str
        
        # Try parsing with dateutil if available
        if parser:
            try:
                return parser.parse(dt_str)
            except:
                pass
        
        # Fallback: try to parse common formats
        for fmt in ['%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S', '%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S']:
            try:
                # Remove microseconds if present and format doesn't support it
                if '.%f' not in fmt and '.' in dt_str:
                    dt_str_clean = dt_str.split('.')[0]
                else:
                    dt_str_clean = dt_str
                return datetime.strptime(dt_str_clean, fmt.split('.')[0] if '.%f' not in fmt else fmt)
            except:
                continue
    return None

def parse_timestamp_units(units_str):
    """Parse timestamp units string to extract the reference datetime.
    Format: 'seconds since YYYY-MM-DD HH:MM:SS...'"""
    # Extract the datetime part after "seconds since "
    match = re.search(r'seconds since (.+)', units_str)
    if match:
        dt_str = match.group(1).strip()
        # Remove timezone info if present
        dt_str = dt_str.replace('Z', '').split('+')[0].split('-')[0] if '-' in dt_str.split()[0] else dt_str
        
        # Try parsing with dateutil if available
        if parser:
            try:
                return parser.parse(dt_str)
            except:
                pass
        
        # Fallback: try to parse common formats
        for fmt in ['%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S', '%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S']:
            try:
                return datetime.strptime(dt_str.split('.')[0] if '.' in dt_str else dt_str, fmt.split('.')[0])
            except:
                continue
    return None

# ============================================================================
# TIMESTAMP TOLERANCE CONFIGURATION
# ============================================================================
# Tolerance for matching timestamps to grid points (in seconds).
# If no timestamp is within this tolerance of a grid point, NaN entries
# will be inserted for that satellite at that index.
# Change this value to adjust the tolerance threshold.
TIMESTAMP_TOLERANCE_SECONDS = 10.0
# ============================================================================

def find_closest_timestamps(timestamps, target_times, tolerance=TIMESTAMP_TOLERANCE_SECONDS):
    """
    Find the timestamp closest to each grid point's target time.
    
    For each grid point in the fixed time grid (0:00:00, 0:05:00, ..., 23:55:00),
    searches for the closest timestamp in the raw data. Always processes every grid point
    in order and always advances to the next grid point, preventing phase drift.
    
    Args:
        timestamps: Sorted array of timestamps (in chronological order)
        target_times: Array of target times (grid points, 5 minutes apart, same for all satellites)
        tolerance: Maximum allowed difference in seconds (default: TIMESTAMP_TOLERANCE_SECONDS)
    
    Returns:
        Tuple of (indices, valid_mask) where:
        - indices: Array of indices into timestamps array (or -1 if no match within tolerance)
        - valid_mask: Boolean array indicating which indices are within tolerance
    """
    indices = []
    valid_mask = []
    
    # Expected sampling rate: ~0.5 seconds, so ~600 indices per 5-minute grid point
    expected_step = 600
    
    for grid_idx, target in enumerate(target_times):
        # Estimate where this grid point's timestamp should be (for efficient search)
        estimated_idx = grid_idx * expected_step
        
        # Search window: check around the estimated position
        # Use a generous window (±10 minutes = ±1200 indices) to handle variations
        search_window = 1200
        start_idx = max(0, estimated_idx - search_window)
        end_idx = min(len(timestamps), estimated_idx + search_window)
        
        if start_idx >= len(timestamps):
            # Past the end of data
            indices.append(-1)
            valid_mask.append(False)
            continue
        
        if end_idx <= start_idx:
            # No data in search window
            indices.append(-1)
            valid_mask.append(False)
            continue
        
        # Binary search to find insertion point for target in the search window
        window_timestamps = timestamps[start_idx:end_idx]
        local_insert_idx = np.searchsorted(window_timestamps, target, side='left')
        global_insert_idx = start_idx + local_insert_idx
        
        # Check the timestamp at insertion point and the one before it (if any)
        candidates = []
        
        # Check timestamp at or after target
        if global_insert_idx < len(timestamps):
            dist = abs(timestamps[global_insert_idx] - target)
            candidates.append((global_insert_idx, dist))
        
        # Check timestamp before target (if within search window)
        if global_insert_idx > start_idx and global_insert_idx > 0:
            dist = abs(timestamps[global_insert_idx - 1] - target)
            candidates.append((global_insert_idx - 1, dist))
        
        if not candidates:
            # No candidates found (shouldn't happen, but handle gracefully)
            indices.append(-1)
            valid_mask.append(False)
            continue
        
        # Find the closest candidate
        closest_idx, min_dist = min(candidates, key=lambda x: x[1])
        
        # Check if within tolerance
        if min_dist <= tolerance:
            indices.append(closest_idx)
            valid_mask.append(True)
        else:
            # No match within tolerance - will use NaN for positions, exact grid time for timestamp
            indices.append(-1)
            valid_mask.append(False)
    
    return np.array(indices), np.array(valid_mask)

def get_timestamps_from_file(filepath):
    """Extract and convert timestamps from a NetCDF file (first pass to determine grid)."""
    try:
        with netCDF4.Dataset(filepath, 'r') as nc:
            timestamp_var = nc.variables['ddm_timestamp_utc']
            timestamp = timestamp_var[:]
            
            # Get the reference datetime
            reference_datetime = None
            
            # Method 1: Use time_coverage_start (most reliable)
            if 'time_coverage_start' in nc.ncattrs():
                try:
                    tc_start = nc.time_coverage_start
                    if 'T' in tc_start:
                        dt_str = tc_start.replace('Z', '').split('.')[0]
                        reference_datetime = datetime.strptime(dt_str, '%Y-%m-%dT%H:%M:%S')
                    else:
                        reference_datetime = datetime.fromisoformat(tc_start.split('.')[0])
                except Exception as e:
                    pass
            
            # Method 2: Parse from units string
            if reference_datetime is None:
                units_str = getattr(timestamp_var, 'units', '')
                if units_str:
                    reference_datetime = parse_timestamp_units(units_str)
            
            # Method 3: Extract date from filename (fallback)
            if reference_datetime is None:
                match = re.search(r's(\d{8})', os.path.basename(filepath))
                if match:
                    date_str = match.group(1)
                    reference_datetime = datetime.strptime(date_str, '%Y%m%d')
                else:
                    return None, None
            
            # Convert to absolute timestamps
            epoch_2000 = datetime(2000, 1, 1, 0, 0, 0)
            absolute_timestamps = np.array([
                (reference_datetime + timedelta(seconds=float(ts)) - epoch_2000).total_seconds()
                for ts in timestamp
            ])
            
            return absolute_timestamps, reference_datetime
            
    except Exception as e:
        print(f"Error reading timestamps from {filepath}: {e}")
        return None, None

def process_netcdf_file(filepath, time_grid=None):
    """
    Extract required fields from a single NetCDF file.
    If time_grid is provided, extracts only data at timestamps closest to grid points.
    """
    print(f"Processing: {filepath}")
    
    try:
        with netCDF4.Dataset(filepath, 'r') as nc:
            # Extract required fields
            timestamp_var = nc.variables['ddm_timestamp_utc']
            timestamp = timestamp_var[:]
            
            # Get the reference datetime - try multiple methods
            reference_datetime = None
            
            # Method 1: Use time_coverage_start (most reliable)
            if 'time_coverage_start' in nc.ncattrs():
                try:
                    tc_start = nc.time_coverage_start
                    # Handle ISO format with or without timezone
                    if 'T' in tc_start:
                        # ISO format: 2025-09-20T00:00:00.499261738Z
                        dt_str = tc_start.replace('Z', '').split('.')[0]  # Remove microseconds and Z
                        reference_datetime = datetime.strptime(dt_str, '%Y-%m-%dT%H:%M:%S')
                    else:
                        reference_datetime = datetime.fromisoformat(tc_start.split('.')[0])
                except Exception as e:
                    print(f"  Warning: Could not parse time_coverage_start: {e}")
            
            # Method 2: Parse from units string
            if reference_datetime is None:
                units_str = getattr(timestamp_var, 'units', '')
                if units_str:
                    reference_datetime = parse_timestamp_units(units_str)
            
            # Method 3: Extract date from filename (fallback)
            if reference_datetime is None:
                # Format: cygXX.ddmi.sYYYYMMDD-...
                match = re.search(r's(\d{8})', os.path.basename(filepath))
                if match:
                    date_str = match.group(1)
                    reference_datetime = datetime.strptime(date_str, '%Y%m%d')
                else:
                    raise ValueError(f"Could not determine reference datetime for {filepath}")
            
            # Convert relative timestamps to absolute UTC timestamps
            epoch_2000 = datetime(2000, 1, 1, 0, 0, 0)
            absolute_timestamps = np.array([
                (reference_datetime + timedelta(seconds=float(ts)) - epoch_2000).total_seconds()
                for ts in timestamp
            ])
            
            pos_x = nc.variables['sc_pos_x'][:]
            pos_y = nc.variables['sc_pos_y'][:]
            pos_z = nc.variables['sc_pos_z'][:]
            
            # Convert positions from meters to kilometers
            pos_x = np.array(pos_x, dtype=np.float64) / 1000.0
            pos_y = np.array(pos_y, dtype=np.float64) / 1000.0
            pos_z = np.array(pos_z, dtype=np.float64) / 1000.0
            
            # Get satellite ID from filename
            satellite_id = extract_satellite_id(os.path.basename(filepath))
            
            # If time_grid is provided, find closest timestamps and sample
            if time_grid is not None:
                # Find indices of timestamps closest to each grid point (with tolerance)
                sample_indices, valid_mask = find_closest_timestamps(absolute_timestamps, time_grid, tolerance=TIMESTAMP_TOLERANCE_SECONDS)
                
                # Create arrays for all grid points
                grid_timestamps = np.full(len(time_grid), np.nan)
                grid_pos_x = np.full(len(time_grid), np.nan)
                grid_pos_y = np.full(len(time_grid), np.nan)
                grid_pos_z = np.full(len(time_grid), np.nan)
                
                # Fill in valid entries with actual matched data (use real timestamps)
                valid_indices = sample_indices[valid_mask]
                grid_timestamps[valid_mask] = absolute_timestamps[valid_indices]  # Use actual matched timestamp
                grid_pos_x[valid_mask] = pos_x[valid_indices]
                grid_pos_y[valid_mask] = pos_y[valid_indices]
                grid_pos_z[valid_mask] = pos_z[valid_indices]
                
                # For invalid matches (NaN entries), use exact grid point time for alignment
                grid_timestamps[~valid_mask] = time_grid[~valid_mask]
                # Positions remain NaN for invalid entries
                
                # Use grid timestamps (which may include NaN for missing points)
                absolute_timestamps = grid_timestamps
                pos_x = grid_pos_x
                pos_y = grid_pos_y
                pos_z = grid_pos_z
                sat_id_array = np.full(len(time_grid), satellite_id)
                
                original_count = len(timestamp)
                valid_count = np.sum(valid_mask)
                nan_count = len(time_grid) - valid_count
                reduction_factor = original_count / valid_count if valid_count > 0 else 1.0
                
                print(f"  Reference datetime: {reference_datetime}")
                print(f"  Original samples: {original_count}, Grid-aligned samples: {valid_count} valid, {nan_count} NaN (reduction: {reduction_factor:.1f}x)")
                if valid_count > 0:
                    valid_times = grid_timestamps[valid_mask]
                    print(f"  Timestamp range: {valid_times[0]:.3f} to {valid_times[-1]:.3f} seconds since 2000-01-01")
            else:
                # No grid provided - return all data (for first pass)
                sat_id_array = np.full(len(timestamp), satellite_id)
                print(f"  Reference datetime: {reference_datetime}")
                print(f"  Total samples: {len(timestamp)}")
                print(f"  Timestamp range: {absolute_timestamps[0]:.3f} to {absolute_timestamps[-1]:.3f} seconds since 2000-01-01")
            
            return {
                'timestamp': absolute_timestamps,
                'pos_x': pos_x,
                'pos_y': pos_y,
                'pos_z': pos_z,
                'satellite_id': sat_id_array,
                'filename': os.path.basename(filepath),
                'reference_datetime': reference_datetime
            }
            
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        import traceback
        traceback.print_exc()
        return None

def create_output_netcdf(data_list, output_path):
    """Create output NetCDF file with all extracted data, organized by timestamp for consistent indexing."""
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
    
    # Sort by timestamp, then by satellite ID to ensure consistent indexing
    # This groups data from the same time together, making it easy to compare across satellites
    sort_indices = np.lexsort((sat_ids, timestamps))
    timestamps = timestamps[sort_indices]
    pos_x = pos_x[sort_indices]
    pos_y = pos_y[sort_indices]
    pos_z = pos_z[sort_indices]
    sat_ids = sat_ids[sort_indices]
    filenames = filenames[sort_indices]
    
    # Create output NetCDF file
    with netCDF4.Dataset(output_path, 'w', format='NETCDF4') as nc:
        # Create dimensions
        nc.createDimension('sample', len(timestamps))
        
        # Create variables
        timestamp_var = nc.createVariable('ddm_timestamp_utc', 'f8', ('sample',))
        timestamp_var.units = 'seconds since 2000-01-01 00:00:00'
        timestamp_var.long_name = 'DDM timestamp in UTC'
        
        pos_x_var = nc.createVariable('sc_pos_x', 'f8', ('sample',))
        pos_x_var.units = 'km'
        pos_x_var.long_name = 'Spacecraft position X component'
        
        pos_y_var = nc.createVariable('sc_pos_y', 'f8', ('sample',))
        pos_y_var.units = 'km'
        pos_y_var.long_name = 'Spacecraft position Y component'
        
        pos_z_var = nc.createVariable('sc_pos_z', 'f8', ('sample',))
        pos_z_var.units = 'km'
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
    
    # ===== FIRST PASS: Determine common time grid =====
    print("\n" + "="*60)
    print("First pass: Determining common time grid")
    print("="*60)
    
    SAMPLE_INTERVAL_SECONDS = 300.0  # 5 minutes
    
    earliest_start = None
    latest_end = None
    
    for nc_file in sorted(nc_files):
        timestamps, ref_dt = get_timestamps_from_file(nc_file)
        if timestamps is not None and len(timestamps) > 0:
            file_start = timestamps[0]
            file_end = timestamps[-1]
            
            if earliest_start is None or file_start < earliest_start:
                earliest_start = file_start
            if latest_end is None or file_end > latest_end:
                latest_end = file_end
    
    if earliest_start is None:
        print("Error: Could not determine time range from any files")
        return 1
    
    # Create time grid: every 5 minutes starting from the earliest start time
    # Round earliest_start down to nearest 5-minute mark
    grid_start = np.floor(earliest_start / SAMPLE_INTERVAL_SECONDS) * SAMPLE_INTERVAL_SECONDS
    grid_end = latest_end
    
    # Generate grid points
    time_grid = np.arange(grid_start, grid_end + SAMPLE_INTERVAL_SECONDS, SAMPLE_INTERVAL_SECONDS)
    
    print(f"Earliest timestamp: {earliest_start:.3f} seconds since 2000-01-01")
    print(f"Latest timestamp: {latest_end:.3f} seconds since 2000-01-01")
    print(f"Grid start: {grid_start:.3f} seconds since 2000-01-01")
    print(f"Grid points: {len(time_grid)} (every {SAMPLE_INTERVAL_SECONDS} seconds)")
    
    # ===== SECOND PASS: Extract data at grid-aligned timestamps =====
    print("\n" + "="*60)
    print("Second pass: Extracting data at grid-aligned timestamps")
    print("="*60)
    
    data_list = []
    for nc_file in sorted(nc_files):
        data = process_netcdf_file(nc_file, time_grid=time_grid)
        data_list.append(data)
    
    # Filter out None results
    valid_data = [d for d in data_list if d is not None]
    
    if not valid_data:
        print("No valid data found")
        return 1
    
    print(f"\nSuccessfully processed {len(valid_data)} files")
    
    # Create output file (default to pos_data directory if output is default name)
    if args.output == 'cygnss_extracted_data.nc' and input_dir.name == 'raw_data':
        output_path = input_dir.parent / 'pos_data' / args.output
    else:
        output_path = input_dir.parent / args.output
    create_output_netcdf(valid_data, output_path)
    
    print(f"\nExtraction complete! Output saved to: {output_path}")
    return 0

if __name__ == '__main__':
    
    exit(main())
