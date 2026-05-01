#!/usr/bin/env python3
"""
Process CYGNSS data to find proximity events and store in a database.

This script:
1. Downloads CYGNSS data for a given date or date range
2. Extracts position data
3. Finds proximity events
4. Appends events to a rolling database NetCDF file

Usage:
    python process_proximity_events.py YYYYMMDD [YYYYMMDD]
    python process_proximity_events.py 20250920
    python process_proximity_events.py 20250920 20250925  # Date range
"""

import os
import sys
import re
import subprocess
import numpy as np
import netCDF4
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta
import argparse
try:
    from dateutil import parser
except ImportError:
    parser = None

# Constants
PROXIMITY_THRESHOLD_KM = 100.0
BASE_URL = "https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2"
SATELLITES = [1, 2, 3, 4, 5, 6, 7, 8]
RAW_DATA_DIR = Path('CYGNSS_data/raw_data')
POS_DATA_DIR = Path('CYGNSS_data/pos_data')
PROXIMITY_DB = Path('CYGNSS_data/proximity_events.nc')

def extract_satellite_id(filename):
    """Extract satellite ID from filename pattern cygXX."""
    match = re.match(r'cyg(\d+)', filename)
    if match:
        return int(match.group(1))
    else:
        raise ValueError(f"Could not extract satellite ID from filename: {filename}")

def parse_timestamp_units(units_str):
    """Parse timestamp units string to extract the reference datetime."""
    match = re.search(r'seconds since (.+)', units_str)
    if match:
        dt_str = match.group(1).strip()
        dt_str = dt_str.replace('Z', '').split('+')[0]
        
        if parser:
            try:
                return parser.parse(dt_str)
            except:
                pass
        
        for fmt in ['%Y-%m-%d %H:%M:%S.%f', '%Y-%m-%d %H:%M:%S', '%Y-%m-%dT%H:%M:%S.%f', '%Y-%m-%dT%H:%M:%S']:
            try:
                dt_str_clean = dt_str.split('.')[0] if '.' in dt_str and '.%f' not in fmt else dt_str
                return datetime.strptime(dt_str_clean, fmt.split('.')[0] if '.%f' not in fmt else fmt)
            except:
                continue
    return None

def download_data_for_date(date_str):
    """Download CYGNSS data for a specific date using the data_import.sh script.
    Returns True if any files were downloaded, False otherwise.
    Continues even if some satellites don't have data (404 errors)."""
    print(f"\n{'='*60}")
    print(f"Downloading data for {date_str}")
    print(f"{'='*60}")
    
    raw_data_dir = RAW_DATA_DIR
    import_script = raw_data_dir / 'data_import.sh'
    
    if not import_script.exists():
        raise FileNotFoundError(f"data_import.sh not found at {import_script}")
    
    # Change to raw_data directory and run the import script
    original_dir = os.getcwd()
    try:
        os.chdir(raw_data_dir)
        result = subprocess.run(
            ['bash', 'data_import.sh', date_str],
            capture_output=False,
            text=True
        )
        # Note: data_import.sh now continues on 404 errors, so non-zero return code
        # might just mean some files weren't available, not a fatal error
    finally:
        os.chdir(original_dir)
    
    # Check if any files were downloaded
    pattern = f"cyg*.ddmi.s{date_str}*.nc"
    downloaded_files = list(raw_data_dir.glob(pattern))
    print(f"Downloaded {len(downloaded_files)} files for {date_str}")
    
    if len(downloaded_files) == 0:
        print(f"Warning: No data files found for {date_str}. This date may not have data available.")
        return False
    
    # Show which satellites were downloaded
    sat_ids = sorted(set([extract_satellite_id(f.name) for f in downloaded_files]))
    print(f"Satellites with data: {sat_ids}")
    
    return True

def extract_data_for_date(date_str):
    """Extract position data using extract_cygnss_data.py script.
    Returns path to extracted data file, or None if extraction failed."""
    print(f"\nExtracting data for {date_str}...")
    
    raw_data_dir = RAW_DATA_DIR
    pos_data_dir = POS_DATA_DIR
    extract_script = pos_data_dir / 'extract_cygnss_data.py'
    
    if not extract_script.exists():
        raise FileNotFoundError(f"extract_cygnss_data.py not found at {extract_script}")
    
    # Use date-specific output filename to avoid conflicts
    output_filename = f'cygnss_extracted_data_{date_str}.nc'
    
    # Run the extraction script
    # The script will place output in pos_data if input is raw_data and output is default name,
    # otherwise in parent of input_dir. We'll check both locations.
    original_dir = os.getcwd()
    try:
        result = subprocess.run(
            ['python3', str(extract_script), str(raw_data_dir), '-o', output_filename],
            capture_output=True,
            text=True
        )
        
        if result.returncode != 0:
            print(f"Error running extraction script:")
            if result.stdout:
                print(result.stdout)
            if result.stderr:
                print(result.stderr)
            return None
        
        # Check possible output locations
        possible_paths = [
            pos_data_dir / output_filename,  # If script logic puts it in pos_data
            raw_data_dir.parent / output_filename,  # If script puts it in parent
        ]
        
        output_path = None
        for path in possible_paths:
            if path.exists():
                output_path = path
                break
        
        if output_path is None:
            print(f"Warning: Extraction script completed but output file not found")
            print(f"  Checked: {possible_paths}")
            if result.stdout:
                print(f"  Script output: {result.stdout}")
            return None
        
        print(f"  Extraction complete: {output_path}")
        return output_path
        
    finally:
        os.chdir(original_dir)

def load_and_process_data(date_str):
    """Load position data from extracted NetCDF file for a specific date.
    Data is already grid-aligned and in kilometers."""
    print(f"\nLoading extracted data for {date_str}...")
    
    # Look for date-specific extracted file (check both possible locations)
    output_filename = f'cygnss_extracted_data_{date_str}.nc'
    possible_paths = [
        POS_DATA_DIR / output_filename,
        RAW_DATA_DIR.parent / output_filename,
    ]
    
    extracted_file = None
    for path in possible_paths:
        if path.exists():
            extracted_file = path
            break
    
    if extracted_file is None:
        print(f"Error: Extracted data file not found for {date_str}")
        print(f"  Checked: {possible_paths}")
        print(f"  Make sure extraction step completed successfully")
        return None
    
    try:
        with netCDF4.Dataset(extracted_file, 'r') as nc:
            timestamps = nc.variables['ddm_timestamp_utc'][:]
            pos_x = nc.variables['sc_pos_x'][:]
            pos_y = nc.variables['sc_pos_y'][:]
            pos_z = nc.variables['sc_pos_z'][:]
            sat_ids = nc.variables['satellite_id'][:]
            
            # Data is already in kilometers and grid-aligned
            df = pd.DataFrame({
                'timestamp': timestamps[:],
                'satellite_id': sat_ids[:],
                'pos_x': pos_x[:],
                'pos_y': pos_y[:],
                'pos_z': pos_z[:]
            })
            
    except Exception as e:
        print(f"  Error reading extracted data file: {e}")
        import traceback
        traceback.print_exc()
        return None
    
    # Check for invalid values
    # NaN values are expected and intentional - they indicate missing data when no timestamp
    # exists within tolerance of a grid point. These will be skipped in distance calculations.
    # Infinite values are errors and should not exist.
    initial_count = len(df)
    nan_count = df[['pos_x', 'pos_y', 'pos_z']].isna().any(axis=1).sum()
    inf_count = (np.isinf(df['pos_x']) | np.isinf(df['pos_y']) | np.isinf(df['pos_z'])).sum()
    
    if nan_count > 0:
        print(f"  Note: Found {nan_count} rows with NaN position data (expected - indicates missing data at grid points)")
        print(f"  These entries will be skipped in proximity calculations to maintain consistent indexing")
    
    if inf_count > 0:
        print(f"  WARNING: Found {inf_count} rows with infinite values - this should not happen!")
        # Filter out infinite values (these are errors)
        df = df[~np.isinf(df['pos_x']) & ~np.isinf(df['pos_y']) & ~np.isinf(df['pos_z'])]
        print(f"  Filtered out {initial_count - len(df)} rows with infinite position data")
    
    # Keep NaN entries for consistent indexing - they will be skipped in distance calculations
    print(f"  Loaded {len(df)} data points ({initial_count - nan_count} with valid positions, {nan_count} with NaN)")
    return df

def find_proximity_events(df, threshold_km=100.0):
    """Find proximity events and group consecutive events."""
    print(f"\nFinding proximity events within {threshold_km} km...")
    
    all_instances = []
    unique_times = sorted(df['timestamp'].unique())
    
    print(f"  Checking {len(unique_times)} unique timestamps...")
    
    for timestamp in unique_times:
        time_data = df[df['timestamp'] == timestamp]
        
        if len(time_data) < 2:
            continue
        
        positions = time_data[['pos_x', 'pos_y', 'pos_z']].values
        sat_ids = time_data['satellite_id'].values
        
        # Check each pair only once
        for i in range(len(positions)):
            for j in range(i + 1, len(positions)):
                pos_i = positions[i]
                pos_j = positions[j]
                
                # Skip if either position has NaN (indicates missing data for that satellite at this timestamp)
                # NaN entries are inserted when no data exists within tolerance of the grid point
                if np.any(np.isnan(pos_i)) or np.any(np.isnan(pos_j)):
                    # Silently skip - this is expected when satellites don't have data at certain timestamps
                    continue
                if np.any(np.isinf(pos_i)) or np.any(np.isinf(pos_j)):
                    print(f"    WARNING: Found infinite value in position data at timestamp {timestamp}")
                    continue
                
                # Calculate distance (positions are already in km, so distance will be in km)
                diff = pos_i - pos_j
                dist_squared = np.sum(diff**2)
                
                # Handle potential floating point errors (shouldn't happen, but be safe)
                if np.isnan(dist_squared):
                    print(f"    WARNING: NaN in distance calculation at timestamp {timestamp}")
                    continue
                if dist_squared < 0:
                    # Floating point error - set to zero
                    print(f"    WARNING: Negative distance squared ({dist_squared}) at timestamp {timestamp}, setting to 0")
                    dist_squared = 0.0
                
                distance_km = np.sqrt(dist_squared)
                
                if distance_km <= threshold_km:
                    sat_1 = int(sat_ids[i])
                    sat_2 = int(sat_ids[j])
                    if sat_1 > sat_2:
                        sat_1, sat_2 = sat_2, sat_1
                    
                    all_instances.append({
                        'timestamp': timestamp,
                        'satellite_1': sat_1,
                        'satellite_2': sat_2,
                        'distance_km': distance_km
                    })
    
    if len(all_instances) == 0:
        print("  No proximity events found")
        return pd.DataFrame(columns=['satellite_1', 'satellite_2', 'start_time', 'end_time', 'min_distance_km'])
    
    # Group consecutive events
    instances_df = pd.DataFrame(all_instances)
    instances_df = instances_df.sort_values(['satellite_1', 'satellite_2', 'timestamp']).reset_index(drop=True)
    
    grouped_events = []
    
    for (sat_1, sat_2), group in instances_df.groupby(['satellite_1', 'satellite_2']):
        group = group.sort_values('timestamp').reset_index(drop=True)
        
        start_time = group.iloc[0]['timestamp']
        end_time = group.iloc[0]['timestamp']
        min_distance = group.iloc[0]['distance_km']
        
        for idx in range(1, len(group)):
            current_time = group.iloc[idx]['timestamp']
            prev_time = group.iloc[idx - 1]['timestamp']
            
            time_gap = current_time - prev_time
            # Data is sampled every 5 minutes (300 seconds), so consecutive timestamps
            # should be ~300 seconds apart. Allow some tolerance for slight variations.
            max_gap = 350.0  # 350 seconds (5 min 50 sec) to account for 5-minute intervals
            
            if time_gap <= max_gap:
                end_time = current_time
                min_distance = min(min_distance, group.iloc[idx]['distance_km'])
            else:
                # Gap is too large - start a new event
                grouped_events.append({
                    'satellite_1': sat_1,
                    'satellite_2': sat_2,
                    'start_time': start_time,
                    'end_time': end_time,
                    'min_distance_km': min_distance
                })
                start_time = current_time
                end_time = current_time
                min_distance = group.iloc[idx]['distance_km']
        
        grouped_events.append({
            'satellite_1': sat_1,
            'satellite_2': sat_2,
            'start_time': start_time,
            'end_time': end_time,
            'min_distance_km': min_distance
        })
    
    proximity_df = pd.DataFrame(grouped_events)
    proximity_df = proximity_df.sort_values('start_time').reset_index(drop=True)
    
    print(f"  Found {len(instances_df)} individual instances")
    print(f"  Grouped into {len(proximity_df)} proximity events")
    return proximity_df

def create_proximity_database(output_path, threshold_km=100.0):
    """Create a new proximity events database NetCDF file."""
    print(f"\nCreating new proximity database: {output_path}")
    
    with netCDF4.Dataset(output_path, 'w', format='NETCDF4') as nc:
        nc.createDimension('event', None)  # Unlimited dimension
        
        # Variables
        sat1_var = nc.createVariable('satellite_1', 'i4', ('event',))
        sat1_var.long_name = 'First satellite ID'
        
        sat2_var = nc.createVariable('satellite_2', 'i4', ('event',))
        sat2_var.long_name = 'Second satellite ID'
        
        start_time_var = nc.createVariable('start_time', 'f8', ('event',))
        start_time_var.units = 'seconds since 2000-01-01 00:00:00 UTC'
        start_time_var.long_name = 'Event start time'
        
        end_time_var = nc.createVariable('end_time', 'f8', ('event',))
        end_time_var.units = 'seconds since 2000-01-01 00:00:00 UTC'
        end_time_var.long_name = 'Event end time'
        
        min_dist_var = nc.createVariable('min_distance_km', 'f8', ('event',))
        min_dist_var.units = 'km'
        min_dist_var.long_name = 'Minimum distance during event'
        
        # Global attributes
        nc.title = 'CYGNSS Proximity Events Database'
        nc.description = 'Rolling database of proximity events between CYGNSS satellites'
        nc.history = f'Created on {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}'
        nc.proximity_threshold_km = threshold_km
        nc.total_events = 0

def append_to_database(proximity_df, db_path, threshold_km=100.0):
    """Append proximity events to the database file."""
    if len(proximity_df) == 0:
        print("  No events to append")
        return
    
    # Create database if it doesn't exist
    if not db_path.exists():
        create_proximity_database(db_path, threshold_km)
    
    print(f"\nAppending {len(proximity_df)} events to database...")
    
    with netCDF4.Dataset(db_path, 'a', format='NETCDF4') as nc:
        current_size = nc.dimensions['event'].size
        
        # Append data
        new_size = current_size + len(proximity_df)
        
        nc.variables['satellite_1'][current_size:new_size] = proximity_df['satellite_1'].values
        nc.variables['satellite_2'][current_size:new_size] = proximity_df['satellite_2'].values
        nc.variables['start_time'][current_size:new_size] = proximity_df['start_time'].values
        nc.variables['end_time'][current_size:new_size] = proximity_df['end_time'].values
        nc.variables['min_distance_km'][current_size:new_size] = proximity_df['min_distance_km'].values
        
        # Update total events
        nc.total_events = new_size
        nc.history = f"{nc.history}; Updated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
    
    print(f"  Database now contains {new_size} total events")

def cleanup_raw_data(date_str):
    """Delete raw data files for a specific date."""
    print(f"\nCleaning up raw data files for {date_str}...")
    
    raw_data_dir = RAW_DATA_DIR
    pattern = f"cyg*.ddmi.s{date_str}*.nc"
    files_to_delete = list(raw_data_dir.glob(pattern))
    
    if files_to_delete:
        for file in files_to_delete:
            file.unlink()
        print(f"  Deleted {len(files_to_delete)} raw data files")
    else:
        print("  No files to delete")

def process_date(date_str, threshold_km=100.0, download=True, cleanup=True):
    """Process a single date: download, extract, find events, append to database.
    Returns True if processing succeeded (even if no proximity events found),
    False if no data was available."""
    print(f"\n{'='*60}")
    print(f"Processing date: {date_str}")
    print(f"{'='*60}")
    
    # Step 1: Download
    if download:
        if not download_data_for_date(date_str):
            print(f"No data available for {date_str} - skipping this date")
            return False
    else:
        print(f"Skipping download (using existing data)")
    
    # Step 2: Extract data using extract_cygnss_data.py
    extracted_file = extract_data_for_date(date_str)
    if extracted_file is None:
        print(f"Extraction failed for {date_str} - skipping this date")
        return False
    
    # Step 3: Load extracted data
    df = load_and_process_data(date_str)
    if df is None or len(df) == 0:
        print(f"No data loaded for {date_str} - skipping this date")
        return False
    
    # Step 4: Find proximity events
    proximity_df = find_proximity_events(df, threshold_km)
    
    # Step 5: Append to database
    if len(proximity_df) > 0:
        append_to_database(proximity_df, PROXIMITY_DB, threshold_km)
    else:
        print("No proximity events to add to database")
    
    # Step 6: Cleanup raw data files
    if cleanup:
        cleanup_raw_data(date_str)
    
    # Also cleanup the extracted file for this date (to save space)
    try:
        extracted_file.unlink()
        print(f"  Cleaned up extracted data file: {extracted_file.name}")
    except Exception as e:
        print(f"  Warning: Could not delete extracted file: {e}")
    
    return True

def date_range(start_date_str, end_date_str):
    """Generate list of dates between start and end (inclusive)."""
    start = datetime.strptime(start_date_str, '%Y%m%d')
    end = datetime.strptime(end_date_str, '%Y%m%d')
    
    dates = []
    current = start
    while current <= end:
        dates.append(current.strftime('%Y%m%d'))
        current += timedelta(days=1)
    
    return dates

def main():
    parser = argparse.ArgumentParser(
        description='Process CYGNSS data to find proximity events',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python process_proximity_events.py 20250920
  python process_proximity_events.py 20250920 20250925
  python process_proximity_events.py 20250920 --no-download --no-cleanup
        """
    )
    parser.add_argument('start_date', help='Start date in YYYYMMDD format')
    parser.add_argument('end_date', nargs='?', help='End date in YYYYMMDD format (optional, defaults to start_date)')
    parser.add_argument('--no-download', action='store_true', help='Skip download step (use existing raw data)')
    parser.add_argument('--no-cleanup', action='store_true', help='Skip cleanup step (keep raw data files)')
    parser.add_argument('--threshold', type=float, default=PROXIMITY_THRESHOLD_KM,
                       help=f'Proximity threshold in km (default: {PROXIMITY_THRESHOLD_KM})')
    
    args = parser.parse_args()
    
    # Validate date format
    if not re.match(r'^\d{8}$', args.start_date):
        print("Error: Start date must be in YYYYMMDD format")
        return 1
    
    if args.end_date and not re.match(r'^\d{8}$', args.end_date):
        print("Error: End date must be in YYYYMMDD format")
        return 1
    
    # Determine dates to process
    if args.end_date:
        dates = date_range(args.start_date, args.end_date)
        print(f"Processing date range: {args.start_date} to {args.end_date} ({len(dates)} days)")
    else:
        dates = [args.start_date]
        print(f"Processing single date: {args.start_date}")
    
    # Get threshold
    threshold_km = args.threshold
    
    # Ensure directories exist
    RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)
    POS_DATA_DIR.mkdir(parents=True, exist_ok=True)
    PROXIMITY_DB.parent.mkdir(parents=True, exist_ok=True)
    
    # Process each date
    success_count = 0
    missing_dates = []
    error_dates = []
    
    for date_str in dates:
        try:
            if process_date(date_str, threshold_km=threshold_km, download=not args.no_download, cleanup=not args.no_cleanup):
                success_count += 1
            else:
                # Date was processed but had no data available
                missing_dates.append(date_str)
        except Exception as e:
            print(f"\nError processing {date_str}: {e}")
            import traceback
            traceback.print_exc()
            error_dates.append(date_str)
            continue
    
    print(f"\n{'='*60}")
    print(f"Processing Summary")
    print(f"{'='*60}")
    print(f"Total dates requested: {len(dates)}")
    print(f"Successfully processed: {success_count}")
    if missing_dates:
        print(f"Dates with no data available: {len(missing_dates)}")
        print(f"  {', '.join(missing_dates)}")
    if error_dates:
        print(f"Dates with errors: {len(error_dates)}")
        print(f"  {', '.join(error_dates)}")
    print(f"\nDatabase location: {PROXIMITY_DB}")
    print(f"{'='*60}\n")
    
    return 0 if success_count == len(dates) else 1

if __name__ == '__main__':
    sys.exit(main())
