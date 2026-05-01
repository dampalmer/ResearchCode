#!/usr/bin/env python3
"""
Test script to check if timestamps align across satellites when using every 600th point.
"""

import numpy as np
import netCDF4
from pathlib import Path
from datetime import datetime, timedelta

def test_timestamp_alignment(extracted_file):
    """Check if timestamps align across satellites in the extracted data."""
    
    print(f"Testing timestamp alignment in: {extracted_file}")
    print("="*60)
    
    with netCDF4.Dataset(extracted_file, 'r') as nc:
        timestamps = nc.variables['ddm_timestamp_utc'][:]
        sat_ids = nc.variables['satellite_id'][:]
        
        # Group timestamps by satellite
        satellites = sorted(set(sat_ids))
        print(f"\nFound {len(satellites)} satellites: {satellites}")
        
        sat_timestamps = {}
        for sat_id in satellites:
            mask = sat_ids == sat_id
            sat_timestamps[sat_id] = timestamps[mask]
            print(f"\nSatellite {sat_id:02d}: {len(sat_timestamps[sat_id])} samples")
            if len(sat_timestamps[sat_id]) > 0:
                print(f"  First timestamp: {sat_timestamps[sat_id][0]:.3f} seconds since 2000-01-01")
                print(f"  Last timestamp:  {sat_timestamps[sat_id][-1]:.3f} seconds since 2000-01-01")
        
        # Check alignment: compare timestamps at each index across satellites
        print("\n" + "="*60)
        print("Timestamp Alignment Check - MISALIGNED INDICES ONLY")
        print("="*60)
        
        # Find minimum length
        min_length = min(len(sat_timestamps[sat_id]) for sat_id in satellites)
        
        misaligned_indices = []
        misaligned_count = 0
        aligned_count = 0
        epoch_2000 = datetime(2000, 1, 1, 0, 0, 0)
        
        # Check all indices
        for i in range(min_length):
            timestamps_at_index = []
            for sat_id in satellites:
                ts = sat_timestamps[sat_id][i]
                timestamps_at_index.append(ts)
            
            # Check if all timestamps are the same (within 0.1 seconds tolerance)
            max_diff = max(timestamps_at_index) - min(timestamps_at_index)
            
            if max_diff < 0.1:
                aligned_count += 1
            else:
                misaligned_count += 1
                misaligned_indices.append((i, timestamps_at_index, max_diff))
        
        # Report only misaligned indices
        if misaligned_indices:
            print(f"\nFound {len(misaligned_indices)} misaligned indices out of {min_length} total:")
            print()
            for i, timestamps_at_index, max_diff in misaligned_indices:
                # Convert to datetime for readability (use first timestamp as reference)
                dt = epoch_2000 + timedelta(seconds=float(timestamps_at_index[0]))
                
                print(f"  Index {i:3d}: ✗ MISALIGNED (max diff: {max_diff:.3f}s)")
                print(f"    Expected time: {dt.strftime('%Y-%m-%d %H:%M:%S')} UTC")
                print(f"    Timestamps:")
                for j, sat_id in enumerate(satellites):
                    ts = timestamps_at_index[j]
                    dt_sat = epoch_2000 + timedelta(seconds=float(ts))
                    diff_from_first = ts - timestamps_at_index[0]
                    print(f"      CYG{sat_id:02d}: {ts:.3f}s ({dt_sat.strftime('%H:%M:%S')} UTC, diff: {diff_from_first:+.3f}s)")
                print()
        else:
            print("\n✓ All timestamps are perfectly aligned!")
        
        print("="*60)
        print("Summary")
        print("="*60)
        print(f"Total samples checked: {min_length}")
        print(f"Aligned samples: {aligned_count}")
        print(f"Misaligned samples: {misaligned_count}")
        
        if misaligned_count == 0:
            print("\n✓ SUCCESS: All timestamps are aligned across satellites!")
        else:
            print(f"\n✗ WARNING: {misaligned_count} samples are misaligned")
            print("  See details above for each misaligned index")
        
        # Check time intervals between consecutive samples (should be ~300s if every 600th point)
        print("\n" + "="*60)
        print("Time Interval Between Consecutive Samples")
        print("(Should be ~300 seconds = 5 minutes if every 600th point is taken)")
        print("="*60)
        all_intervals = []
        for sat_id in satellites:
            if len(sat_timestamps[sat_id]) > 1:
                intervals = np.diff(sat_timestamps[sat_id])
                avg_interval = np.mean(intervals)
                std_interval = np.std(intervals)
                min_interval = np.min(intervals)
                max_interval = np.max(intervals)
                all_intervals.extend(intervals)
                print(f"CYG{sat_id:02d}: avg = {avg_interval:.3f}s, std = {std_interval:.3f}s, range = [{min_interval:.3f}, {max_interval:.3f}]s")
        
        if all_intervals:
            overall_avg = np.mean(all_intervals)
            overall_std = np.std(all_intervals)
            print(f"\nOverall average interval: {overall_avg:.3f}s (expected: ~300.0s for every 600th point)")
            print(f"Overall std: {overall_std:.3f}s")
            if overall_avg < 10:
                print("  ⚠ WARNING: Intervals are too small - sampling may not be applied correctly!")
            elif abs(overall_avg - 300.0) > 10:
                print(f"  ⚠ WARNING: Average interval ({overall_avg:.1f}s) differs from expected 300s")
        
        # Check average difference between timestamps at same index across satellites
        print("\n" + "="*60)
        print("Average Timestamp Difference Across Satellites at Same Index")
        print("(Should be ~0 seconds if timestamps are aligned)")
        print("="*60)
        cross_sat_diffs = []
        for i in range(min_length):
            timestamps_at_index = [sat_timestamps[sat_id][i] for sat_id in satellites]
            max_ts = max(timestamps_at_index)
            min_ts = min(timestamps_at_index)
            diff = max_ts - min_ts
            cross_sat_diffs.append(diff)
        
        if cross_sat_diffs:
            avg_cross_diff = np.mean(cross_sat_diffs)
            max_cross_diff = np.max(cross_sat_diffs)
            std_cross_diff = np.std(cross_sat_diffs)
            print(f"Average difference: {avg_cross_diff:.3f}s")
            print(f"Maximum difference: {max_cross_diff:.3f}s")
            print(f"Std deviation: {std_cross_diff:.3f}s")
            if avg_cross_diff < 0.1:
                print("  ✓ Timestamps are well-aligned across satellites")
            else:
                print(f"  ⚠ WARNING: Timestamps differ by {avg_cross_diff:.3f}s on average")
                print("    This suggests satellites have different start times or sampling rates")
        
        return misaligned_count == 0

if __name__ == '__main__':
    extracted_file = Path('CYGNSS_data/pos_data/cygnss_extracted_data.nc')
    
    if not extracted_file.exists():
        print(f"Error: {extracted_file} does not exist")
        print("Please extract data first using extract_cygnss_data.py")
        exit(1)
    
    test_timestamp_alignment(extracted_file)
