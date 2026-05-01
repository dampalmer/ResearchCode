#!/usr/bin/env python3
"""
Create a small test NetCDF file with known proximity events for debugging.
"""

import numpy as np
import netCDF4
from pathlib import Path

# Create test data with known proximity events
# Format: timestamp, satellite_id, pos_x, pos_y, pos_z

# Test scenario:
# Time 0: Sat 1 at (0, 0, 0), Sat 2 at (50km, 0, 0) = 50km apart (within threshold)
# Time 0: Sat 3 at (0, 0, 0), Sat 4 at (150km, 0, 0) = 150km apart (outside threshold)
# Time 10: Sat 1 at (0, 0, 0), Sat 2 at (80km, 0, 0) = 80km apart (within threshold)
# Time 10: Sat 3 at (0, 0, 0), Sat 4 at (0, 90km, 0) = 90km apart (within threshold)
# Time 20: Sat 1 at (0, 0, 0), Sat 2 at (200km, 0, 0) = 200km apart (outside threshold)

test_data = [
    # Time 0.0 seconds
    (0.0, 1, 0.0, 0.0, 0.0),           # Sat 1 at origin
    (0.0, 2, 50000.0, 0.0, 0.0),       # Sat 2 at 50km (within 100km threshold) ✓
    (0.0, 3, 0.0, 0.0, 0.0),           # Sat 3 at origin
    (0.0, 4, 150000.0, 0.0, 0.0),      # Sat 4 at 150km (outside threshold)
    
    # Time 10.0 seconds
    (10.0, 1, 0.0, 0.0, 0.0),          # Sat 1 at origin
    (10.0, 2, 80000.0, 0.0, 0.0),      # Sat 2 at 80km (within threshold) ✓
    (10.0, 3, 0.0, 0.0, 0.0),          # Sat 3 at origin
    (10.0, 4, 0.0, 90000.0, 0.0),      # Sat 4 at 90km (within threshold) ✓
    
    # Time 20.0 seconds
    (20.0, 1, 0.0, 0.0, 0.0),          # Sat 1 at origin
    (20.0, 2, 200000.0, 0.0, 0.0),    # Sat 2 at 200km (outside threshold)
    (20.0, 3, 0.0, 0.0, 0.0),          # Sat 3 at origin
    (20.0, 4, 0.0, 0.0, 95000.0),      # Sat 4 at 95km (within threshold) ✓
]

# Expected proximity events:
# Time 0.0: Sat 1 - Sat 2 (50km) ✓
# Time 0.0: Sat 1 - Sat 3 (0km, same position) ✓
# Time 10.0: Sat 1 - Sat 2 (80km) ✓
# Time 10.0: Sat 1 - Sat 3 (0km, same position) ✓
# Time 10.0: Sat 3 - Sat 4 (90km) ✓
# Time 20.0: Sat 1 - Sat 3 (0km, same position) ✓
# Time 20.0: Sat 3 - Sat 4 (95km) ✓

output_file = Path('test_proximity_data.nc')

# Create NetCDF file
with netCDF4.Dataset(output_file, 'w', format='NETCDF4') as nc:
    n_samples = len(test_data)
    
    # Create dimensions
    nc.createDimension('sample', n_samples)
    
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
    
    # Extract data
    timestamps = [d[0] for d in test_data]
    sat_ids = [d[1] for d in test_data]
    pos_x = [d[2] for d in test_data]
    pos_y = [d[3] for d in test_data]
    pos_z = [d[4] for d in test_data]
    
    # Assign data
    timestamp_var[:] = timestamps
    pos_x_var[:] = pos_x
    pos_y_var[:] = pos_y
    pos_z_var[:] = pos_z
    sat_id_var[:] = sat_ids
    
    # Add global attributes
    nc.title = 'Test CYGNSS Data with Known Proximity Events'
    nc.description = 'Small test dataset for debugging proximity detection'

print(f"Created test data file: {output_file}")
print(f"\nExpected Proximity Events (within 100km):")
print("="*60)
print("Time 0.0s:  Sat 1 - Sat 2 (50km)")
print("Time 0.0s:  Sat 1 - Sat 3 (0km - same position)")
print("Time 10.0s: Sat 1 - Sat 2 (80km)")
print("Time 10.0s: Sat 1 - Sat 3 (0km - same position)")
print("Time 10.0s: Sat 3 - Sat 4 (90km)")
print("Time 20.0s: Sat 1 - Sat 3 (0km - same position)")
print("Time 20.0s: Sat 3 - Sat 4 (95km)")
print("="*60)
print(f"\nTotal expected events: 7")
