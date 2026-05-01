#!/usr/bin/env python3
"""
CYGNSS Proximity Events Dashboard

Displays proximity events from the proximity events database.
The database is created and maintained by process_proximity_events.py
"""

import numpy as np
import netCDF4
import dash
from dash import dcc, html, dash_table
import pandas as pd
from pathlib import Path
from datetime import datetime, timedelta

# Constants
PROXIMITY_DB = Path('CYGNSS_data/proximity_events.nc')

def load_proximity_events():
    """Load proximity events from the database NetCDF file.
    Returns (dataframe, threshold_km)."""
    if not PROXIMITY_DB.exists():
        return pd.DataFrame(columns=['satellite_1', 'satellite_2', 'start_time', 'end_time', 'min_distance_km']), 100.0
    
    print(f"Loading proximity events from {PROXIMITY_DB}...")
    
    with netCDF4.Dataset(PROXIMITY_DB, 'r') as nc:
        # Check if database has any events
        if 'event' not in nc.dimensions or nc.dimensions['event'].size == 0:
            print("Database is empty")
            threshold = getattr(nc, 'proximity_threshold_km', 100.0)
            return pd.DataFrame(columns=['satellite_1', 'satellite_2', 'start_time', 'end_time', 'min_distance_km']), threshold
        
        # Read all events
        sat1 = nc.variables['satellite_1'][:]
        sat2 = nc.variables['satellite_2'][:]
        start_time = nc.variables['start_time'][:]
        end_time = nc.variables['end_time'][:]
        min_distance = nc.variables['min_distance_km'][:]
        
        # Get threshold from attributes if available
        threshold = getattr(nc, 'proximity_threshold_km', 100.0)
        total_events = getattr(nc, 'total_events', len(sat1))
    
    df = pd.DataFrame({
        'satellite_1': sat1,
        'satellite_2': sat2,
        'start_time': start_time,
        'end_time': end_time,
        'min_distance_km': min_distance
    })
    
    print(f"Loaded {len(df)} proximity events (threshold: {threshold} km)")
    return df, threshold

# Load proximity events
print("="*60)
print("CYGNSS Proximity Events Dashboard")
print("="*60)

proximity_df, threshold_km = load_proximity_events()

# Create Dash app
app = dash.Dash(__name__)
app.title = "CYGNSS Proximity Events"

# Prepare table data
epoch = datetime(2000, 1, 1, 0, 0, 0)

if len(proximity_df) > 0:
    display_df = proximity_df.copy()
    display_df['Satellites'] = display_df.apply(
        lambda row: f"CYG{int(row['satellite_1']):02d} - CYG{int(row['satellite_2']):02d}", 
        axis=1
    )
    
    # Convert seconds to UTC datetime strings
    display_df['Start Time (UTC)'] = display_df['start_time'].apply(
        lambda s: (epoch + timedelta(seconds=float(s))).strftime('%Y-%m-%d %H:%M:%S UTC')
    )
    display_df['End Time (UTC)'] = display_df['end_time'].apply(
        lambda s: (epoch + timedelta(seconds=float(s))).strftime('%Y-%m-%d %H:%M:%S UTC')
    )
    display_df['Min Distance (km)'] = display_df['min_distance_km'].round(2)
    
    table_data = display_df[['Start Time (UTC)', 'End Time (UTC)', 'Satellites', 'Min Distance (km)']].to_dict('records')
    table_columns = [{"name": col, "id": col} for col in ['Start Time (UTC)', 'End Time (UTC)', 'Satellites', 'Min Distance (km)']]
else:
    table_data = []
    table_columns = []

app.layout = html.Div([
    html.H1("CYGNSS Satellite Proximity Events", 
            style={'textAlign': 'center', 'marginBottom': '20px'}),
    
    html.Div([
        html.P(f"Proximity Threshold: {threshold_km} km"),
        html.P(f"Total Events in Database: {len(proximity_df)}"),
        html.P(f"Database Location: {PROXIMITY_DB}"),
        html.P([
            "Note: To add new events, run ",
            html.Code("process_proximity_events.py YYYYMMDD"),
            " to process data for specific dates."
        ], style={'fontSize': '12px', 'color': '#666', 'marginTop': '10px'})
    ], style={'padding': '20px', 'backgroundColor': '#f0f0f0', 'borderRadius': '5px', 'marginBottom': '20px'}),
    
    html.Div([
        dash_table.DataTable(
            data=table_data,
            columns=table_columns,
            style_table={'overflowX': 'auto'},
            style_cell={'textAlign': 'left', 'padding': '10px'},
            style_header={'backgroundColor': 'rgb(230, 230, 230)', 'fontWeight': 'bold'},
            page_size=50,
            sort_action='native',
            filter_action='native',
            export_format='csv',
            export_headers='display'
        ) if len(table_data) > 0 else html.Div([
            html.P("No proximity events found in database."),
            html.P("Run process_proximity_events.py to process data and populate the database.")
        ], style={'padding': '20px', 'textAlign': 'center'})
    ]),
])

if __name__ == '__main__':
    print(f"\nStarting server on http://127.0.0.1:8050")
    print("="*60 + "\n")
    app.run(debug=True, host='127.0.0.1', port=8050)
