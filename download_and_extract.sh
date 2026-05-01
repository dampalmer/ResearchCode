#!/bin/bash

# Script to download CYGNSS data for a given date, extract it, and clean up raw data
# Usage: ./download_and_extract.sh YYYYMMDD
# Example: ./download_and_extract.sh 20250920

set -e  # Exit on error

if [ -z "$1" ]; then
    echo "Usage: $0 YYYYMMDD"
    echo "Example: $0 20250920"
    exit 1
fi

DATE_STR="$1"

# Validate date format
if ! [[ "$DATE_STR" =~ ^[0-9]{8}$ ]]; then
    echo "Error: Date must be in YYYYMMDD format"
    exit 1
fi

# Get script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RAW_DATA_DIR="${SCRIPT_DIR}/CYGNSS_data/raw_data"
POS_DATA_DIR="${SCRIPT_DIR}/CYGNSS_data/pos_data"
EXTRACT_SCRIPT="${POS_DATA_DIR}/extract_cygnss_data.py"

echo "=========================================="
echo "CYGNSS Data Download and Extraction"
echo "Date: ${DATE_STR}"
echo "=========================================="
echo

# Step 1: Download data
echo "Step 1: Downloading data for ${DATE_STR}..."
cd "${RAW_DATA_DIR}"
if [ ! -f "data_import.sh" ]; then
    echo "Error: data_import.sh not found in ${RAW_DATA_DIR}"
    exit 1
fi

chmod +x data_import.sh
./data_import.sh "${DATE_STR}"

if [ $? -ne 0 ]; then
    echo "Error: Download failed"
    exit 1
fi

# Count downloaded files
DOWNLOADED_COUNT=$(ls -1 cyg*.nc 2>/dev/null | wc -l)
echo "Downloaded ${DOWNLOADED_COUNT} files"
echo

# Step 2: Extract data
echo "Step 2: Extracting position data..."
if [ ! -f "${EXTRACT_SCRIPT}" ]; then
    echo "Error: extract_cygnss_data.py not found in ${POS_DATA_DIR}"
    exit 1
fi

cd "${POS_DATA_DIR}"
python3 extract_cygnss_data.py "${RAW_DATA_DIR}" --pattern "*s${DATE_STR}*.nc" -o "cygnss_extracted_data.nc"

if [ $? -ne 0 ]; then
    echo "Error: Extraction failed"
    exit 1
fi

# Check if extraction was successful
if [ ! -f "cygnss_extracted_data.nc" ]; then
    echo "Error: Extracted data file was not created"
    exit 1
fi

EXTRACTED_SIZE=$(du -h cygnss_extracted_data.nc | cut -f1)
echo "Extracted data file: ${EXTRACTED_SIZE}"
echo

# Step 3: Delete raw data files
echo "Step 3: Cleaning up raw data files..."
cd "${RAW_DATA_DIR}"
RAW_SIZE=$(du -sh . 2>/dev/null | cut -f1)
echo "Raw data size before cleanup: ${RAW_SIZE}"

# Delete only the files for this date
DELETED_COUNT=0
for file in cyg*.ddmi.s${DATE_STR}*.nc; do
    if [ -f "$file" ]; then
        rm -f "$file"
        DELETED_COUNT=$((DELETED_COUNT + 1))
    fi
done

echo "Deleted ${DELETED_COUNT} raw data files for ${DATE_STR}"

# Show remaining size
if [ -d "${RAW_DATA_DIR}" ]; then
    REMAINING_SIZE=$(du -sh "${RAW_DATA_DIR}" 2>/dev/null | cut -f1)
    echo "Remaining raw data size: ${REMAINING_SIZE}"
fi

echo
echo "=========================================="
echo "Process complete!"
echo "Extracted data: ${POS_DATA_DIR}/cygnss_extracted_data.nc"
echo "=========================================="
