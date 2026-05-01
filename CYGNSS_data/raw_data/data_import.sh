#!/bin/bash

# Modified data_import.sh to accept a date parameter
# Usage: ./data_import.sh YYYYMMDD
# Example: ./data_import.sh 20250920

GREP_OPTIONS=''

# Get date from command line argument or use default
if [ -z "$1" ]; then
    echo "Usage: $0 YYYYMMDD"
    echo "Example: $0 20250920"
    exit 1
fi

DATE_STR="$1"
# Validate date format (YYYYMMDD)
if ! [[ "$DATE_STR" =~ ^[0-9]{8}$ ]]; then
    echo "Error: Date must be in YYYYMMDD format"
    exit 1
fi

# Format: sYYYYMMDD-000000-eYYYYMMDD-235959
DATE_PATTERN="s${DATE_STR}-000000-e${DATE_STR}-235959"

cookiejar=$(mktemp cookies.XXXXXXXXXX)
netrc=$(mktemp netrc.XXXXXXXXXX)
chmod 0600 "$cookiejar" "$netrc"
function finish {
  rm -rf "$cookiejar" "$netrc"
}

trap finish EXIT
WGETRC="$wgetrc"

check_existing_netrc() {
    # Check if ~/.netrc exists and has Earthdata credentials
    home_netrc="$HOME/.netrc"
    if [ -f "$home_netrc" ] && grep -q "urs.earthdata.nasa.gov" "$home_netrc" 2>/dev/null; then
        # Copy Earthdata credentials to temp netrc file for testing
        grep "urs.earthdata.nasa.gov" "$home_netrc" >> "$netrc" 2>/dev/null
        
        # Test if the credentials work (match the same validation logic as detect_app_approval)
        test_url="https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg01.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc"
        test_status=$(curl -s -b "$cookiejar" -c "$cookiejar" -L --max-redirs 5 --netrc-file "$netrc" "$test_url" -w '\n%{http_code}' 2>/dev/null | tail -1)
        
        # Only accept 200, 301, or 302 (same as detect_app_approval function)
        if [ "$test_status" = "200" ] || [ "$test_status" = "301" ] || [ "$test_status" = "302" ]; then
            echo "Using existing Earthdata credentials from ~/.netrc"
            return 0
        else
            # Credentials don't work - remove from temp netrc and return failure
            # Recreate temp netrc without the Earthdata line
            grep -v "urs.earthdata.nasa.gov" "$netrc" > "${netrc}.tmp" 2>/dev/null || true
            mv "${netrc}.tmp" "$netrc" 2>/dev/null || touch "$netrc"
            return 1
        fi
    fi
    return 1
}

prompt_credentials() {
    echo "Enter your Earthdata Login or other provider supplied credentials"
    read -p "Username (dmpalmer): " username
    username=${username:-dmpalmer}
    read -s -p "Password: " password
    echo
    echo "machine urs.earthdata.nasa.gov login $username password $password" >> $netrc
    
    # Ask if user wants to save credentials to ~/.netrc for future use
    read -p "Save credentials to ~/.netrc for future use? (y/n): " save_creds
    if [ "$save_creds" = "y" ] || [ "$save_creds" = "Y" ]; then
        home_netrc="$HOME/.netrc"
        # Remove existing Earthdata entry if present
        if [ -f "$home_netrc" ]; then
            grep -v "urs.earthdata.nasa.gov" "$home_netrc" > "${home_netrc}.tmp" 2>/dev/null || true
            mv "${home_netrc}.tmp" "$home_netrc" 2>/dev/null || true
        fi
        # Add new entry
        echo "machine urs.earthdata.nasa.gov login $username password $password" >> "$home_netrc"
        chmod 0600 "$home_netrc"
        echo "Credentials saved to ~/.netrc"
    fi
    echo
}

exit_with_error() {
    echo
    echo "Unable to Retrieve Data"
    echo
    echo $1
    echo
    echo "Example URL: https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg01.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc"
    echo
    exit 1
}

# Try to use existing credentials first, otherwise prompt
if ! check_existing_netrc; then
    prompt_credentials
fi
  detect_app_approval() {
    # Use satellite 01 for testing
    test_url="https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg01.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc"
    approved=`curl -s -b "$cookiejar" -c "$cookiejar" -L --max-redirs 5 --netrc-file "$netrc" "$test_url" -w '\n%{http_code}' | tail  -1`
    if [ "$approved" -ne "200" ] && [ "$approved" -ne "301" ] && [ "$approved" -ne "302" ]; then
        # User didn't approve the app. Direct users to approve the app in URS
        exit_with_error "Please ensure that you have authorized the remote application by visiting the link below "
    fi
}

setup_auth_curl() {
    # Firstly, check if it require URS authentication
    test_url="https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg01.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc"
    status=$(curl -s -z "$(date)" -w '\n%{http_code}' "$test_url" | tail -1)
    if [[ "$status" -ne "200" && "$status" -ne "304" ]]; then
        # URS authentication is required. Now further check if the application/remote service is approved.
        detect_app_approval
    fi
}

setup_auth_wget() {
    # The safest way to auth via curl is netrc. Note: there's no checking or feedback
    # if login is unsuccessful
    touch ~/.netrc
    chmod 0600 ~/.netrc
    credentials=$(grep 'machine urs.earthdata.nasa.gov' ~/.netrc)
    if [ -z "$credentials" ]; then
        cat "$netrc" >> ~/.netrc
    fi
}

fetch_urls() {
  if command -v curl >/dev/null 2>&1; then
      setup_auth_curl
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        # Download with error handling - continue on 404, exit on other errors
        http_code=$(curl -f -b "$cookiejar" -c "$cookiejar" -L --netrc-file "$netrc" -g -o "$stripped_query_params" -w '%{http_code}' -- "$line" 2>&1 | tail -1)
        if [ "$http_code" = "200" ] || [ "$http_code" = "301" ] || [ "$http_code" = "302" ]; then
            echo "Downloaded: $stripped_query_params"
        elif [ "$http_code" = "404" ]; then
            echo "Warning: File not found (404) - $stripped_query_params (satellite may not have data for this date)"
            # Remove empty file if created
            [ -f "$stripped_query_params" ] && rm -f "$stripped_query_params"
        else
            echo "Error: Failed to download $stripped_query_params (HTTP $http_code)"
            [ -f "$stripped_query_params" ] && rm -f "$stripped_query_params"
            # Don't exit on individual file failures, but note the error
        fi
      done;
  elif command -v wget >/dev/null 2>&1; then
      # We can't use wget to poke provider server to get info whether or not URS was integrated without download at least one of the files.
      echo
      echo "WARNING: Can't find curl, use wget instead."
      echo "WARNING: Script may not correctly identify Earthdata Login integrations."
      echo
      setup_auth_wget
      while read -r line; do
        # Get everything after the last '/'
        filename="${line##*/}"

        # Strip everything after '?'
        stripped_query_params="${filename%%\?*}"

        # Download with error handling - continue on 404, note other errors
        if wget --load-cookies "$cookiejar" --save-cookies "$cookiejar" --output-document "$stripped_query_params" --keep-session-cookies -- "$line" 2>&1 | grep -q "404"; then
            echo "Warning: File not found (404) - $stripped_query_params (satellite may not have data for this date)"
            [ -f "$stripped_query_params" ] && rm -f "$stripped_query_params"
        elif wget --load-cookies "$cookiejar" --save-cookies "$cookiejar" --output-document "$stripped_query_params" --keep-session-cookies -- "$line" 2>&1; then
            echo "Downloaded: $stripped_query_params"
        else
            echo "Error: Failed to download $stripped_query_params"
            [ -f "$stripped_query_params" ] && rm -f "$stripped_query_params"
        fi
      done;
  else
      exit_with_error "Error: Could not find a command-line downloader.  Please install curl or wget"
  fi
}

# Generate URLs for all satellites (01-08)
echo "Generating URLs for date: ${DATE_STR}"
echo "Satellites: 01-08"

# Use unquoted heredoc delimiter to allow variable expansion
fetch_urls <<EDSCEOF
https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg01.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc
https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg02.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc
https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg03.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc
https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg04.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc
https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg05.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc
https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg06.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc
https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg07.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc
https://archive.podaac.earthdata.nasa.gov/podaac-ops-cumulus-protected/CYGNSS_L1_V3.2/cyg08.ddmi.${DATE_PATTERN}.l1.power-brcs.a32.d33.nc
EDSCEOF
