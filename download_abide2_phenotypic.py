"""
Download ABIDE2 Composite Phenotypic Data
Fetches phenotypic CSV file from ABIDE project repository
"""

import os
import sys
import urllib.request
import urllib.error
from pathlib import Path
import json

# ============================================================================
# CONFIG
# ============================================================================
SCRIPT_DIR = Path(__file__).parent
OUTPUT_DIR = SCRIPT_DIR / 'ABIDE_2/raw'
OUTPUT_FILE = OUTPUT_DIR / 'ABIDEII_Composite_Phenotypic.csv'

# Download sources (in priority order)
DOWNLOAD_SOURCES = [
    {
        'name': 'AWS S3 ABIDE Data',
        'url': 'https://fcp-indi.s3.amazonaws.com/data/Projects/ABIDE/Phenotypic_V1.0b.csv',
        'description': 'AWS S3 repository with ABIDE phenotypic data'
    },
    {
        'name': 'GitHub ABIDE Repository',
        'url': 'https://raw.githubusercontent.com/ChildMindInstitute/ABIDE/master/ABIDEII_Composite_Phenotypic.csv',
        'description': 'GitHub repository mirror'
    },
    {
        'name': 'NITRC FCP-INDI Direct',
        'url': 'http://fcon_1000.projects.nitrc.org/indi/abide/ABIDEII_Composite_Phenotypic.csv',
        'description': 'Original NITRC source'
    },
]

# ============================================================================
# DOWNLOAD FUNCTIONS
# ============================================================================

def download_file(url, output_path, source_name):
    """Download file from URL with progress tracking."""
    print(f"\n  Attempting: {source_name}")
    print(f"  URL: {url}")
    
    try:
        # Create custom headers to avoid 403 errors
        headers = {
            'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
        }
        
        # Download with progress
        print(f"  Downloading...", end='', flush=True)
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req) as response:
            with open(output_path, 'wb') as out_file:
                out_file.write(response.read())
        
        file_size = os.path.getsize(output_path)
        file_size_mb = file_size / (1024 * 1024)
        
        print(f" ✓ Complete!")
        print(f"  File size: {file_size_mb:.2f} MB")
        return True
        
    except urllib.error.HTTPError as e:
        print(f" ✗ HTTP Error {e.code}: {e.reason}")
        if os.path.exists(output_path):
            os.remove(output_path)
        return False
    except urllib.error.URLError as e:
        print(f" ✗ URL Error: {e.reason}")
        if os.path.exists(output_path):
            os.remove(output_path)
        return False
    except Exception as e:
        print(f" ✗ Error: {e}")
        if os.path.exists(output_path):
            os.remove(output_path)
        return False


def verify_csv(filepath):
    """Verify that downloaded file is a valid CSV."""
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            # Read first line to check if it's CSV
            header = f.readline()
            if not header or ',' not in header:
                print(f"  WARNING: File may not be valid CSV (no commas in header)")
                return False
            
            # Count rows
            row_count = sum(1 for _ in f)
            print(f"  ✓ Valid CSV with {row_count} data rows")
            return True
    except Exception as e:
        print(f"  ✗ Error verifying file: {e}")
        return False


def main():
    print("=" * 70)
    print("DOWNLOADING ABIDE2 PHENOTYPIC DATA")
    print("=" * 70)
    
    # Create output directory
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    
    # Check if file already exists
    if OUTPUT_FILE.exists():
        file_size_mb = OUTPUT_FILE.stat().st_size / (1024 * 1024)
        print(f"\n⚠ File already exists: {OUTPUT_FILE}")
        print(f"  Size: {file_size_mb:.2f} MB")
        
        response = input("  Overwrite? (y/n): ").strip().lower()
        if response != 'y':
            print("  Skipping download.")
            return True
        
        os.remove(OUTPUT_FILE)
        print("  File removed.")
    
    # Try each download source
    success = False
    for i, source in enumerate(DOWNLOAD_SOURCES, 1):
        print(f"\n[{i}/{len(DOWNLOAD_SOURCES)}] {source['name']}")
        print(f"    {source['description']}")
        
        if download_file(source['url'], OUTPUT_FILE, source['name']):
            # Verify the file
            if verify_csv(OUTPUT_FILE):
                success = True
                break
            else:
                os.remove(OUTPUT_FILE)
                print("  File verification failed, trying next source...")
    
    # Summary
    print("\n" + "=" * 70)
    if success:
        print("✓ SUCCESS! ABIDE2 phenotypic data downloaded.")
        print(f"  Saved to: {OUTPUT_FILE}")
        print("\nNext steps:")
        print("  1. Run: python phenotypic_eda.py")
        print("  2. Run: python comparability.py")
        print("  3. Run: python generate_report.py")
        return True
    else:
        print("✗ FAILED! Could not download ABIDE2 phenotypic data.")
        print("\nAlternative: Manual download")
        print("  1. Visit: http://fcon_1000.projects.nitrc.org/indi/abide/")
        print("  2. Download: ABIDEII_Composite_Phenotypic.csv")
        print(f"  3. Save to: {OUTPUT_FILE}")
        return False


if __name__ == '__main__':
    try:
        success = main()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        print("\n\nDownload cancelled by user.")
        sys.exit(1)
    except Exception as e:
        print(f"\nFatal error: {e}")
        sys.exit(1)
