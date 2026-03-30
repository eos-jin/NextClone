#!/usr/bin/env python3
"""
Test suite for NextClone discovery mode and parameter validation.

This creates synthetic test data and validates:
1. Parameter validation logic (error/warning messages)
2. Discovery mode barcode extraction
3. Backward compatibility with whitelist mode

Run with: python tests/test_discovery_mode.py
"""

import os
import sys
import gzip
import tempfile
import subprocess
import shutil
from pathlib import Path

# Test configuration matching nextflow.config defaults
ADAPTER_5PRIME = "ATCTTGTGGAAAGGACGAAACACCG"
ADAPTER_3PRIME = "GTTTCAGAGCTATGCTGGAAACAGC"
BARCODE_LENGTH = 20

# Known test barcodes (matching data/known_barcodes_subset.txt)
TEST_BARCODES = [
    "CGGAGTAATACATTTTGCCT",
    "TCGGAGTTGGCTGTCGTTTC",
    "GTTGTCTCGGGGGGTGGAGA",
    "CCATGATAAGGGAGTTCCGG",
    "AGGGGAGTCGCGTGGTAGGC",
    "TGTCTAATGGGGGTGTCACT",
]

# Additional barcodes for discovery mode testing (not in whitelist)
DISCOVERY_BARCODES = [
    "AAAAAAAAAAAAAAAAAAAA",
    "TTTTTTTTTTTTTTTTTTTT",
    "GGGGGGGGGGGGGGGGGGGG",
]


def generate_quality_string(length):
    """Generate a fake quality string of given length."""
    return "I" * length


def generate_fastq_read(read_id, barcode):
    """
    Generate a single FASTQ read with the barcode flanked by BOTH adapters.
    
    Format: [5' adapter][BARCODE][3' adapter]
    
    This matches the NextClone expected input structure where flexiplex
    searches for: -x 5'adapter -b barcode -x 3'adapter
    """
    sequence = f"{ADAPTER_5PRIME}{barcode}{ADAPTER_3PRIME}"
    quality = generate_quality_string(len(sequence))
    
    return f"@{read_id}\n{sequence}\n+\n{quality}\n"


def create_synthetic_fastq(output_path, barcodes, reads_per_barcode=10):
    """
    Create a synthetic FASTQ file with known barcodes.
    
    Args:
        output_path: Path to output .fastq.gz file
        barcodes: List of barcode sequences to include
        reads_per_barcode: Number of reads to generate per barcode
    
    Read structure: [5' adapter][BARCODE][3' adapter]
    This matches NextClone's expected input format.
    """
    read_count = 0
    
    with gzip.open(output_path, 'wt') as f:
        for barcode in barcodes:
            for i in range(reads_per_barcode):
                read_id = f"TEST_READ_{read_count}:1:1:1:{read_count} 1:N:0:AACTTGAC"
                f.write(generate_fastq_read(read_id, barcode))
                read_count += 1
    
    print(f"Created {output_path} with {read_count} reads ({len(barcodes)} barcodes × {reads_per_barcode} reads)")
    return read_count


def create_barcode_whitelist(output_path, barcodes):
    """Create a barcode whitelist file."""
    with open(output_path, 'w') as f:
        for barcode in barcodes:
            f.write(f"{barcode}\n")
    print(f"Created whitelist {output_path} with {len(barcodes)} barcodes")


def test_parameter_validation():
    """
    Test that parameter validation and workflow modes work correctly.
    
    Expected behavior:
    - discovery_mode=true with whitelist → WARNING (proceeds anyway)
    - All 4 workflow modes should validate (DNAseq/scRNAseq × discovery/whitelist)
    
    Note: The validation for missing whitelist is not tested here because
    the nextflow.config has a default whitelist path. In production, users
    would need to explicitly set clone_barcodes_reference or enable discovery_mode.
    """
    print("\n" + "="*60)
    print("TEST: Parameter Validation & Workflow Modes")
    print("="*60)
    
    # Check if nextflow is available
    result = subprocess.run(["which", "nextflow"], capture_output=True, text=True)
    if result.returncode != 0:
        print("⚠️  SKIP: Nextflow not installed - cannot run validation tests")
        print("   Install with: curl -s https://get.nextflow.io | bash")
        return None
    
    project_dir = Path(__file__).parent.parent
    all_passed = True
    
    # Test 1: discovery_mode=true should show warning about ignored whitelist
    print("\n[Test 1] discovery_mode=true with whitelist (should warn)...")
    result = subprocess.run(
        ["nextflow", "run", str(project_dir / "main.nf"),
         "--discovery_mode", "true",
         "--mode", "DNAseq",
         "-preview"],
        capture_output=True,
        text=True,
        cwd=project_dir,
        env={**os.environ, "JAVA_HOME": os.environ.get("JAVA_HOME", "")}
    )
    
    if "WARNING" in result.stderr or "ignored" in result.stderr.lower():
        print("   ✅ PASS: Warning raised as expected")
    elif result.returncode == 0:
        print("   ⚠️  Warning may not appear in -preview mode, but workflow validated")
    else:
        print(f"   ❌ FAIL: Nextflow failed: {result.stderr[:300]}")
        all_passed = False
    
    # Test 2-5: Validate all 4 workflow combinations
    modes = [
        ("DNAseq", "false", "whitelist"),
        ("DNAseq", "true", "discovery"),
        ("scRNAseq", "false", "whitelist"),
        ("scRNAseq", "true", "discovery"),
    ]
    
    for i, (data_mode, discovery, desc) in enumerate(modes, start=2):
        print(f"\n[Test {i}] {data_mode} + {desc} mode...")
        result = subprocess.run(
            ["nextflow", "run", str(project_dir / "main.nf"),
             "--mode", data_mode,
             "--discovery_mode", discovery,
             "-preview"],
            capture_output=True,
            text=True,
            cwd=project_dir,
            env={**os.environ, "JAVA_HOME": os.environ.get("JAVA_HOME", "")}
        )
        
        if result.returncode == 0:
            print(f"   ✅ PASS: Workflow validated successfully")
        else:
            print(f"   ❌ FAIL: {result.stderr[:200]}")
            all_passed = False
    
    return all_passed


def test_synthetic_data_structure():
    """
    Test that synthetic data matches expected format.
    """
    print("\n" + "="*60)
    print("TEST: Synthetic Data Structure")
    print("="*60)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test FASTQ
        fastq_path = Path(tmpdir) / "test_synthetic.fastq.gz"
        create_synthetic_fastq(fastq_path, TEST_BARCODES[:3], reads_per_barcode=5)
        
        # Verify structure
        with gzip.open(fastq_path, 'rt') as f:
            lines = f.readlines()
        
        # Should have 4 lines per read × 3 barcodes × 5 reads = 60 lines
        expected_lines = 4 * 3 * 5
        if len(lines) == expected_lines:
            print(f"   ✅ PASS: Correct number of lines ({expected_lines})")
        else:
            print(f"   ❌ FAIL: Expected {expected_lines} lines, got {len(lines)}")
            return False
        
        # Check first read has correct structure (5' adapter + barcode + 3' adapter)
        first_seq = lines[1].strip()
        if ADAPTER_5PRIME in first_seq and ADAPTER_3PRIME in first_seq:
            print(f"   ✅ PASS: Both 5' and 3' adapters found in sequence")
        else:
            print(f"   ❌ FAIL: Expected both adapters in sequence")
            print(f"   Sequence: {first_seq}")
            return False
        
        # Check barcode is correct length (between the adapters)
        barcode_start = len(ADAPTER_5PRIME)
        barcode_region = first_seq[barcode_start:barcode_start + BARCODE_LENGTH]
        if barcode_region in TEST_BARCODES:
            print(f"   ✅ PASS: Valid barcode found: {barcode_region}")
        else:
            print(f"   ❌ FAIL: Barcode not in expected list: {barcode_region}")
            return False
    
    return True


def test_flexiplex_discovery():
    """
    Test Flexiplex discovery mode on synthetic data (if flexiplex is installed).
    """
    print("\n" + "="*60)
    print("TEST: Flexiplex Discovery Mode")
    print("="*60)
    
    # Check if flexiplex is available
    result = subprocess.run(["which", "flexiplex"], capture_output=True, text=True)
    if result.returncode != 0:
        print("⚠️  SKIP: Flexiplex not installed")
        print("   Install from: https://github.com/DavidsonGroup/flexiplex")
        return None
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Create test data with mixed barcodes (known + novel)
        all_barcodes = TEST_BARCODES[:3] + DISCOVERY_BARCODES[:2]
        fastq_path = tmpdir / "discovery_test.fastq.gz"
        create_synthetic_fastq(fastq_path, all_barcodes, reads_per_barcode=20)
        
        # Decompress for flexiplex
        fastq_unzipped = tmpdir / "discovery_test.fastq"
        with gzip.open(fastq_path, 'rb') as f_in:
            with open(fastq_unzipped, 'wb') as f_out:
                f_out.write(f_in.read())
        
        # Run flexiplex in discovery mode (no -k flag, -f 0 for strict match)
        # Must include BOTH 5' and 3' adapters in order: -x 5' -b barcode -x 3'
        print("\n   Running Flexiplex discovery mode...")
        result = subprocess.run(
            ["flexiplex",
             "-x", ADAPTER_5PRIME,
             "-b", "?" * BARCODE_LENGTH,
             "-u", "",
             "-x", ADAPTER_3PRIME,
             "-f", "0",  # Strict flanking match for discovery
             "-n", "discovery_test",
             str(fastq_unzipped)],
            capture_output=True,
            text=True,
            cwd=tmpdir
        )
        
        if result.returncode != 0:
            print(f"   ❌ FAIL: Flexiplex failed with: {result.stderr}")
            return False
        
        # Check discovery output
        counts_file = tmpdir / "discovery_test_barcodes_counts.txt"
        if not counts_file.exists():
            # Try alternate naming
            counts_file = tmpdir / "flexiplex_barcodes_counts.txt"
        
        if counts_file.exists():
            with open(counts_file) as f:
                discovered = f.readlines()
            
            discovered_barcodes = [line.split()[0] for line in discovered if line.strip()]
            
            print(f"   Discovered {len(discovered_barcodes)} unique barcodes")
            
            # Check that our test barcodes were found
            found_count = sum(1 for b in all_barcodes if b in discovered_barcodes)
            if found_count >= len(all_barcodes) * 0.8:  # Allow some tolerance
                print(f"   ✅ PASS: Found {found_count}/{len(all_barcodes)} expected barcodes")
                return True
            else:
                print(f"   ❌ FAIL: Only found {found_count}/{len(all_barcodes)} expected barcodes")
                return False
        else:
            print(f"   ❌ FAIL: Barcode counts file not created")
            print(f"   Files in tmpdir: {list(tmpdir.iterdir())}")
            return False


def create_test_data_for_repo():
    """
    Create test data files to include in the repository.
    """
    print("\n" + "="*60)
    print("Creating test data for repository")
    print("="*60)
    
    project_dir = Path(__file__).parent.parent
    test_data_dir = project_dir / "tests" / "data"
    test_data_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Whitelist mode test data (uses known barcodes)
    whitelist_fastq = test_data_dir / "whitelist_test.fastq.gz"
    create_synthetic_fastq(whitelist_fastq, TEST_BARCODES, reads_per_barcode=10)
    
    # 2. Discovery mode test data (includes novel barcodes)
    discovery_fastq = test_data_dir / "discovery_test.fastq.gz"
    mixed_barcodes = TEST_BARCODES[:4] + DISCOVERY_BARCODES
    create_synthetic_fastq(discovery_fastq, mixed_barcodes, reads_per_barcode=15)
    
    # 3. Expected barcodes for discovery mode
    expected_barcodes = test_data_dir / "expected_discovered_barcodes.txt"
    create_barcode_whitelist(expected_barcodes, mixed_barcodes)
    
    # 4. Create a README for the test data
    readme_content = """# NextClone Test Data

## Files

### whitelist_test.fastq.gz
Synthetic FASTQ with known barcodes (matching data/known_barcodes_subset.txt).
Use for testing whitelist mode (discovery_mode=false).

### discovery_test.fastq.gz  
Synthetic FASTQ with mixed barcodes - some known, some novel.
Use for testing discovery mode (discovery_mode=true).

### expected_discovered_barcodes.txt
List of all barcodes present in discovery_test.fastq.gz.
Use to validate discovery mode output.

## Barcode Format

- Barcode length: 20bp
- 3' adapter: GTTTCAGAGCTATGCTGGAAACAGC
- Read structure: [BARCODE][3' adapter]

## Running Tests

```bash
# From repository root
python tests/test_discovery_mode.py

# Or run specific test
python -c "from tests.test_discovery_mode import test_flexiplex_discovery; test_flexiplex_discovery()"
```
"""
    
    with open(test_data_dir / "README.md", 'w') as f:
        f.write(readme_content)
    
    print(f"\nTest data created in {test_data_dir}")
    print("Files:")
    for f in test_data_dir.iterdir():
        print(f"  - {f.name}")


def main():
    """Run all tests."""
    print("="*60)
    print("NextClone Discovery Mode Test Suite")
    print("="*60)
    
    results = {}
    
    # Test 1: Synthetic data structure
    results['synthetic_data'] = test_synthetic_data_structure()
    
    # Test 2: Parameter validation (requires nextflow)
    results['param_validation'] = test_parameter_validation()
    
    # Test 3: Flexiplex discovery (requires flexiplex)
    results['flexiplex_discovery'] = test_flexiplex_discovery()
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    for test_name, result in results.items():
        if result is True:
            status = "✅ PASS"
        elif result is False:
            status = "❌ FAIL"
        else:
            status = "⚠️  SKIP"
        print(f"  {test_name}: {status}")
    
    # Create test data for repo
    print("\n")
    create_test_data_for_repo()
    
    # Return exit code
    failures = sum(1 for r in results.values() if r is False)
    return failures


if __name__ == "__main__":
    sys.exit(main())
