# NextClone Test Data

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
