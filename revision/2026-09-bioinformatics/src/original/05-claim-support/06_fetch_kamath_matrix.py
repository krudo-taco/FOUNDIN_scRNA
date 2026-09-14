#!/usr/bin/env python3
"""Download the public GEO count matrix and verify its header against supplied features/barcodes."""
import gzip
import hashlib
import json
import shutil
import time
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXT = HERE / 'external'
NAME = 'GSE178265_Homo_matrix.mtx.gz'
URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE178nnn/GSE178265/suppl/' + NAME
path = EXT / NAME
part = EXT / (NAME + '.partial')

if not path.exists():
    start = part.stat().st_size if part.exists() else 0
    request = urllib.request.Request(URL, headers={'Range': 'bytes=%d-' % start} if start else {})
    with urllib.request.urlopen(request, timeout=60) as response:
        status = response.status
        if start and status == 206:
            assert response.headers['Content-Range'].startswith('bytes %d-' % start)
            mode = 'ab'
        else:
            assert status == 200
            mode, start = 'wb', 0
        expected = start + int(response.headers['Content-Length'])
        assert shutil.disk_usage(EXT).free > expected - start + 2 * 1024 ** 3
        total, last = start, time.time()
        print('Expected bytes', expected, 'resuming at', start, flush=True)
        with part.open(mode) as output:
            while True:
                block = response.read(4 * 1024 ** 2)
                if not block:
                    break
                output.write(block)
                total += len(block)
                if time.time() - last > 20:
                    print('Downloaded GiB', round(total / 1024 ** 3, 3), 'of', round(expected / 1024 ** 3, 3), flush=True)
                    last = time.time()
        assert total == expected, (total, expected)
    part.rename(path)

with gzip.open(path, 'rt') as f:
    header = f.readline().strip()
    assert header.lower() == '%%matrixmarket matrix coordinate integer general', header
    line = f.readline()
    while line.startswith('%'):
        line = f.readline()
    genes, cells, entries = map(int, line.split())
with gzip.open(EXT / 'GSE178265_Homo_features.tsv.gz', 'rt') as f:
    n_features = sum(1 for _ in f)
with gzip.open(EXT / 'GSE178265_Homo_bcd.tsv.gz', 'rt') as f:
    n_barcodes = sum(1 for _ in f)
assert (genes, cells) == (n_features, n_barcodes)
h = hashlib.sha256()
with path.open('rb') as f:
    for block in iter(lambda: f.read(8 * 1024 ** 2), b''):
        h.update(block)
audit = {'source': URL, 'bytes': path.stat().st_size, 'sha256': h.hexdigest(),
         'matrix_market_header': header, 'matrix_genes': genes, 'matrix_cells': cells,
         'matrix_nonzero_entries': entries, 'dimensions_match_feature_and_barcode_files': True,
         'full_stream_and_gzip_CRC_verification': 'Outside the download check; see kamath_count_extraction_audit.json for full-stream verification',
         'cell_annotations_and_disease_effects': 'Not verified by this download'}
(HERE / 'inputs/kamath_matrix_download_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
print(json.dumps(audit, indent=2), flush=True)
