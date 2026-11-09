#!/usr/bin/env python3
"""Explicit, cached acquisition; no package-load or test-time downloads."""
import argparse, hashlib, json, pathlib, urllib.request
FILES = {
    'README.txt': 'c3678e5249f9b91acfd1b282a1f5314f088bc1727b8d5f2fa27b01eddd48d801',
    'Ti_combined_data_oct24_2023_ace_order_3_degree_20_rcut_6.0_cw_moreV.json': 'aa4c9f27bbb17709ce31ab0bd7786b9f363b458f6e0ed2cfcbfe9a7a57834e8b',
    'ti64_combined_data_Apr08_rcut_6.0_order_3_degree_15_custom_reg.json': '652218003294aa3a582ace53c86f8c90749e643e9027d4b99a06cadb50d1d277',
}
def acquire(dest):
    dest.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen('https://zenodo.org/api/records/15114121') as response:
        metadata = response.read()
    record = json.loads(metadata)
    if record['metadata']['license']['id'] != 'cc-by-4.0':
        raise ValueError('Unexpected archive license; inspect before acquisition')
    (dest/'record.json').write_bytes(metadata)
    for name, expected in FILES.items():
        target = dest/name
        if not target.exists():
            url = next(f['links']['self'] for f in record['files'] if f['key'] == name)
            temporary = target.with_suffix(target.suffix+'.part')
            try:
                with urllib.request.urlopen(url) as response, temporary.open('wb') as stream:
                    while block := response.read(1024*1024):
                        stream.write(block)
                if hashlib.sha256(temporary.read_bytes()).hexdigest() != expected:
                    raise ValueError('Checksum mismatch: '+name)
                temporary.replace(target)
            finally:
                temporary.unlink(missing_ok=True)
        if hashlib.sha256(target.read_bytes()).hexdigest() != expected:
            raise ValueError('Checksum mismatch in cached file: '+name)
        print(name, expected)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--destination', type=pathlib.Path, required=True)
    acquire(parser.parse_args().destination)
