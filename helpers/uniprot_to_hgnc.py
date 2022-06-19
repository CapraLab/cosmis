#!/usr/bin/env python3

import sys
import json
import urllib.parse
import urllib.request


url = 'https://www.uniprot.org/uploadlists/'

with open(sys.argv[1], 'rt') as ipf:
    queries = [l.strip() for l in ipf]

params = {
    'from': 'ACC+ID',
    'to': 'GENENAME',
    'format': 'tab',
    'query': ' '.join(queries)
}

data = urllib.parse.urlencode(params)
data = data.encode('utf-8')
req = urllib.request.Request(url, data)
with urllib.request.urlopen(req) as f:
   response = f.read()

response_str = response.decode('utf-8').split('\n')
# remove the header line
uniprot_to_hgnc = {}
for line in response_str[1:]:
    if line.strip():
        u, e = line.split('\t')
        if u not in uniprot_to_hgnc:
            uniprot_to_hgnc[u] = e

# print the total number of mapped uniprot ids
print(
    '{} out of {} UniProt IDs mapped.'.format(
        len(uniprot_to_hgnc.keys()), len(queries)
    )
)
print(
    'UniProt IDs not mapped:\n{}'.format(
        queries - uniprot_to_hgnc.keys()
    )
)

# write to disk file in json format
with open(sys.argv[2], 'wt') as opf:
    json.dump(uniprot_to_hgnc, opf, indent=4)
