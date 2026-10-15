#!python
import os
from monty.serialization import loadfn, dumpfn

lat = 'bcc'

dirs = os.listdir('.')
results = {}
for diri in dirs:
    if os.path.isdir(diri):
        fullname = os.path.join(diri, '{}-{}-elastic.json'.format(lat, diri))
        data = loadfn(fullname)
        results.update({diri: data})
dumpfn(results, '{}-elastic.json'.format(lat))
