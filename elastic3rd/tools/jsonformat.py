#!python
#
#
import json

filename = 'PureEle'
ext = '.json'

with open(filename + ext, 'r') as f:
    data = json.load(f)

with open('{}.json'.format(filename + "_formated" + ext), 'w') as fp:    
    json.dump(data, fp, indent=4, sort_keys=False)