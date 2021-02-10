#Organize experimental data into a databank
#Amount of each molecule in the membrane is given as a ratio. For a membrane of single molecule it's 1.

import os
import numpy as np

import yaml
import json

exp_information="""
#BEGIN

@TEMPERATURE=298
@LIPIDS=POPS
@AMOUNTS=1

#END
"""

save_dir = "/media/akiirikk/DATADRIVE1/tietokanta/expDATABANK/Data/POPS/T298K"

save_info = {}
lipids = []
amounts = []

for line in exp_information.split("\n"):
     if line.startswith('@'):
        key, value = line.split("=")

        if key.strip('@') == "LIPIDS":
           if "," in key:
               lipids = value.split(",")
           else:
               lipids.append(value) 
        elif key.strip('@') == "AMOUNTS":
            if "," in key:
                amounts = value.split(",")
            else:
                amounts.append(value)
        else:
            save_info[key.strip('@')] = value

print(lipids)

MOLECULE_AMOUNTS = {}

for i,lipid in enumerate(lipids):
    print(lipids)
    MOLECULE_AMOUNTS[lipid] = amounts[i]

save_info['MOLECULE_AMOUNTS'] = MOLECULE_AMOUNTS
print(save_info)

outfileDICT = str(save_dir) + '/' + 'README.yaml'
    
with open(outfileDICT, 'w') as f:
    yaml.dump(save_info,f, sort_keys=False)



