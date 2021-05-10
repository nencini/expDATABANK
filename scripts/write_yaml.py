#Organize experimental data into a databank
#Amount of each molecule in the membrane is given as a ratio. For a membrane of single molecule it's 1.

import os
import numpy as np

import yaml
import json

save_dir = "/media/akiirikk/DATADRIVE1/tietokanta/expDATABANK/Data/POPS/T298K"

#experiment information

save_info = {}

save_info['TEMPERATURE'] = 298
save_info['MOLECULE_FRACTIONS'] = {'POPS':1}
print(save_info)

outfileDICT = str(save_dir) + '/' + 'README.yaml'
    
with open(outfileDICT, 'w') as f:
    yaml.dump(save_info,f, sort_keys=False)


