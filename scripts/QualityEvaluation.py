import os
import yaml
import json
import matplotlib.pyplot as plt
import numpy as np
import math

from matplotlib import cm
from scipy.stats import norm

lipid_numbers_list = ['POPC', 'POPG', 'POPS', 'POPE', 'CHOL', 'DPPC', 'DHMDMAB', 'DMPC', 'POPI'] # should contain all lipid names
#################################
class Simulation:
    def __init__(self, readme, data, indexingPath):
        self.readme = readme
        self.data = data #dictionary where key is the lipid type and value is order parameter file
        self.indexingPath = indexingPath
        
    def getLipids(self, molecules=lipid_numbers_list):
        lipids = []
        for key in molecules:
            try:
                if self.readme['N'+key] != [0,0]: 
                    lipids.append(key)
            except KeyError:
                continue
        return lipids     
        
class Experiment:
    pass   

#################################
#Quality evaluation of simulated data
#assume that an order parameter calculated S from a simulation is normally distributed
# OP_sd = sqrt(N)*STEM, where N is number of lipids and STEM is standard error of mean
#def OP_sd(N, STEM):
#    op_sd = math.sqrt(N)*STEM
#    
#    return op_sd

# op_sd = STEM
    
# P: what is the probability that S_exp +/- 0.02 is in g(s) where g(s) is the probability density function of normal distribution N(s, S_sim, S_sim_sd)

def prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_sd):
    #normal distribution N(s, OP_sim, op_sim_sd)
    a = OP_exp - exp_error
    b = OP_exp + exp_error
    
    P_S = norm.cdf(b, loc=OP_sim, scale=op_sim_sd) - norm.cdf(a, loc=OP_sim, scale=op_sim_sd)
    
    return P_S
    
# quality of simulated order parameter
def OPquality(P_S, op_sim_err):
    
    quality = P_S / op_sim_err
    
    return quality
    
###################################################################################################
simulations = []
for subdir, dirs, files in os.walk(r'../../NMRLipids_Databank/Databank/Data/Simulations/'): 
    for filename1 in files:
        filepath = subdir + os.sep + filename1
        
        if filepath.endswith("README.yaml"):
            READMEfilepathSimulation = subdir + '/README.yaml'
            with open(READMEfilepathSimulation) as yaml_file_sim:
                readmeSim = yaml.load(yaml_file_sim, Loader=yaml.FullLoader)
                indexingPath = "/".join(filepath.split("/")[6:10])
                try:
                    if readmeSim['EXPERIMENT']:
                        simOPdata = {} #order parameter files for each type of lipid
                        for filename2 in files:
                            if filename2.endswith('OrderParameters.json'):
                                key_data1 = filename2.replace('OrderParameters.json', '')
                                OPfilepath = subdir + '/' + filename2
                                with open(OPfilepath) as json_file:
                                    simOPdata[key_data1] = json.load(json_file)
                                    json_file.close()
                        simulations.append(Simulation(readmeSim, simOPdata, indexingPath))
                        yaml_file_sim.close()
                except KeyError:
                    print("No matching experimental data for system " + readmeSim['SYSTEM'] + " in directory " + indexingPath)
                    continue
                    
os.system('mkdir ../Data/QualityEvaluation')
os.system('mkdir ../Data/QualityEvaluation/')

for simulation in simulations:
    sub_dirs = simulation.indexingPath.split("/")
    os.system('mkdir ../Data/QualityEvaluation/' + sub_dirs[0])
    os.system('mkdir ../Data/QualityEvaluation/' + sub_dirs[0] + '/' + sub_dirs[1])
    os.system('mkdir ../Data/QualityEvaluation/' + sub_dirs[0] + '/' + sub_dirs[1] + '/' + sub_dirs[2])
    
    os.system('mkdir ../Data/QualityEvaluation/' + sub_dirs[0] + '/' + sub_dirs[1] + '/' + sub_dirs[2] + '/' + sub_dirs[3])
    DATAdir = '../Data/QualityEvaluation/' + str(sub_dirs[0]) + '/' + str(sub_dirs[1]) + '/' + str(sub_dirs[2]) + '/' + str(sub_dirs[3])
    
    # measured values do not exist for all CH!!!!! need to take mapping name and find matching CH from simulation
    # read existing experimental values and mapping names to match with simulated CH bonds

    for lipid1 in simulation.getLipids():
        print(lipid1)
        print(simulation.indexingPath)
        print(simulation.data.keys())
        
        #print(OP_data_lipid)
        
        experimentFilepath = simulation.readme['EXPERIMENT']
        READMEfilepathExperiment  = experimentFilepath + 'README.yaml'
        experiment = Experiment()
        with open(READMEfilepathExperiment) as yaml_file_exp:
            readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
            experiment.readme = readmeExp
            expOPdata = {}
            for filename2 in experimentFilepath:
                if filename2.endswith('Order_Parameters.json'):
                   key_data2 = filename2.replace('_Order_Parameters.json','')
                   #print(key_data2)
                   dataPath = experimentFilepath + '/' + filename2 
                   with open(dataPath) as json_file:
                       expOPdata[key_data2] = json.load(json_file)
                       json_file.close()
            experiment.data = expOPdata
#            experiment = Experiment(readmeExp, expOPdata, exp_subdir)
            yaml_file_exp.close()
        
        try:
            exp_OP_data = experiment.data[lipid1]
        except KeyError:
            print("Experimental order parameter data do not exist for lipid " + lipid1 + ".")
            continue
        exp_error = 0.02

        OP_qual_data = {}
        
        for key, value in exp_OP_data.items():
            if value[0] is not 'NaN':
                OP_array = OP_data_lipid[key].astype(np.float) # convert elements to float because these can be string in some files
                OP_exp = value[0][0]
                OP_sim = OP_array[0][0]
           #     print(OP_data_lipid[key])
                op_sim_sd = OP_array[0][2] #standard error of mean
                print(key)
                print(op_sim_STEM) 
                op_sim_err = OP_array[0][1]
                op_sim_sd =  OP_sd(N, op_sim_STEM)

                S_prob = prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_sd)
            
                OP_qual_data[key] = OPquality(S_prob, op_sim_err)
                
        outfile = DATAdir + '/' + lipid1 + '_OP_quality.json'
        
        with open(outfile, 'w') as f:
            json.dump(OP_qual_data,f)

        f.close()                    
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
