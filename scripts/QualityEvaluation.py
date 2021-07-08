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
def OPquality(P_S): #, op_sim_STEM):
    
    quality = P_S #/ (op_sim_STEM*op_sim_STEM)
    quality_float = quality.item()
    
    return quality_float
    
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
#    DATAdir = '../../NMRLipids_Data/Databank/Data/Simulations/' + str(sub_dirs[0]) + '/' + str(sub_dirs[1]) + '/' + str(sub_dirs[2]) + '/' + str(sub_dirs[3])
    DATAdir = '../Data/QualityEvaluation/' + str(sub_dirs[0]) + '/' + str(sub_dirs[1]) + '/' + str(sub_dirs[2]) + '/' + str(sub_dirs[3])
    # measured values do not exist for all CH!!!!! need to take mapping name and find matching CH from simulation
    # read existing experimental values and mapping names to match with simulated CH bonds

    for lipid1 in simulation.getLipids():
        print(lipid1)
        print(simulation.indexingPath)
        print(simulation.data.keys())
        
        OP_data_lipid = simulation.data[lipid1]
        
        #print(OP_data_lipid)

        # get readme file of the experiment
        experimentFilepath = simulation.readme['EXPERIMENT']
        READMEfilepathExperiment  = experimentFilepath + 'README.yaml'
        experiment = Experiment()
        with open(READMEfilepathExperiment) as yaml_file_exp:
            readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
            experiment.readme = readmeExp
            #print(experiment.readme)
        yaml_file_exp.close()

        lipidExpOPdata = {}
        try:
            exp_OP_filepath = simulation.readme['EXPERIMENT'] + lipid1 + '_Order_Parameters.json'
        #print(exp_OP_filepath)
            with open(exp_OP_filepath) as json_file:
                lipidExpOPdata = json.load(json_file)
            json_file.close()
            
            simulationREADMEsave = DATAdir + '/README.yaml'
            with open(simulationREADMEsave, 'w') as f:
                yaml.dump(simulation.readme,f, sort_keys=False)
            f.close()
        except FileNotFoundError:
            print("Experimental order parameter data do not exist for lipid " + lipid1 + ".")
            continue

        exp_error = 0.02

        OP_qual_data = {}
        
        for key, value in lipidExpOPdata.items():
            if lipidExpOPdata[key][0][0] is not 'NaN':
                OP_array = [float(x) for x in OP_data_lipid[key][0]] #convert elements to float because in some files the elements are strings
                #print(OP_array)
                #print(type(OP_array))
                OP_exp = value[0][0]
                OP_sim = OP_array[0]
                op_sim_sd = OP_array[2] #standard error of mean
                op_sim_STEM = OP_array[2] 

                S_prob = prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_sd)
             
                op_quality = OPquality(S_prob) #, op_sim_STEM) #numpy float must be converted to float
                #print(type(op_quality))
                OP_array.append(op_quality)
                #print(OP_array)
                
                OP_qual_data[key] = OP_array
                
#            else:
#                OP_qual_data[key] = value.append('nan')
                
#        outfile = DATAdir + '/' + lipid1 + '_OP_quality.json'
        
#        with open(outfile, 'w') as f:
#            json.dump(OP_qual_data,f)

#        f.close()  
        
        
        # quality data should be added into the OrderParameters.json file of the simulation                  
        outfile = DATAdir + '/' + lipid1 + '_OrderParameters.json'
        
        with open(outfile, 'w') as f:
            json.dump(OP_qual_data,f)
        f.close()
        
       
        
        
        
      #  print(OP_qual_data)                        
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
                
