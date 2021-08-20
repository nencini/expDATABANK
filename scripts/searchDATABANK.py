import os
import yaml
import json
import matplotlib.pyplot as plt
import numpy as np
import math

from matplotlib import cm
from scipy.stats import norm



databank_path = '../../Databank/Data/Simulations'


lipid_numbers_list = ['POPC', 'POPG', 'POPS', 'POPE', 'CHOL', 'DPPC', 'DHMDMAB', 'DMPC', 'POPI'] # should contain all lipid names
ions_list = ['POT', 'SOD', 'CLA', 'CAL'] # should contain names of all ions


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
    
    def getIons(self, ions):
        simIons = []
        
        for key in ions:
           try:
               if self.readme['N'+key] != 0:
                   simIons.append(key)
           except KeyError:
               continue
                
        return simIons

    #fraction of each lipid with respect to total amount of lipids      
    def molarFraction(self, molecule,molecules=lipid_numbers_list): #only for lipids
        sum_lipids = 0
        for key in molecules:
            try:
                #print(readme[key])
                sum_lipids += sum(self.readme['N'+key])
            except KeyError:
                continue
        #print("N" + molecule +": " + str(self.readme['N'+molecule]))
        
        number = sum(self.readme['N'+molecule]) 
    
        return number / sum_lipids
        
        
    # concentration of other molecules than lipids 
    # change name to ionConcentration()
    def ionConcentration(self, molecule, exp_counter_ions):
        lipids1 = self.getLipids()
        c_water = 55.5
        N_water = self.readme['NSOL']
        try:
            N_molecule = self.readme['N'+molecule] #number of ions
        except KeyError:
            N_molecule = 0
        

        lipids2 = []
        if exp_counter_ions and N_molecule != 0:
            for lipid in lipids1:
                if molecule in exp_counter_ions.keys() and lipid == exp_counter_ions[molecule]:
                    N_lipid = self.readme['N'+lipid]
                 #   print(molecule + " " + lipid)
                 #   print(self.readme)
                    lipids2.append(sum(N_lipid))
        
        N_molecule = N_molecule - sum(lipids2)
       # print(N_molecule)
        
        c_molecule = (N_molecule * c_water) / N_water
        #print(c_molecule)
        
        return c_molecule
        
    def totalLipidConcentration(self):
        lipids = self.getLipids()
        c_water = 55.5
        N_water = self.readme['NSOL']
        N_lipids = 0
        for lipid in lipids:
            N_lipids += sum(self.readme['N'+lipid])
        try:
            if (N_water / N_lipids) > 25 :
                tot_lipid_c = 'full hydration'
              #  print('full hydration')
            else:
                tot_lipid_c = (N_lipids * c_water) / N_water
        except ZeroDivisionError:
            print(self.readme)    
        return tot_lipid_c
        
##################
class Experiment:
    def __init__(self, readme, data, dataPath):
        self.readme = readme
        self.data = data #dictionary where key is the lipid type and value is order parameter file
        self.dataPath = dataPath
        
    def getLipids(self, molecules=lipid_numbers_list):
        lipids = []
        for key in molecules:
            try:
                if key in self.readme['MOLAR_FRACTIONS'].keys():
                    lipids.append(key)
            except KeyError:
                continue
        return lipids
    
    def getIons(self, ions):
        expIons = []
        
        for key in ions:
            try:
                if self.readme['ION_CONCENTRATIONS'][key] != 0: 
                    expIons.append(key)
            except KeyError:
                continue    
            try:
                if key in self.readme['COUNTER_IONS'].keys():
                    expIons.append(key)
            except AttributeError:
                continue
        return expIons
 
       
####################plot data #########################
def plotData(simulation, experiment):
    lipids = simulation.getLipids()
    for lipid in lipids:
        fig= plt.figure(figsize=(10,7))

        for key,value in simulation.data[lipid].items():
        #print (key,value[0][0],value[0][2])
            plt.gca().invert_yaxis()
                       
            if lipid == 'POPG' and 'M_G3C6_M' in key:
                plt.plot(0,value[0][0],"s", color='red',marker=".", markersize=10)  #color=colors[i],
                plt.errorbar(0,value[0][0],color='red',yerr=value[0][2])
            if 'M_G3N6' in key:
           #     print(key)
           #     print(value)
           #     print(simulation.indexingPath)
                plt.plot(0,float(value[0][0]),"s",color='red',marker=".", markersize=10)
                plt.errorbar(0,float(value[0][0]),color='red',yerr=float(value[0][2]))
            if 'M_G3C5_M' in key:
                plt.plot(1,value[0][0],"s",color='red',marker=".", markersize=10)
                plt.errorbar(1,value[0][0],color='red',yerr=value[0][2])
            if 'M_G3C4_M' in key:
                plt.plot(2,value[0][0],"s",color='red', marker=".", markersize=10)
                plt.errorbar(2,value[0][0],color='red',yerr=value[0][2])
            if 'M_G3_M' in key:
                plt.plot(3,value[0][0],"s",color='red',marker=".", markersize=10)
                plt.errorbar(3,value[0][0],color='red',yerr=value[0][2])
            if 'M_G2_M' in key:
                plt.plot(4,value[0][0],"s",label=simulation.readme.get('SYSTEM')+" "+simulation.readme.get('FF'),color='red', marker=".", markersize=10)
                plt.errorbar(4,value[0][0],color='red',yerr=value[0][2])
            if 'M_G1_M' in key:
                plt.plot(5,value[0][0],"s",color='red',marker=".", markersize=10)
                plt.errorbar(5,value[0][0],color='red',yerr=value[0][2])

        dataFile = experiment.data
 #   print(dataFile)
        for key,value in simulation.data[lipid].items():
            if lipid == 'POPG' and 'M_G3C6_M' in key:
                plt.plot(0,value[0][0],"s", color='blue',marker=".", markersize=10)  
                plt.errorbar(0,value[0][0], color='blue',yerr=value[0][1])            
            if 'M_G3N6' in key:
                plt.plot(0,value[0][0],"s",color='blue',marker=".", markersize=10)   
                plt.errorbar(0,value[0][0],color='blue',yerr=value[0][1])
            if 'M_G3C5_M' in key:
                plt.plot(1,value[0][0],"s",color='blue',marker=".", markersize=10)
                plt.errorbar(1,value[0][0],color='blue',yerr=value[0][1])
            if 'M_G3C4_M' in key:
                plt.plot(2,value[0][0],"s",color='blue', marker=".", markersize=10)
                plt.errorbar(2,value[0][0],color='blue',yerr=value[0][1])
            if 'M_G3_M' in key:
                plt.plot(3,value[0][0],"s",color='blue',marker=".", markersize=10)
                plt.errorbar(3,value[0][0],color='blue',yerr=value[0][1])
            if 'M_G2_M' in key:
                plt.plot(4,value[0][0],"s",label="experiment " +lipid +" "+str(experiment.readme.get('TEMPERATURE')) + "K",color='blue', marker=".", markersize=10)
                plt.errorbar(4,value[0][0],color='blue',yerr=value[0][1])
            if 'M_G1_M' in key:
                plt.plot(5,value[0][0],"s",color='blue',marker=".", markersize=10)
                plt.errorbar(5,value[0][0],color='blue',yerr=value[0][1])
            
        plt.legend(loc='lower center',ncol=2, fontsize=15, bbox_to_anchor=(0.5, 1.01))
        plt.ylabel('S_CH', fontsize=20)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
#save figures
# make path like in the simulation databank!!
       # plt.savefig('./figures/' + simulation.readme.get('SYSTEM') + '_' + simulation.readme.get('FF') + '.png', bbox_inches='tight')
        save_plot_path = '../Data/QualityEvaluation/OrderParameters/' + simulation.indexingPath + '/'
   #     plt.savefig('./figures/' + simulation.readme.get('SYSTEM') + '.png', bbox_inches='tight')
        plt.savefig(save_plot_path + simulation.readme.get('SYSTEM') + '.png', bbox_inches='tight')
        plt.close()
####################################################################################################
#Quality evaluation
#assume that an order parameter calculated S from a simulation is normally distributed
# OP_sd = sqrt(N)*STEM, where N is number of lipids and STEM is standard error of mean
#def OP_sd(N, STEM):
#    op_sd = math.sqrt(N)*STEM
#    
#    return op_sd

# op_sd = STEM
    
# P: what is the probability that S_exp +/- 0.02 is in g(s) where g(s) is the probability density function of normal distribution N(s, S_sim, S_sim_sd)

#def prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_sd):
#    #normal distribution N(s, OP_sim, op_sim_sd)
#    a = OP_exp - exp_error
#    b = OP_exp + exp_error
    
#    P_S = norm.cdf(b, loc=OP_sim, scale=op_sim_sd) - norm.cdf(a, loc=OP_sim, scale=op_sim_sd)
    
#    return P_S
    
# quality of simulated order parameter
#def OPquality(P_S, op_sim_err):
    
#    quality = P_S / op_sim_err
    
#    return quality

##############################################
#loop over the simulations in the simulation databank and read simulation readme and order parameter files into objects
simulations = []

for subdir, dirs, files in os.walk(r'../../Databank/Data/Simulations/'): 
    for filename1 in files:
        filepath = subdir + os.sep + filename1
        
        if filepath.endswith("README.yaml"):
            READMEfilepathSimulation = subdir + '/README.yaml'
            with open(READMEfilepathSimulation) as yaml_file_sim:
                readmeSim = yaml.load(yaml_file_sim, Loader=yaml.FullLoader)
                indexingPath = "/".join(filepath.split("/")[5:9])
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
                
#loop over the experiment entries in the experiment databank and read experiment readme and order parameter files into objects 
experiments = []                    
for exp_subdir, exp_dirs, exp_files in os.walk(r'../Data/experiments'):
    for filename1 in exp_files:
        filepath_exp = exp_subdir + os.sep + filename1
        if filepath_exp.endswith("README.yaml"):
            READMEfilepathExperiment = exp_subdir + '/README.yaml'
       #     print(READMEfilepathExperiment)
            with open(READMEfilepathExperiment) as yaml_file_exp:
                readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
                expOPdata = {}
                for filename2 in exp_files:
                    if filename2.endswith('Order_Parameters.json'):
                        key_data2 = filename2.replace('_Order_Parameters.json','')
                        #print(key_data2)
                        dataPath = exp_subdir + '/' + filename2 #json!
                        with open(dataPath) as json_file:
                            expOPdata[key_data2] = json.load(json_file)
                            json_file.close()
                experiments.append(Experiment(readmeExp, expOPdata, exp_subdir))
                yaml_file_exp.close()

#Pair each simulation with an experiment with the closest matching temperature and composition
pairs = []

for experiment in experiments: 
 #    print(experiment.readme)
    # check lipid composition matches the simulation
    exp_lipids = experiment.getLipids() 
      #  print(experiment.getLipids())

    exp_total_lipid_concentration = experiment.readme['TOTAL_LIPID_CONCENTRATION']
    exp_ions = experiment.getIons(ions_list)
   # print('experiment ions: ' + str(exp_ions)) 
    exp_counter_ions = experiment.readme['COUNTER_IONS']

    for simulation in simulations:
        sim_lipids = simulation.getLipids()
        # continue if lipids are same
        if set(sim_lipids) == set(exp_lipids):
            sim_ions = simulation.getIons(ions_list)
          #  print("simulation ions " + str(sim_ions)) 
            sim_total_lipid_concentration = simulation.totalLipidConcentration() 
            #print(lipids)
            t_sim = simulation.readme['TEMPERATURE']
            sim_molar_fractions = {}
            #calculate molar fractions from simulation
            for lipid in sim_lipids:
                sim_molar_fractions[lipid] = simulation.molarFraction(lipid)
  #          print(simulation.readme['SYSTEM'])
  #          print(sim_molar_fractions)
  #          print("total lipid concentration: " + str(sim_total_lipid_concentration))

            #calculate concentrations of other molecules
            sim_concentrations = {}
            for molecule in ions_list:
                sim_concentrations[molecule] = simulation.ionConcentration(molecule, exp_counter_ions)
   #         print(sim_concentrations)
    

        
            mf_ok = 0
            for key in sim_lipids:
                if (experiment.readme['MOLAR_FRACTIONS'][key] >= sim_molar_fractions[key] - 0.05) and (experiment.readme['MOLAR_FRACTIONS'][key] <= sim_molar_fractions[key]+ 0.05):
                    mf_ok +=1 

            c_ok = 0


            if set(sim_ions) == set(exp_ions):
              #  print(sim_ions)
              #  print(exp_ions)
                for key in sim_ions:
                    if (experiment.readme['ION_CONCENTRATIONS'][key] >= sim_concentrations[key] - 0.05) and (experiment.readme['ION_CONCENTRATIONS'][key] <= sim_concentrations[key] + 0.05): # onko simulaation ja kokeen ionikonsentraatio tarpeeksi lähellä
                        c_ok += 1 


            switch = 0
            
            if (type(exp_total_lipid_concentration) == float) and (type(sim_total_lipid_concentration) == float): 
                if ((exp_total_lipid_concentration >= sim_total_lipid_concentration - 0.1) and (exp_total_lipid_concentration <= sim_total_lipid_concentration + 0.1)):
                    switch = 1
            elif type(exp_total_lipid_concentration) == str and type(sim_total_lipid_concentration) == str:
                if exp_total_lipid_concentration == sim_total_lipid_concentration:
                    switch = 1
                        
            if switch == 1: 
            #check temperature +/- 2 degrees
                t_exp = experiment.readme['TEMPERATURE']
                if (mf_ok == len(sim_lipids)) and (c_ok == len(sim_ions)) and (t_exp >= float(t_sim) - 2.0) and (t_exp <= float(t_sim) + 2.0):
                    #  print(simulation.indexingPath)
                    pairs.append([simulation, experiment])
                    print(simulation.readme['SYSTEM'])
                    print(simulation.indexingPath)
                    print(experiment.dataPath)
                    #Add path to experiment into simulation README.yaml
                    simulation.readme['EXPERIMENT'] = "/".join(experiment.dataPath.split("/")[3:6])         #"/".join(filepath.split("/")[6:10])
                    print(simulation.readme['EXPERIMENT'])
                    outfileDICT = '../../Databank/Data/Simulations/'+ simulation.indexingPath + '/README.yaml'
    
                    with open(outfileDICT, 'w') as f:
                        yaml.dump(simulation.readme,f, sort_keys=False)
print("Found " + str(len(pairs)) + " pairs")  
for pair in pairs:
    print(pair[0].readme)
    print(pair[1].readme)
                       

########################ORDER PARAMETER QUALITY ANALYSIS#######################################                
                
#make file paths for saving quality evaluation 
#os.system('mkdir ../Data/QualityEvaluation')
#os.system('mkdir ../Data/QualityEvaluation/OrderParameters')
#os.system('mkdir ../Data/QualityEvaluation/')

#for pair in pairs:
#    sub_dirs = pair[0].indexingPath.split("/")
#    os.system('mkdir ../Data/QualityEvaluation/' + sub_dirs[0])
#    os.system('mkdir ../Data/QualityEvaluation/' + sub_dirs[0] + '/' + sub_dirs[1])
#    os.system('mkdir ../Data/QualityEvaluation/' + sub_dirs[0] + '/' + sub_dirs[1] + '/' + sub_dirs[2])
    
#    os.system('mkdir ../Data/QualityEvaluation/' + sub_dirs[0] + '/' + sub_dirs[1] + '/' + sub_dirs[2] + '/' + sub_dirs[3])
#    DATAdir = '../Data/QualityEvaluation/' + str(sub_dirs[0]) + '/' + str(sub_dirs[1]) + '/' + str(sub_dirs[2]) + '/' + str(sub_dirs[3])
    
    # measured values do not exist for all CH!!!!! need to take mapping name and find matching CH from simulation
    # read existing experimental values and mapping names to match with simulated CH bonds
#    print(pair[0].indexingPath)
#    for lipid1 in pair[0].getLipids():
#        print(lipid1)
#        print(pair[1].dataPath)
#        print(pair[1].data.keys())
#        OP_data_lipid = pair[0].data[lipid1]
        #print(OP_data_lipid)
#        N = float(sum(pair[0].readme['N'+lipid1])) # number of lipids per lipid type
#        print(N)
#        try:
#            exp_OP_data = pair[1].data[lipid1]
#        except KeyError:
#            print("Experimental order parameter data do not exist for lipid " + lipid1 + ".")
#            continue
#        exp_error = 0.02

#        OP_qual_data = {}
        
#        for key, value in exp_OP_data.items():
#            if value[0] is not 'NaN':
#                OP_exp = value[0][0]
#                OP_sim = OP_data_lipid[key][0][0]
#           #     print(OP_data_lipid[key])
#                op_sim_sd = OP_data_lipid[key][0][2] #standard error of mean
#                print(key)
#                print(op_sim_STEM) 
#                op_sim_err = OP_data_lipid[key][0][1]
#                op_sim_sd =  OP_sd(N, op_sim_STEM)

#                S_prob = prob_S_in_g(OP_exp, exp_error, OP_sim, op_sim_sd)
            
#                OP_qual_data[key] = OPquality(S_prob, op_sim_err)
                
#        outfile = DATAdir + '/' + lipid1 + '_OP_quality.json'
        
#        with open(outfile, 'w') as f:
#            json.dump(OP_qual_data,f)

#        f.close()
                
        
                
                
            
            
            
            
            










