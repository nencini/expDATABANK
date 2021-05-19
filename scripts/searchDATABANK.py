import os
import yaml
import json
import matplotlib.pyplot as plt

from matplotlib import cm


databank_path = '../../NMRLipids_Databank_updated/Databank/Data/Simulations'


lipid_numbers_list = [ 'POPC', 'POPG', 'POPS', 'POPE', 'CHOL'] # should contain all lipid names
other_molecules = ['POT', 'SOD', 'CLA', 'CAL', 'SOL'] # should contain names of all the other molecules than lipids

class Simulation:
    def __init__(self, readme, data, indexingPath):
        self.readme = readme
        self.data = data
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
##################
class Experiment:
    def __init__(self, readme, data, dataPath):
        self.readme = readme
        self.data = data
        self.dataPath = dataPath
        
    def getLipids(self, molecules=lipid_numbers_list):
        lipids = []
        for key in molecules:
            try:
                if key in self.readme['MOLECULE_FRACTIONS'].keys():
                    lipids.append(key)
            except KeyError:
                continue
        return lipids
##################

def molecularFraction(readme, molecule,molecules=lipid_numbers_list):
    sum_molecules = 0
    for key in molecules:
        try:
            #print(readme[key])
            sum_molecules = sum_molecules + sum(readme['N'+key])
        except KeyError:
            continue
    
    number = sum(readme.get('N'+molecule)) 
    return number / sum_molecules

def getLipids(readme,molecules=lipid_numbers_list):
    lipids = []
    for key in molecules:
        try:
            if readme['N'+key] != [0,0]: 
                lipids.append(key)
        except KeyError:
            continue

    return lipids
       
####################plot data #########################
def plotData(simulation, experiment):
    lipids = getLipids(simulation.readme)
    for lipid in lipids:
        fig= plt.figure(figsize=(10,7))

        for key,value in simulation.data.items():
        #print (key,value[0][0],value[0][2])
            plt.gca().invert_yaxis()
                        
            if lipid == 'POPG' and 'M_G3C6_M' in key:
                plt.plot(0,value[0][0],"s", color='red',marker=".", markersize=10)  #color=colors[i],
                plt.errorbar(0,value[0][0],color='red',yerr=value[0][2])
            if 'M_G3N6' in key:
           #     print(key)
           #     print(value)
                plt.plot(0,value[0][0],"s",color='red',marker=".", markersize=10)
                plt.errorbar(0,value[0][0],color='red',yerr=value[0][2])
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
        for key,value in simulation.data.items():
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

##############################################
#loop over the simulations in the simulation databank
simulations = []

for subdir, dirs, files in os.walk(r'../../NMRLipids_Databank_updated/Databank/Data/Simulations/'): 
    for filename1 in files:
        filepath = subdir + os.sep + filename1
        
        if filepath.endswith("README.yaml"):
            READMEfilepathSimulation = subdir + '/README.yaml'
            with open(READMEfilepathSimulation) as yaml_file_sim:
                readmeSim = yaml.load(yaml_file_sim, Loader=yaml.FullLoader)
                indexingPath = "/".join(filepath.split("/")[6:10])
                for filename2 in files:
                    #print(filename2)
                    if filename2.endswith('OrderParameters.json'):
                        OPfilepath = subdir + '/' + filename2
                        with open(OPfilepath) as json_file:
                            simOPdata = json.load(json_file) #modify to contain OP data of each lipid in the membrane thst it works for mixtures
                            simulations.append(Simulation(readmeSim, simOPdata, indexingPath))
                            json_file.close()
                yaml_file_sim.close()
                    
experiments = []                    
for exp_subdir, exp_dirs, exp_files in os.walk(r'../Data/experiments'):
    for filename1 in exp_files:
        filepath_exp = exp_subdir + os.sep + filename1
        if filepath_exp.endswith("README.yaml"):
            READMEfilepathExperiment = exp_subdir + '/README.yaml'
       #     print(READMEfilepathExperiment)
       #     with open(READMEfilepathExperiment) as yaml_file_exp:
       #         readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
       #         dataPath = exp_subdir + '/Order_Parameters.dat' #json!
       #         experiments.append(Experiment(readmeExp, dataPath))
       #         yaml_file_exp.close()
            with open(READMEfilepathExperiment) as yaml_file_exp:
                readmeExp = yaml.load(yaml_file_exp, Loader=yaml.FullLoader)
                for filename2 in exp_files:
                    if filename2.endswith('Order_Parameters.json'):
                        dataPath = exp_subdir + '/Order_Parameters.json' #json!
                        with open(dataPath) as json_file:
                            expOPdata = json.load(json_file)
                            experiments.append(Experiment(readmeExp, expOPdata, exp_subdir)) #???
                            json_file.close()
                yaml_file_exp.close()

#Pair simulations with a matching experiment
pairs = []

for simulation in simulations:
    sim_lipids = simulation.getLipids()
    #print(lipids)
    t_sim = simulation.readme['TEMPERATURE']
    sim_mol_fractions = {}
    for lipid in sim_lipids:
        sim_mol_fractions[lipid] = molecularFraction(simulation.readme, lipid)
    #print(sim_mol_fractions)
    
    for experiment in experiments: 
    #    print(experiment.readme)
    # check lipids match
        exp_lipids = experiment.getLipids()
        if set(sim_lipids) == set(exp_lipids):
    # check molecular fractions +/-5% lipids, ions, water checked in the same loop 
            for molecule in sim_mol_fractions.keys():
             #   print(sim_mol_fractions[molecule])
             #   print(experiment.readme['MOLECULE_FRACTIONS'][molecule])
                mf_ok = 0
                if (experiment.readme['MOLECULE_FRACTIONS'][molecule] >= sim_mol_fractions[molecule] - 0.05) and (experiment.readme['MOLECULE_FRACTIONS'][molecule] <= sim_mol_fractions[molecule] + 0.05):
                    mf_ok += 1
            #check temperature +/- 2 degrees
            t_exp = experiment.readme['TEMPERATURE']
        #    print(t_exp)
        #    print(simulation.indexingPath)
        #    print(t_sim)
            if (mf_ok == len(sim_mol_fractions.keys())) and (t_exp >= float(t_sim) - 2.0) and (t_exp <= float(t_sim) + 2.0):
              #  print(simulation.indexingPath)
                pairs.append([simulation, experiment])
                print(simulation.readme['SYSTEM'])
                print(experiment.dataPath)
                #Add path to experiment into simulation README.yaml
                simulation.readme['EXPERIMENT'] = experiment.dataPath +'/'
                
                outfileDICT = '../../NMRLipids_Databank_updated/Databank/Data/Simulations/'+ simulation.indexingPath + '/README.yaml'
    
                with open(outfileDICT, 'w') as f:
                    yaml.dump(simulation.readme,f, sort_keys=False)
                
#make file paths for saving quality evaluation plots
os.system('mkdir ../Data/QualityEvaluation')     
os.system('mkdir ../Data/QualityEvaluation/OrderParameters')


for pair in pairs:
    sub_dirs = pair[0].indexingPath.split("/")
    os.system('mkdir ../Data/QualityEvaluation/OrderParameters/' + sub_dirs[0])
    os.system('mkdir ../Data/QualityEvaluation/OrderParameters/' + sub_dirs[0] + '/' + sub_dirs[1])
    os.system('mkdir ../Data/QualityEvaluation/OrderParameters/' + sub_dirs[0] + '/' + sub_dirs[1] + '/' + sub_dirs[2])
    os.system('mkdir ../Data/QualityEvaluation/OrderParameters/' + sub_dirs[0] + '/' + sub_dirs[1] + '/' + sub_dirs[2] + '/' + sub_dirs[3])
    
    #plot order parameters
    plotData(pair[0],pair[1])
            

print(len(pairs))








