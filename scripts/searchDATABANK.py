import os
import yaml
import json
import matplotlib.pyplot as plt

from matplotlib import cm


databank_path = '../../NMRLipids_Databank/Databank/Data/Simulations'


lipid_numbers_list = ['POPC', 'POPG', 'POPS', 'POPE', 'CHOL', 'DPPC', 'DHMDMAB', 'DMPC', 'POPI'] # should contain all lipid names
other_molecules = ['POT', 'SOD', 'CLA', 'CAL'] # should contain names of all the other molecules than lipids


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
    def moleculeConcentration(self, molecule):
        c_water = 55.5
        N_water = self.readme['NSOL']
        try:
            N_molecule = self.readme['N'+molecule]
        except KeyError:
            N_molecule = 0
            
        c_molecule = (N_molecule * c_water) / N_water
        
        return c_molecule
        
    def totalLipidConcentration(self):
        lipids = self.getLipids()
        c_water = 55.5
        N_water = self.readme['NSOL']
        N_lipids = 0
        for lipid in lipids:
            N_lipids += sum(self.readme['N'+lipid])
        try:
            if (N_water / N_lipids) > 30 :
                tot_lipid_c = 'full hydration'
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
##################



#???????
def getMolecules(readme,molecules=lipid_numbers_list):
    molecules = []
    for key in molecules:
        try:
            if readme['N'+key] != [0,0] or readme['N'+key] != 0: 
                molecules.append(key)
        except KeyError:
            continue

    return molecules
    
 
       
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

##############################################
#loop over the simulations in the simulation databank
simulations = []

for subdir, dirs, files in os.walk(r'../../NMRLipids_Databank/Databank/Data/Simulations/'): 
    for filename1 in files:
        filepath = subdir + os.sep + filename1
        
        if filepath.endswith("README.yaml"):
            READMEfilepathSimulation = subdir + '/README.yaml'
            with open(READMEfilepathSimulation) as yaml_file_sim:
                readmeSim = yaml.load(yaml_file_sim, Loader=yaml.FullLoader)
                indexingPath = "/".join(filepath.split("/")[6:10])
                simOPdata = {} #order parameter files for each lipid
                for filename2 in files:
                    #print(filename2)
                    if filename2.endswith('OrderParameters.json'):
                        key_data1 = filename2.replace('OrderParameters.json', '')
                        OPfilepath = subdir + '/' + filename2
                        with open(OPfilepath) as json_file:
                            simOPdata[key_data1] = json.load(json_file) #modify to contain OP data of each lipid in the membrane thst it works for mixtures
                            json_file.close()
                simulations.append(Simulation(readmeSim, simOPdata, indexingPath)) #toimiiks?
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
                expOPdata = {}
                for filename2 in exp_files:
                    if filename2.endswith('Order_Parameters.json'):
                        key_data2 = filename2.replace('Order_Parameters.json','')
                        dataPath = exp_subdir + '/Order_Parameters.json' #json!
                        with open(dataPath) as json_file:
                            expOPdata[key_data2] = json.load(json_file)
                            json_file.close()
                experiments.append(Experiment(readmeExp, expOPdata, exp_subdir))
                yaml_file_exp.close()

#Pair each simulation with an experiment with the closest matching temperature and composition
pairs = []

for simulation in simulations:
    sim_lipids = simulation.getLipids()
    sim_molecules = getMolecules(simulation.readme, molecules= other_molecules)
    sim_total_lipid_concentration = simulation.totalLipidConcentration() 
    #print(lipids)
    t_sim = simulation.readme['TEMPERATURE']
    sim_molar_fractions = {}
    #calculate molar fractions from simulation
    for lipid in sim_lipids:
        sim_molar_fractions[lipid] = simulation.molarFraction(lipid)
    print(simulation.readme['SYSTEM'])
    print(sim_molar_fractions)
    print("total lipid concentration: " + str(sim_total_lipid_concentration))

    #calculate concentrations of other molecules
    sim_concentrations = {}
    for molecule in other_molecules:
        if (molecule != 'SOL'):
            sim_concentrations[molecule] = simulation.moleculeConcentration(molecule)
    print(sim_concentrations)
    
    for experiment in experiments: 
    #    print(experiment.readme)
    # check lipid composition matches the simulation
        exp_lipids = experiment.getLipids() ### MUUT MOLEKYYLIT?
      #  print(experiment.getLipids())

        exp_total_lipid_concentration = experiment.readme['TOTAL_LIPID_CONCENTRATION']

        if set(sim_lipids) == set(exp_lipids):
            mf_ok = 0
            for key in sim_lipids:
                if (experiment.readme['MOLAR_FRACTIONS'][key] >= sim_molar_fractions[key] - 0.05) and (experiment.readme['MOLAR_FRACTIONS'][key] <= sim_molar_fractions[key]+ 0.05):
                    mf_ok +=1 

            c_ok = 0

            for key in sim_molecules: 
                if (experiment.readme['ION_CONCENTRATIONS'][key] >= sim_concentrations[key] - 0.05) and (experiment.readme['ION_CONCENTRATIONS'][key] <= sim_concentrations[key] + 0.05):
                    c_ok += 1
               #     print("TOIMII") #TOIMII

            switch = 0
      #      print("koe:")
      #      print(type(exp_total_lipid_concentration))
      #      print("simulaatio:")
      #      print(type(sim_total_lipid_concentration))
            
            if (type(exp_total_lipid_concentration) == float) and (type(sim_total_lipid_concentration) == float): 
                if ((exp_total_lipid_concentration >= sim_total_lipid_concentration - 0.1) and (exp_total_lipid_concentration <= sim_total_lipid_concentration + 0.1)):
                    switch = 1
            elif type(exp_total_lipid_concentration) == str and type(sim_total_lipid_concentration) == str:
                if exp_total_lipid_concentration == sim_total_lipid_concentration:
                    switch = 1
                        
            if switch == 1: 
         #       print("toimii")
            #check temperature +/- 2 degrees
                t_exp = experiment.readme['TEMPERATURE']
                if (mf_ok == len(sim_lipids)) and (c_ok == len(sim_molecules)) and (t_exp >= float(t_sim) - 2.0) and (t_exp <= float(t_sim) + 2.0):
                    print("Onnistui jee!!")
                    #  print(simulation.indexingPath)
                    pairs.append([simulation, experiment])
                    print(simulation.readme['SYSTEM'])
                    print(experiment.dataPath)
                    #Add path to experiment into simulation README.yaml
                    simulation.readme['EXPERIMENT'] = experiment.dataPath +'/'
                
                    outfileDICT = '../../NMRLipids_Databank/Databank/Data/Simulations/'+ simulation.indexingPath + '/README.yaml'
    
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
#    plotData(pair[0],pair[1])
            
    print(pair[0].readme)
    print(pair[1].readme)
print(len(pairs))








