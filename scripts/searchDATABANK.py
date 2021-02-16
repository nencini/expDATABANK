#Search and compare simulation order parameter data with experimental data
import os
import yaml
import json
import matplotlib.pyplot as plt
from matplotlib import cm

#simulation data directory
#simPath = "/media/akiirikk/DATADRIVE1/tietokanta/NMRlipidsIVPEandPG/Data/Simulations/"

#search parameters
#T_RANGE is range of temperatures to be searched
#mixture of lipids????
searchParameters = {'LIPID' : 'POPC', 'T_RANGE': [290,320] ,'AMOUNT' : 1}

#dictionaries
lipid_numbers_list = [ 'NPOPC', 'NPOPG', 'NPOPS', 'NPOPE', 'NCHOL'] 


class Simulation:
    def __init__(self, readme, data):
        self.readme = readme
        self.data = data
    


class Experiment:
    def __init__(self, readme, data):
        self.readme = readme
        self.data = data




#returns the relative amount of molecule in the entire membrane not per single leaflet
def relativeAmount(readme, molecule,molecules=lipid_numbers_list):
    sum_molecules = 0
    for key in molecules:
        try:
            #print(readme[key])
            sum_molecules = sum_molecules + sum(readme[key])
        except KeyError:
            continue
    
    number = sum(readme.get('N'+molecule))
    return number / sum_molecules

#search simulations 
#returns list of simulations that match to search parameters
def search_simulations(searchParameters):
    simulations = []
    for subdir, dirs, files in os.walk(r'../../NMRlipidsIVPEandPG/Data/Simulations/'): 
        for filename1 in files:
            filepath = subdir + os.sep + filename1
            if filepath.endswith("README.yaml"):
                READMEfilepath = subdir + '/README.yaml'
                with open(READMEfilepath) as yaml_file:
                    readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
                    lipid = searchParameters['LIPID']
                    lipid_amount = relativeAmount(readme, lipid)
                    #print(readme.get('SYSTEM'))
                    #print(filepath)
                    #print(lipid_amount)
                    tsim = int(readme.get('TEMPERATURE'))
                    tmin = searchParameters['T_RANGE'][0]
                    tmax = searchParameters['T_RANGE'][1]
                    #check that lipid, its amount and temperature match the search parameters 
                    if (readme.get('N'+lipid) != 0) and (lipid_amount == searchParameters['AMOUNT']) and (tsim >= tmin and tsim <= tmax):
                        print(readme.get('SYSTEM'))
                        print(filepath)
                        print(lipid_amount)
                        print(readme)
                    #
                        for filename2 in files:
                            #print(filename2)
                            if filename2.endswith('OrderParameters.json'):
                                OPfilepath = subdir + '/' + filename2
                                #print(OPfilepath)
                                with open(OPfilepath) as json_file:
                                    data = json.load(json_file)
                                simulations.append(Simulation(readme, data))

    return simulations
                
###########

#search experimental data
def search_exp(searchParameters):
    experiments = []
    for subdir, dirs, files in os.walk(r'../Data/'):
       for filename in files:
           filepath = subdir + os.sep + filename
           if filepath.endswith("README.yaml"):
               READMEfilepath = subdir + '/README.yaml'
               with open(READMEfilepath) as yaml_file:
                   readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
                   lipid = searchParameters['LIPID']
                   #print(readme)
                   molecules = readme['MOLECULE_AMOUNTS']
                  # print(molecules.keys())
                   texp = int(readme.get('TEMPERATURE'))
                   tmin = searchParameters['T_RANGE'][0]
                   tmax = searchParameters['T_RANGE'][1]
                   try:
                       if float(molecules[lipid])==searchParameters['AMOUNT'] and (texp >= tmin and texp <= tmax):
                           print(readme)
                           dataPath = subdir + '/Order_Parameters.dat'
                           experiments.append(Experiment(readme,dataPath)) #saves the path to data
                   except KeyError:
                       continue
    return experiments

sims = search_simulations(searchParameters)
print("simulations: " + str(len(sims)))


exp = search_exp(searchParameters)

print("experiments: " + str(len(exp))) 

#define n colors for plotting simulation 
#n = len(sims) + len(exp)
#colormap = plt.cm.get_cmap('hsv')
#colors = [colormap(i) for i in range(0,n)]

#fig= plt.figure(figsize=(10,7))
#plots each simulation against experiment in its own graph
#plotting works for one experiment if there are more experimental data matching with the simulations i don't know what to do
experiment = exp[0]
for i,simulation in enumerate(sims):
    fig= plt.figure(figsize=(10,7))
    lipid = searchParameters['LIPID']
    for key,value in simulation.data.items():
    #print (key,value[0][0],value[0][2])
        plt.gca().invert_yaxis()
                        
        if lipid == 'POPG' and 'M_G3C6_M' in key:
            plt.plot(0,value[0][0],"s", color='red',marker=".", markersize=10)  #color=colors[i],
            plt.errorbar(0,value[0][0],color='red',yerr=value[0][2])
        if 'M_G3N6' in key:
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
    

        
    # plot experimental data   
    #lipid = searchParameters['LIPID']
    dataFile = experiment.data
 #   print(dataFile)
    with open(dataFile, 'r') as f:
        for line in f:
            if ('#' not in line) and (line is not '\n') and ('nan' not in line):
                #print(line)
                line = line.replace('\n','')
                values = line.split(" ")[-2:] #takes value and error value from the end of line. works for the head group but error values of the tails are missing
                
                values = [float(i) for i in values]
                if lipid == 'POPG' and 'M_G3C6_M' in line:
                    plt.plot(0,values[0],"s", color='blue',marker=".", markersize=10)  
                    plt.errorbar(0,values[0], color='blue',yerr=values[1])            
                if 'M_G3N6' in line:
                    plt.plot(0,values[0],"s",color='blue',marker=".", markersize=10)   
                    plt.errorbar(0,values[0],color='blue',yerr=values[1])
                if 'M_G3C5_M' in line:
                    plt.plot(1,values[0],"s",color='blue',marker=".", markersize=10)
                    plt.errorbar(1,values[0],color='blue',yerr=values[1])
                if 'M_G3C4_M' in line:
                    plt.plot(2,values[0],"s",color='blue', marker=".", markersize=10)
                    plt.errorbar(2,values[0],color='blue',yerr=values[1])
                if 'M_G3_M' in line:
                    plt.plot(3,values[0],"s",color='blue',marker=".", markersize=10)
                    plt.errorbar(3,values[0],color='blue',yerr=values[1])
                if 'M_G2_M' in line:
                    plt.plot(4,values[0],"s",label="experiment " +lipid +" "+experiment.readme.get('TEMPERATURE') + "K",color='blue', marker=".", markersize=10)
                    plt.errorbar(4,values[0],color='blue',yerr=values[1])
                if 'M_G1_M' in line:
                    plt.plot(5,values[0],"s",color='blue',marker=".", markersize=10)
                    plt.errorbar(5,values[0],color='blue',yerr=values[1])
            
    plt.legend(loc='lower center',ncol=2, fontsize=15, bbox_to_anchor=(0.5, 1.01))
    plt.ylabel('S_CH', fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
#save figures
    plt.savefig('./figures/' + simulation.readme.get('SYSTEM') + '_' + simulation.readme.get('FF') + '.png', bbox_inches='tight')
    plt.close()







