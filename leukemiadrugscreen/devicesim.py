import numpy as np
import matplotlib as plt
import configparser
import pandas as pd
import leukemiadrugscreen as leukdev
import os
import shutil
from leukemiadrugscreen.doseresponsesurf import greco3D as greco
from statsmodels.stats.multicomp import pairwise_tukeyhsd


'''
The function of this script is to generate simulated data to use with
devicezonetest.py analysis script. The structure of the resulting data
is a set of excel files contained within a folder.
Requires the concentration map in .csv format.
'''

# Load in map file
def load_map():
    """Loads in the concentration map and randomly samples
    n=num_cells from it.
    Columns x, y, c1, c2, c3
    Returns pandas dataframe"""
    concentration_map_path = os.path.join(os.path.dirname(leukdev.__file__), "concentration_map.csv")
    concMap = pd.read_csv(concentration_map_path)
    return concMap

def sample_cells(concMap, num_cells):
    cells = concMap.sample(num_cells)
    return cells

def import_config():
    config = configparser.ConfigParser()
    config.read('configurations.conf')
    return config

def export_config():
    ftpUrl = "demoftp.codeteddy.com"
    userName = "codeteddy"
    password = "my#supersecret#password"

def generate_config():
    config_file = configparser.ConfigParser()

    config_file.add_section("FTPSettings")
    config_file.set("FTPSettings", "ftpUrl", "demoftp.codeteddy.com")
    config_file.set("FTPSettings", "userName", "codeteddy")
    config_file.set("FTPSettings", "password", "my#supersecret#password")

    config_file["Drug1"] = {
        "Name":"Daunorubicin",
        "Units":"nM",
        "EC50":50,
        "Emax":0.2,
        "Hill":2
    }
    
    config_file["Drug2"] = {
        "Name":"Vincristine",
        "Units":"nM",
        "EC50":100,
        "Emax":0.3,
        "Hill":3
    }


    with open(r"configurations.conf", 'w') as configfileObj:
        config_file.write(configfileObj)
        configfileObj.flush()
        configfileObj.close()

    print("Test to see configuration file written")

def musyc(conc1, conc2, conc3, h1, h2, h3, E0, E1, E2, E3, C1, C2, C3, alpha, beta, gamma):
    '''Cell viability based on MuSyC model assuming no 3-way drug interactions.
    alpha must be provided as a 6-long numpy array.
    beta and gamma must be provided as 3-long numpy arrays.'''

    r1, r2, r3 = 1, 1, 1 # arbitrary number
    rn1 = (C1**h1)*r1
    rn2 = (C2**h2)*r1
    rn3 = (C3**h3)*r3

    Ed = np.empty(len(conc1))

    for i in range(len(Ed)):
        d1 = conc1[i]
        d2 = conc2[i]
        d3 = conc3[i]

        E12 = -beta[0]*min(E1,E2)+min(E1,E2)
        E13 = -beta[1]*min(E1,E3)+min(E1,E3)
        E23 = -beta[2]*min(E2,E3)+min(E2,E3)
        E123 = min([E12,E13,E23])

        Ev = np.array([[E0, E1, E2, E3, E12, E13, E23, E123]])
        Y = np.zeros((Ev.shape[1], Ev.shape[1]))

        Y[0,0] = - r1*d1**h1 - r2*d2**h2 - r3*d3**h3
        Y[0,1] = rn1
        Y[0,2] = rn2
        Y[0,3] = rn3

        Y[1,0] = r1*d1**h1
        Y[1,1] = - rn1 - r2*alpha[0]*d2**h2 - r3*alpha[2]*d3**h3
        Y[1,4] = rn2
        Y[1,5] = rn3

        Y[2,0] = r2*d2**h2
        Y[2,2] = - r1*alpha[1]*d1**h1 - rn2 - r3*alpha[4]*d3**h3
        Y[2,4] = rn1
        Y[2,6] = rn3

        Y[3,0] = r3*d3**h3
        Y[3,3] = - r1*alpha[3]*d1**h1 - r2*alpha[5]*d2**h2 - rn3;
        Y[3,5] = rn1
        Y[3,6] = rn2
        
        Y[4,1] = r2*alpha[0]*d2**h2
        Y[4,2] = r1*alpha[1]*d1**h1
        Y[4,4] = - rn1 - rn2 - r3*gamma[0]*d3**h3
        Y[4,7] = rn3
        
        Y[5,1] = r3*alpha[2]*d3**h3
        Y[5,3] = r1*alpha[3]*d1**h1
        Y[5,5] = - rn1 - r2*gamma[1]*d2**h2 - rn3
        Y[5,6] = 1 # is this technically correct? Need to show this
        Y[5,7] = rn2

        Y[6,2] = r3*alpha[4]*d3**h3
        Y[6,3] = r2*alpha[5]*d2**h2
        Y[6,6] = - r1*gamma[2]*d1**h1 - rn2 - rn3
        Y[6,7] = rn1
        
        Y[7,:] = 1

        b = np.array([[0,0,0,0,0,0,0,1]]).T

        Ed[i] = np.matmul(np.matmul(Ev, np.linalg.inv(Y)), b)
    return(Ed)

def musyc_from_excel(filename):
    '''Imports Excel file with simulation conditions and runs synergy analysis'''
    data = pd.read_excel(filename)

    syn12 = np.empty(len(data))
    syn13 = np.empty(len(data))
    syn23 = np.empty(len(data))

    for i in range(len(data)): # for each simulation
                
        dose1 = np.array(eval(data.at[i, 'drug1.dose']))
        dose2 = np.array(eval(data.at[i, 'drug2.dose']))
        dose3 = np.array(eval(data.at[i, 'drug3.dose']))
    
        device_list = []
        num_devices = len(dose1)

        for j in range(num_devices): # for each device in experiment
            input_dose_1 = dose1[j]
            input_dose_2 = dose2[j]
            input_dose_3 = dose3[j]

            num_cells = data.at[i, 'num_cells']

            model = data.at[i, 'model']

            h1 = data.at[i, 'drug1.hill']
            h2 = data.at[i, 'drug2.hill']
            h3 = data.at[i, 'drug3.hill']
            E0 = data.at[i, 'E0']
            E1 = data.at[i, 'drug1.Emax']
            E2 = data.at[i, 'drug2.Emax']
            E3 = data.at[i, 'drug3.Emax']
            C1 = data.at[i, 'drug1.EC50']
            C2 = data.at[i, 'drug2.EC50']
            C3 = data.at[i, 'drug3.EC50']
            alpha12 = data.at[i, 'alpha12']
            alpha13 = data.at[i, 'alpha13']
            alpha23 = data.at[i, 'alpha23']
            beta12 = data.at[i, 'beta12']
            beta13 = data.at[i, 'beta13']
            beta23 = data.at[i, 'beta23']
            
            alpha = np.array([alpha12, alpha12, alpha23, alpha23, alpha13, alpha13])
            beta = np.array([beta12, beta23, beta13])
            gamma = np.array([1, 1, 1])

            concMap = load_map()
            cells = concMap.sample(num_cells)

            conc1 = cells[["c1"]].to_numpy()
            conc2 = cells[["c2"]].to_numpy()
            conc3 = cells[["c3"]].to_numpy()

            if model == "musyc":
                Ed = musyc(conc1*input_dose_1, conc2*input_dose_2, conc3*input_dose_3,\
                    h1, h2, h3, E0, E1, E2, E3, C1, C2, C3, alpha, beta, gamma)
            
            rand = np.random.random_sample((cells.shape[0]))
            live = (rand < Ed)
            cells["live"] = live.astype(int)

            device1 = leukdev.Device(cells, ["drug 1", "drug 2", "drug 3"], [input_dose_1, input_dose_2, input_dose_3], ["nM", "nM", "nM"])
            device1.print_summary()
            device_list.append(device1)

        stack = leukdev.DeviceStack(device_list)
        print(stack)

        syn12[i] = stack.diamond(1,2)
        syn13[i] = stack.diamond(1,3)
        syn23[i] = stack.diamond(2,3)


    # add the syn12, syn13, and syn23 columns to the dataframe if it doesn't exist
    #if not "syn12" in data.columns:
    data["syn12"] = syn12
    data["syn13"] = syn13
    data["syn23"] = syn23

    # replace the dataframe as excel
    data.to_excel(filename, index = False)

    
    print(f"syn12 = {syn12}")
    print(f"syn13 = {syn13}")
    print(f"syn23 = {syn23}")

def greco_from_excel(filename):
    '''Imports Excel file with simulation conditions and runs synergy analysis'''
    data = pd.read_excel(filename)

    syn12 = np.empty(len(data))
    syn13 = np.empty(len(data))
    syn23 = np.empty(len(data))

    for i in range(len(data)): # for each simulation
                
        dose1 = np.array(eval(data.at[i, 'drug1.dose']))
        dose2 = np.array(eval(data.at[i, 'drug2.dose']))
        dose3 = np.array(eval(data.at[i, 'drug3.dose']))
    
        device_list = []
        num_devices = len(dose1)

        for j in range(num_devices): # for each device in experiment
            input_dose_1 = dose1[j]
            input_dose_2 = dose2[j]
            input_dose_3 = dose3[j]

            num_cells = data.at[i, 'num_cells']

            model = data.at[i, 'model']

            h1 = data.at[i, 'drug1.hill']
            h2 = data.at[i, 'drug2.hill']
            h3 = data.at[i, 'drug3.hill']
            E0 = data.at[i, 'E0']
            Emax = data.at[i, 'Emax']
            C1 = data.at[i, 'drug1.EC50']
            C2 = data.at[i, 'drug2.EC50']
            C3 = data.at[i, 'drug3.EC50']
            beta_D12 = data.at[i, 'beta_D12']
            beta_D13 = data.at[i, 'beta_D13']
            beta_D23 = data.at[i, 'beta_D23']
            

            concMap = load_map()
            cells = concMap.sample(num_cells)

            conc1 = cells[["c1"]].to_numpy()
            conc2 = cells[["c2"]].to_numpy()
            conc3 = cells[["c3"]].to_numpy()

            if model == "greco":
                Ed = greco(conc1*input_dose_1, conc2*input_dose_2, conc3*input_dose_3,\
                    beta_D12 = beta_D12, beta_D13 = beta_D13, beta_D23 = beta_D23,\
                    IC50_1 = C1, IC50_2 = C2, IC50_3 = C3,\
                    E0 = E0, Emax = Emax, alpha_m1 = h1, alpha_m2 = h2, alpha_m3 = h3)
            
            rand = np.random.random_sample((cells.shape[0]))
            live = (rand < Ed)
            cells["live"] = live.astype(int)

            device1 = leukdev.Device(cells, ["drug 1", "drug 2", "drug 3"], [input_dose_1, input_dose_2, input_dose_3], ["nM", "nM", "nM"])
            device1.print_summary()
            device_list.append(device1)

        stack = leukdev.DeviceStack(device_list)

        syn12[i] = stack.diamond(1,2)
        syn13[i] = stack.diamond(1,3)
        syn23[i] = stack.diamond(2,3)


    # add the syn12, syn13, and syn23 columns to the dataframe if it doesn't exist
    #if not "syn12" in data.columns:
    data["syn12"] = syn12
    data["syn13"] = syn13
    data["syn23"] = syn23

    # replace the dataframe as excel
    data.to_excel(filename, index = False)


def export_simulation(device_stack, simulation_name):
    working_dir = os.getcwd()
    sim_path = os.path.join(working_dir, simulation_name)
    # check to see if simulation folder exists
    if os.path.exists(sim_path) == False:
        os.mkdir(sim_path)
    else:
        confirm_delete = input(f"Directory '{simulation_name}' already exists. Over write files? (y/n): ")
        if confirm_delete == "y":
            shutil.rmtree(sim_path)
            os.mkdir(sim_path)
        else:
            raise ValueError("Selected 'n' for overwrite files. Recheck directory names and try again.")

    for i in range(device_stack.get_num_devices()):
        data = device_stack.devices[i].dataTable
        filename = simulation_name + "_Dev" + str(i) + ".xlsx"
        path = os.path.join(working_dir, simulation_name, filename)
        data.to_excel(path, index = False, header = False)
    
    # config file
    config_file = configparser.ConfigParser()
    config_file["Drug1"] = {
        "name":device_stack.drug_names[0],
        "doses":device_stack.dose[:,0],
        "units":device_stack.dose_units[0],
    }
    config_file.add_section("Drug2")
    config_file.add_section("Drug3")
    config_file.add_section("Synergy")

    # info.txt
    line1 = f"{device_stack.drug_names[0]}; {device_stack.dose[:,0].tolist()}; {device_stack.dose_units[0]}"
    line2 = f"{device_stack.drug_names[1]}; {device_stack.dose[:,1].tolist()}; {device_stack.dose_units[1]}"
    line3 = f"{device_stack.drug_names[2]}; {device_stack.dose[:,2].tolist()}; {device_stack.dose_units[2]}"
    f = open(os.path.join(sim_path, "info.txt"), 'w')
    f.write(line1)
    f.write('\n')
    f.write(line2)
    f.write('\n')
    f.write(line3)

    print("Simulation results export successful.")


def main(model = "musyc"):
    num_devices = 8
    doses_1 = np.array([400, 200, 100, 50, 25, 12.5, 6.25, 0])
    doses_2 = np.array([400, 200, 100, 50, 25, 12.5, 6.25, 0])
    doses_3 = np.array([400, 200, 100, 50, 25, 12.5, 6.25, 0])

    #generate_config()
    config = import_config()
    #drug1_hill = config['Drug2']['Hill']
    #drug1_hill = config.get('Drug2', 'Hill')

    device_list = [] # arbitrary list
    for i in range(num_devices):
        input_dose_1 = doses_1[i]
        input_dose_2 = doses_2[i]
        input_dose_3 = doses_3[i]

        # make device from this
        concMap = load_map()
        cells = sample_cells(concMap, 6000)
        
        drug1_name = "Daunorubicin"
        drug2_name = "Vincristine"
        drug3_name = "Prednisolone"

        conc1 = cells[["c1"]].to_numpy()
        conc2 = cells[["c2"]].to_numpy()
        conc3 = cells[["c3"]].to_numpy()

        h1 = 2
        h2 = 2
        h3 = 2
        E0 = 1
        E1 = 0.2
        E2 = 0.2
        E3 = 0.2
        C1 = 50
        C2 = 50
        C3 = 50

        alpha12 = 100
        alpha23 = 0
        alpha13 = 1

        beta12 = 0
        beta23 = 0
        beta13 = 0

        alpha = np.array([alpha12, alpha12, alpha23, alpha23, alpha13, alpha13])
        beta = np.array([beta12, beta23, beta13])
        gamma = np.array([1, 1, 1])

        if model == "musyc":
            Ed = musyc(conc1*input_dose_1, conc2*input_dose_2, conc3*input_dose_3,\
                h1, h2, h3, E0, E1, E2, E3, C1, C2, C3, alpha, beta, gamma)
        elif model == "greco":
            Ed = greco(conc1*input_dose_1, conc2*input_dose_2, conc3*input_dose_3, beta_D12 = 0, beta_D13 = 2, beta_D23 = -2)
        rand = np.random.random_sample((cells.shape[0]))
        live = (rand < Ed)

        cells["live"] = live.astype(int)
        
        device1 = leukdev.Device(cells, [drug1_name, drug2_name, drug3_name], [input_dose_1, input_dose_2, input_dose_3], ["nM", "nM", "nM"])
        device1.print_summary()
        device_list.append(device1)

    devicestack1 = leukdev.DeviceStack(device_list)
    #devicestack1.singledrugfits()
    export_simulation(devicestack1, "sim_ant_syn")


def generate_template(filename, model = "musyc"):
    '''Generates Excel template for multi-condition simulations.
    Current choices of models include greco and musyc'''

    columns_musyc = ['simID', 'model', 'description', 'num_cells', 'E0',
        'drug1.Emax', 'drug1.EC50', 'drug1.hill',
        'drug2.Emax', 'drug2.EC50', 'drug2.hill',
        'drug3.Emax', 'drug3.EC50', 'drug3.hill',
        'drug1.dose', 'drug2.dose', 'drug3.dose',
        'alpha12', 'alpha13', 'alpha23',
        'beta12', 'beta13', 'beta23']

    columns_greco = ['simID', 'model', 'description', 'num_cells', 'E0',
        'Emax', 'drug1.EC50', 'drug1.hill',
        'drug2.EC50', 'drug2.hill',
        'drug3.EC50', 'drug3.hill',
        'drug1.dose', 'drug2.dose', 'drug3.dose',
        'beta_D12', 'beta_D13', 'beta_D23']


    example_data_musyc = [[1, model, "example", 1000, 1,
        0, 1, 2, 0, 1, 2, 0, 1, 2,
        "[0, 1, 2, 3]", "[0, 1, 2, 3]", "[0, 1, 2, 3]",
        0, 0, 0, 0, 0, 0]]

    example_data_greco = [[1, model, "example", 1000, 1, 0,
        1, 2, 1, 2, 1, 2,
        "[0, 1, 2, 3]", "[0, 1, 2, 3]", "[0, 1, 2, 3]",
        0, 0, 0]]

    match model:
        case "musyc":
            data = pd.DataFrame(example_data_musyc, columns = columns_musyc)
        case "greco":
            data = pd.DataFrame(example_data_greco, columns = columns_greco)

    
    data.to_excel(filename + ".xlsx", index = False)
    print(f"Template saved to {filename}.xlsx")

def sim_expand():
    # expands the number of simulations by conditions
    # for given number of replicates

    filename = "conditions_greco_sim.xlsx"
    reps = 10

    df = pd.read_excel(filename)

    col_names = ["beta_D12", "beta_D13", "beta_D23"]
    #data = [[0, 0, 0]]

    data = np.empty([reps*len(df), 3])
    conditions = []

    n = 0
    for i in range(len(df)):
        for j in range(reps):
            condition1 = df.at[i, "drug12"]
            condition2 = df.at[i, "drug13"]
            condition3 = df.at[i, "drug23"]

            conditions.append(str(condition1) + "_" + str(condition2) + "_" + str(condition3))
            
            match condition1:
                case "#": score = 0
                case "*": score = -2
                case "**": score = -4
                case "!": score = 2
                case "!!": score = 4
            data[n, 0] = score

            match condition2:
                case "#": score = 0
                case "*": score = -2
                case "**": score = -4
                case "!": score = 2
                case "!!": score = 4
            data[n, 1] = score

            match condition3:
                case "#": score = 0
                case "*": score = -2
                case "**": score = -4
                case "!": score = 2
                case "!!": score = 4
            data[n, 2] = score
            n += 1


    output = pd.DataFrame(data, columns = col_names)
    output["condition"] = conditions

    output.to_excel("expanded.xlsx", index = False)
    print(df)

    print(output)


def sim_expand_analysis():
    '''Expands the simulations and computes p values for average synergy'''
    condition_list = pd.read_excel("conditions_greco_sim.xlsx")
    results = pd.read_excel("greco_study_expand.xlsx")

    # preallocate
    f_value = np.empty(len(condition_list))
    p_value = np.empty(len(condition_list))
    syn12_avg = np.empty(len(condition_list))
    syn13_avg = np.empty(len(condition_list))
    syn23_avg = np.empty(len(condition_list))

    for i, condition in enumerate(condition_list["condition"]):
        subset = results.loc[results["description"]==condition]
        f_value[i], p_value[i] = f_oneway(subset["syn12"], subset["syn13"], subset["syn23"])
        syn12_avg[i] = subset["syn12"].mean()
        syn13_avg[i] = subset["syn13"].mean()
        syn23_avg[i] = subset["syn23"].mean()

        transformed = subset[["syn12", "syn13", "syn23"]]
        transformed = pd.melt(transformed, value_vars = ["syn12", "syn13", "syn23"])
        # results in 2 columns: variable and value
        #tukey = pairwise_tukeyhsd(endog = transformed["value"], groups = transformed["variable"], alpha = 0.05)
        # tukey takes a long time...

    condition_list["F-statistic"] = f_value
    condition_list["P-value"] = p_value
    condition_list["syn12_avg"] = syn12_avg
    condition_list["syn13_avg"] = syn13_avg
    condition_list["syn23_avg"] = syn23_avg

    print(condition_list)
    condition_list.to_excel("condition_stats2.xlsx", index = False)


if __name__ == "__main__":
    #main("musyc")
    #main()
    #musyc_from_excel("test.xlsx")
    #greco_from_excel("greco_study_expand.xlsx")
    #musyc_from_excel("musyc_dvn.xlsx")
    greco_from_excel("dvp_equivalency.xlsx")

