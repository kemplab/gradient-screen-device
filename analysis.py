import leukemiadrugscreen as leukdev
import os
import sys
import EC50_tern
import numpy as np
from synergy.single.hill import Hill
import ternary
from ternary_diagram import TernaryDiagram
import matplotlib.pyplot as plt
import pandas as pd

def load_folder(folder = "empty"):
    if folder == "empty":
        folder = input("Type name of folder in same directory, or path: ")
    
    working_dir = os.getcwd()
    folder_path = os.path.join(working_dir, folder)
    if os.path.exists(folder_path) == False:
        raise ImportError(f"No folder with name '{folder}' exists.")

    return leukdev.DeviceStack.from_folder(folder)

def main():
    stack = load_folder()
    print(stack)
    #print(stack.single_drug_fit(1, False))
    #print(stack.single_drug_fit(2, False))
    #print(stack.single_drug_fit(3, False))
    print(f"{stack.diamond(1,2):.2f}")
    print(f"{stack.diamond(1,3):.2f}")
    print(f"{stack.diamond(2,3):.2f}")


def bliss_test():
    stack = load_folder("rep1")
    print(stack)
    


def main2():
    stack = load_folder()
    print(stack.devices[0].get_zone(11, True))

def diamond3_test():
    stack = load_folder("rep2")
    print(stack)
    FIC2_xy = stack.diamond(1,2, plot=True)
    FIC2_xz = stack.diamond(1,3, plot = True)
    FIC2_yz = stack.diamond(2,3, plot = True)

    print(f"FIC2_xy = {FIC2_xy:.2f}")
    print(f"FIC2_xz = {FIC2_xz:.2f}")
    print(f"FIC2_yz = {FIC2_yz:.2f}")

    #stack.drug3_fit(0.4, 0.2, 0.3)

def hexagon_zone(foldername):
    '''Plots ternary plots for EC50 and Emax with hexagonal zones'''
    from drugfit import Hill

    stack = load_folder(foldername)
    
    tern_grid = EC50_tern.gen_tern_grid(11)
    
    viability = np.empty([len(tern_grid), stack.get_num_devices()])
    num_cells = np.empty_like(viability)

    for i in range(len(tern_grid)):
        for x in range(stack.get_num_devices()):
            tern_point = (tern_grid["a"][i], tern_grid["b"][i], tern_grid["c"][i])
            zone = leukdev.Zone.from_point(tern_point, stack.devices[x])
            viability[i, x] = zone.cells["live"].mean()
            num_cells[i, x] = len(zone.cells)

    # fit and normalize to EC50
    EC50_A = stack.singledrugfits[0].C
    EC50_B = stack.singledrugfits[1].C
    EC50_C = stack.singledrugfits[2].C

    print(f"EC50_A = {EC50_A:.2f}")
    print(f"EC50_B = {EC50_B:.2f}")
    print(f"EC50_C = {EC50_C:.2f}")

    EC50 = np.empty(len(tern_grid))
    Emax = np.empty_like(EC50)

    dose_inputs_A = stack.dose[:,0].transpose()
    dose_inputs_B = stack.dose[:,1].transpose()
    dose_inputs_C = stack.dose[:,2].transpose()

    for i in range(len(EC50)):
        dose_A = tern_grid["a"][i] * dose_inputs_A / EC50_A
        dose_B = tern_grid["b"][i] * dose_inputs_B / EC50_B
        dose_C = tern_grid["c"][i] * dose_inputs_C / EC50_C
        dose_combo = dose_A + dose_B + dose_C
        
        x = dose_combo
        y = viability[i, :]

        if (np.isnan(np.sum(x)) == True or np.isnan(np.sum(y)) == True):
            print("Found NaN")
            EC50[i] = np.nan
            Emax[i] = np.nan
            continue

        hill_model = Hill.fit(x, y)
        EC50[i] = hill_model.C
        Emax[i] = hill_model.Emax # E0, Emax, hill, EC50

    # filter out bad EC50
    nan_idx = np.isnan(Emax)
    EC50[nan_idx] = 1

    # set each edge to where nan is to single_drug_emax
    drug1_grid_idx = tern_grid[(tern_grid["a"]>=0.8) & (tern_grid["b"]<0.2)].index.to_list()
    drug2_grid_idx = tern_grid[(tern_grid["b"]>=0.8) & (tern_grid["a"]<0.2)].index.to_list()
    drug3_grid_idx = tern_grid[(tern_grid["c"]>=0.8) & (tern_grid["a"]<0.2)].index.to_list()
    
    #Emax[drug1_grid_idx] = stack.singledrugfits[0].Emax
    
    Emax[drug1_grid_idx] = stack.singledrugfits[0].Emax
    Emax[drug2_grid_idx] = stack.singledrugfits[1].Emax
    Emax[drug3_grid_idx] = stack.singledrugfits[2].Emax

    # log scale the EC50
    EC50 = np.log2(EC50)

    # plotting EC50
    plot_data = dict()
    for x in range(len(tern_grid)):
        i = int(tern_grid["a"][x]*10)
        j = int(tern_grid["b"][x]*10)
        plot_data[(i,j)] = EC50[x]
    


    figure, tax = ternary.figure(scale = 10)
    #tax.gridlines(color = "black", multiple = 1)
    tax.heatmap(plot_data, style = "hexagonal", cmap = "coolwarm_r", vmax = 2, vmin = -2)
    #tax.heatmap(plot_data, style = "hexagonal", cmap = "coolwarm_r", vmax = 2, vmin = 0)
    fontsize = 12
    #tax.right_corner_label(stack.drug_names[0], fontsize = fontsize)
    #tax.top_corner_label(stack.drug_names[1], fontsize = fontsize)
    #tax.left_corner_label(stack.drug_names[2], fontsize = fontsize)

    tax.gridlines(color = "grey", multiple = 1)
    tax.get_axes().axis('off')
    #tax.ticks(axis = 'lbr', multiple = 1, linewidth = 1, offset = 0.025)
    tax.clear_matplotlib_ticks()
    tax.boundary(linewidth = 2.0)
    
    tax.show()


    # plotting Emax
    plot_data = dict()
    for x in range(len(tern_grid)):
        i = int(tern_grid["a"][x]*10)
        j = int(tern_grid["b"][x]*10)
        plot_data[(i,j)] = Emax[x]
    


    figure, tax = ternary.figure(scale = 10)
    #tax.gridlines(color = "black", multiple = 1)
    tax.heatmap(plot_data, style = "hexagonal", cmap = "viridis_r", vmax = 1, vmin = 0)
    #tax.heatmap(plot_data, style = "hexagonal", cmap = "coolwarm_r", vmax = 2, vmin = 0)
    fontsize = 12
    #tax.right_corner_label(stack.drug_names[0], fontsize = fontsize)
    #tax.top_corner_label(stack.drug_names[1], fontsize = fontsize)
    #tax.left_corner_label(stack.drug_names[2], fontsize = fontsize)

    tax.gridlines(color = "grey", multiple = 1)
    tax.get_axes().axis('off')
    #tax.ticks(axis = 'lbr', multiple = 1, linewidth = 1, offset = 0.025)
    tax.clear_matplotlib_ticks()
    tax.boundary(linewidth = 2.0)
    
    tax.show()

    # alternate plotting package, maybe prettier?
    tern_grid["z"] = EC50
    labels = tern_grid.columns[0:3].tolist()

if __name__ == "__main__":
    #diamond3_test()
    hexagon_zone("DVP_Rep1")