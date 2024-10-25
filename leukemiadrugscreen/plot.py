import numpy as np
import ternary
from . import ternary_tools
import leukemiadrugscreen as leukdev
from .drugfit import Hill
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt

def _color_list():
    colors = ["tab:red", "tab:green", "tab:blue", "y", "m", "c", "#a0b1ba", "#a6761d"]
    return colors

def _color(str):
    # get matplotlib color code based on string
    colors = _color_list()

    match str:
        case "drug1":
            color_code = colors[0]
        case "drug2":
            color_code = colors[1]
        case "drug3":
            color_code = colors[2]
        case "combo12":
            color_code = colors[3]
        case "combo13":
            color_code = colors[4]
        case "combo23":
            color_code = colors[5]
        case "combo123":
            color_code = colors[6]
    
    return color_code

def ternary_plots(stack):
    '''Plots ternary plots for EC50 and Emax with hexagonal zones'''
    # this should be changed to utilize device stack instead
    # maybe store in the ternary tools page instead
    
    tern_grid = ternary_tools.gen_tern_grid(11)
    
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
        y=y/y[-1]  #normalize to control device viability

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
    fontsize = 11
    tax.right_corner_label(stack.abbrevs[0], fontsize = fontsize)
    tax.top_corner_label(stack.abbrevs[1], fontsize = fontsize)
    tax.left_corner_label(stack.abbrevs[2], fontsize = fontsize)

    tax.gridlines(color = "grey", multiple = 1)
    tax.get_axes().axis('off')
    #tax.ticks(axis = 'lbr', multiple = 1, linewidth = 1, offset = 0.025)
    tax.clear_matplotlib_ticks()
    tax.boundary(linewidth = 2.0)
    
    print('Normalized EC50 Ternary Plot')
    
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
    fontsize = 11
    tax.right_corner_label(stack.abbrevs[0], fontsize = fontsize)
    tax.top_corner_label(stack.abbrevs[1], fontsize = fontsize)
    tax.left_corner_label(stack.abbrevs[2], fontsize = fontsize)

    tax.gridlines(color = "grey", multiple = 1)
    tax.get_axes().axis('off')
    #tax.ticks(axis = 'lbr', multiple = 1, linewidth = 1, offset = 0.025)
    tax.clear_matplotlib_ticks()
    tax.boundary(linewidth = 2.0)
    
    print('Emax Ternary Plot')
    tax.show()

    #plotting beta
    Emax_A = stack.singledrugfits[0].Emax
    Emax_B = stack.singledrugfits[1].Emax
    Emax_C = stack.singledrugfits[2].Emax
    plot_data = dict()
    for x in range(len(tern_grid)):
        i = int(tern_grid["a"][x]*10)
        j = int(tern_grid["b"][x]*10)
        plot_data[(i,j)] = ((min(Emax_A,Emax_B, Emax_C)- Emax[x])/(1-min(Emax_A,Emax_B, Emax_C)))
        
    figure, tax = ternary.figure(scale = 10)
    #tax.gridlines(color = "black", multiple = 1)
    #tax.heatmap(plot_data, style = "hexagonal", cmap = "viridis_r", vmax = 1, vmin = 0)
    tax.heatmap(plot_data, style = "hexagonal", cmap = "coolwarm", vmax = 0.5, vmin = -0.5)
    fontsize = 11
    tax.right_corner_label(stack.abbrevs[0], fontsize = fontsize)
    tax.top_corner_label(stack.abbrevs[1], fontsize = fontsize)
    tax.left_corner_label(stack.abbrevs[2], fontsize = fontsize)

    tax.gridlines(color = "grey", multiple = 1)
    tax.get_axes().axis('off')
    #tax.ticks(axis = 'lbr', multiple = 1, linewidth = 1, offset = 0.025)
    tax.clear_matplotlib_ticks()
    tax.boundary(linewidth = 2.0)
    
    print('Beta Ternary Plot')
    tax.show()

    # alternate plotting package, maybe prettier?
    tern_grid["z"] = EC50
    labels = tern_grid.columns[0:3].tolist()

def ternary_plots_GR(stack):
    '''Plots ternary plots for EC50 and Emax with hexagonal zones'''
    # this should be changed to utilize device stack instead
    # maybe store in the ternary tools page instead
    
    tern_grid = ternary_tools.gen_tern_grid(11)
    
    GR = np.empty([len(tern_grid), stack.get_num_devices()])
    num_cells = np.empty_like(GR)

    for i in range(len(tern_grid)):
        tern_point = (tern_grid["a"][i], tern_grid["b"][i], tern_grid["c"][i])
        region = leukdev.Region.from_point(tern_point, stack)
        GR[i, :] = region.get_GR()

    # fit and normalize to EC50
    GR50_A = stack.singledrugfitsGR[0].C
    GR50_B = stack.singledrugfitsGR[1].C
    GR50_C = stack.singledrugfitsGR[2].C

    print(f"GR50_A = {GR50_A:.2f}")
    print(f"GR50_B = {GR50_B:.2f}")
    print(f"GR50_C = {GR50_C:.2f}")

    GR50 = np.empty(len(tern_grid))
    GRmax = np.empty_like(GR50)

    dose_inputs_A = stack.dose[:,0].transpose()
    dose_inputs_B = stack.dose[:,1].transpose()
    dose_inputs_C = stack.dose[:,2].transpose()

    for i in range(len(GR50)):
        dose_A = tern_grid["a"][i] * dose_inputs_A / GR50_A
        dose_B = tern_grid["b"][i] * dose_inputs_B / GR50_B
        dose_C = tern_grid["c"][i] * dose_inputs_C / GR50_C
        dose_combo = dose_A + dose_B + dose_C
        
        x = dose_combo
        y = GR[i, :]

        if (np.isnan(np.sum(x)) == True or np.isnan(np.sum(y)) == True):
            print("Found NaN")
            GR50[i] = np.nan
            GRmax[i] = np.nan
            continue

        hill_model = Hill.fit(x, y)
        GR50[i] = hill_model.C
        GRmax[i] = hill_model.Emax # E0, Emax, hill, EC50

    # filter out bad EC50
    nan_idx = np.isnan(GRmax)
    GR50[nan_idx] = 1

    # set each edge to where nan is to single_drug_emax
    drug1_grid_idx = tern_grid[(tern_grid["a"]>=0.8) & (tern_grid["b"]<0.2)].index.to_list()
    drug2_grid_idx = tern_grid[(tern_grid["b"]>=0.8) & (tern_grid["a"]<0.2)].index.to_list()
    drug3_grid_idx = tern_grid[(tern_grid["c"]>=0.8) & (tern_grid["a"]<0.2)].index.to_list()
    
    #Emax[drug1_grid_idx] = stack.singledrugfits[0].Emax
    
    GRmax[drug1_grid_idx] = stack.singledrugfitsGR[0].Emax
    GRmax[drug2_grid_idx] = stack.singledrugfitsGR[1].Emax
    GRmax[drug3_grid_idx] = stack.singledrugfitsGR[2].Emax

    # log scale the EC50
    GR50 = np.log2(GR50)

    # plotting EC50
    plot_data = dict()
    for x in range(len(tern_grid)):
        i = int(tern_grid["a"][x]*10)
        j = int(tern_grid["b"][x]*10)
        plot_data[(i,j)] = GR50[x]
    


    figure, tax = ternary.figure(scale = 10)
    #tax.gridlines(color = "black", multiple = 1)
    tax.heatmap(plot_data, style = "hexagonal", cmap = "coolwarm_r", vmax = 2, vmin = -2)
    #tax.heatmap(plot_data, style = "hexagonal", cmap = "coolwarm_r", vmax = 2, vmin = 0)
    fontsize = 11
    tax.right_corner_label(stack.abbrevs[0], fontsize = fontsize)
    tax.top_corner_label(stack.abbrevs[1], fontsize = fontsize)
    tax.left_corner_label(stack.abbrevs[2], fontsize = fontsize)

    tax.gridlines(color = "grey", multiple = 1)
    tax.get_axes().axis('off')
    #tax.ticks(axis = 'lbr', multiple = 1, linewidth = 1, offset = 0.025)
    tax.clear_matplotlib_ticks()
    tax.boundary(linewidth = 2.0)
    
    print('GR50 Ternary Plot')
    tax.show()


    # plotting Emax
    plot_data = dict()
    for x in range(len(tern_grid)):
        i = int(tern_grid["a"][x]*10)
        j = int(tern_grid["b"][x]*10)
        plot_data[(i,j)] = GRmax[x]
    


    figure, tax = ternary.figure(scale = 10)
    #tax.gridlines(color = "black", multiple = 1)
    tax.heatmap(plot_data, style = "hexagonal", cmap = "viridis_r", vmax = 1, vmin = -1)
    #tax.heatmap(plot_data, style = "hexagonal", cmap = "coolwarm_r", vmax = 2, vmin = 0)
    fontsize = 11
    tax.right_corner_label(stack.abbrevs[0], fontsize = fontsize)
    tax.top_corner_label(stack.abbrevs[1], fontsize = fontsize)
    tax.left_corner_label(stack.abbrevs[2], fontsize = fontsize)

    tax.gridlines(color = "grey", multiple = 1)
    tax.get_axes().axis('off')
    #tax.ticks(axis = 'lbr', multiple = 1, linewidth = 1, offset = 0.025)
    tax.clear_matplotlib_ticks()
    tax.boundary(linewidth = 2.0)
    
    print('GRmax Ternary Plot')
    tax.show()



def zone_plots(stack):
    '''Plots the single drug, equipotent, and 3-drug regions of a DeviceStack'''
    
    # set regions
    d12_point = (0.5, 0.5, 0)
    d13_point = (0.5, 0, 0.5)
    d23_point = (0, 0.5, 0.5)

    region12 = leukdev.Region.from_point(d12_point, stack)
    region13 = leukdev.Region.from_point(d13_point, stack)
    region23 = leukdev.Region.from_point(d23_point, stack)


    d123_point = (1/3, 1/3, 1/3)
    region123 = leukdev.Region.from_point(d123_point, stack)
    
    fig = plt.figure(constrained_layout = True, figsize = (14, 6))
    gs = GridSpec(3, 7, figure = fig)

    colors = _color_list()


    ax1 = fig.add_subplot(gs[:, :-4])
    ax1.set_aspect('equal')
    ax2 = fig.add_subplot(gs[0, 3:5])
    ax3 = fig.add_subplot(gs[1, 3:5])
    ax4 = fig.add_subplot(gs[2, 3:5])

    ax5 = fig.add_subplot(gs[0, 5:7])
    ax6 = fig.add_subplot(gs[1, 5:7])
    ax7 = fig.add_subplot(gs[2, 5:7])

    def plot_tern_bg(zone, ax, color):
        ternx, terny = leukdev.ternary_tools.terncoords(zone.device.dataTable[["c1", "c2", "c3"]])
        ternx_sel, terny_sel = leukdev.ternary_tools.terncoords(zone.cells[["c1", "c2", "c3"]])
        ax.scatter(ternx, terny, c = color) 

    def plot_zone_tern(zone, ax, color):
        ternx, terny = leukdev.ternary_tools.terncoords(zone.device.dataTable[["c1", "c2", "c3"]])
        ternx_sel, terny_sel = leukdev.ternary_tools.terncoords(zone.cells[["c1", "c2", "c3"]])
        #ax.scatter(ternx, terny, c = '0.8')
        ax.scatter(ternx_sel, terny_sel, c = color)

    plot_tern_bg(stack.singledrugregions[0].zones[0], ax1, "0.8")

    plot_zone_tern(stack.singledrugregions[0].zones[0], ax1, colors[0])
    plot_zone_tern(stack.singledrugregions[1].zones[0], ax1, colors[1])
    plot_zone_tern(stack.singledrugregions[2].zones[0], ax1, colors[2])

    plot_zone_tern(region12.zones[0], ax1, colors[3])
    plot_zone_tern(region13.zones[0], ax1, colors[4])
    plot_zone_tern(region23.zones[0], ax1, colors[5])

    plot_zone_tern(region123.zones[0], ax1, colors[6])

    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.spines.clear()

    # plot single drug fits
    def plot_single_drug(drug_idx, ax, color):
        y = stack.singledrugfits[drug_idx-1].y
        x = stack.singledrugfits[drug_idx-1].x
        hillfit = stack.singledrugfits[drug_idx-1]
        xplot = np.linspace(0, max(x), 100)
        yplot = hillfit.E(xplot)
        ax.scatter(x, y, c = color)
        ax.plot(xplot, yplot, c = color)
        ax.set_ylim(0, 1)
        ax.text(max(x)-max(x)/4, 0.8+0.05, stack.drug_names[drug_idx-1])
        ax.text(max(x)-max(x)/4, 0.65+0.05, f"EC50: {hillfit.C:.2f}")
        ax.text(max(x)-max(x)/4, 0.5+0.05, f"Emax: {hillfit.Emax:.2f}")
        #ax.text(max(x)-max(x)/4, 0.5-0.10, f"hill: {hillfit.h:.2f}")
        ax.set_xlabel("Concentration")
        ax.set_ylabel("Viability")

    plot_single_drug(1, ax2, colors[0])
    plot_single_drug(2, ax3, colors[1])
    plot_single_drug(3, ax4, colors[2])


    # ratios and normalization
    dose1 = stack.dose[:,0]
    dose2 = stack.dose[:,1]
    dose3 = stack.dose[:,2]
    EC50_1 = stack.singledrugfits[0].C
    EC50_2 = stack.singledrugfits[1].C
    EC50_3 = stack.singledrugfits[2].C

    print(dose1)
    print(dose2)
    print(dose3)
    print(f"EC50_1: = {EC50_1:.2f}, EC50_2: = {EC50_2:.2f}, EC50_3: = {EC50_3:.2f}")

    dose1_norm = (dose1/EC50_1)
    dose2_norm = (dose2/EC50_2)
    dose3_norm = (dose3/EC50_3)

    ratio12 = dose1_norm[1]/dose2_norm[1]
    print(f"Ratio 1:2 = {ratio12:.2f}")
    ratio13 = dose1_norm[1]/dose3_norm[1]
    print(f"Ratio 1:3 = {ratio13:.2f}")
    ratio23 = dose2_norm[1]/dose3_norm[1]
    print(f"Ratio 2:3 = {ratio23:.2f}")

    def plot_pair_combo(region, drug_idx1, drug_idx2, ax, color):
        y = region.get_viability()
        y= y/y[-1] #normalize to control device viability
        x1 = region.stack.dose[:, drug_idx1-1] / 2
        x2 = region.stack.dose[:, drug_idx2-1] / 2
        EC50 = [EC50_1, EC50_2, EC50_3]
        x = (x1/EC50[drug_idx1-1]) + (x2/EC50[drug_idx2-1])
        hillfit = Hill.fit(x, y)
        xplot = np.linspace(0, max(x), 100)
        yplot = hillfit.E(xplot)
        ax.scatter(x, y, c = color)
        ax.plot(xplot, yplot, c = color)
        ax.set_ylim(0, 1)
        ax.text(max(x)-max(x)/4, 0.8+0.05, region.stack.drug_names[drug_idx1-1][0:4] + "+" + region.stack.drug_names[drug_idx2-1][0:4])
        ax.text(max(x)-max(x)/4, 0.65+0.05, f"EC50: {hillfit.C:.2f}")
        ax.text(max(x)-max(x)/4, 0.5+0.05, f"Emax: {hillfit.Emax:.2f}")
        #ax.text(max(x)-max(x)/4, 0.5-0.1, f"hill: {hillfit.h:.2f}")
        ax.set_xlabel("Normalized EC50")
        ax.set_ylabel("Viability")

    plot_pair_combo(region12, 1, 2, ax5, "y")
    plot_pair_combo(region13, 1, 3, ax6, "c")
    plot_pair_combo(region23, 2, 3, ax7, "m")
    plt.show()

    fig2 = plt.figure(constrained_layout = True, figsize = (4, 3))
    ax = fig2.add_subplot(111)
    y = region123.get_viability()
    y = y/y[-1] #normalize to control device viability
    x1 = region123.stack.dose[:, 0]/3
    x2 = region123.stack.dose[:, 1]/3
    x3 = region123.stack.dose[:, 2]/3
    x = (x1/EC50_1) + (x2/EC50_2) + (x3/EC50_3)
    hillfit = Hill.fit(x, y)
    xplot = np.linspace(0, max(x), 100)
    yplot = hillfit.E(xplot)
    ax.scatter(x, y, c = colors[7])
    ax.plot(xplot, yplot, c = colors[7])
    ax.set_ylim(0, 1)

    text = region123.stack.drug_names[0][0] + "+" + region123.stack.drug_names[1][0] + "+" + region123.stack.drug_names[2][0] + "\n" +\
        f"EC50: {hillfit.C:.2f}\n" +\
        f"Emax: {hillfit.Emax:.2f}\n"
    ax.text(max(x)-max(x)/4, 0.7, text)

    ax.set_xlabel("Normalized EC50")
    ax.set_ylabel("Viability")
    plt.show()

def zone_plots_GR(stack):
    '''Plots the single drug, equipotent, and 3-drug regions of a DeviceStack'''
    
    # set regions
    d12_point = (0.5, 0.5, 0)
    d13_point = (0.5, 0, 0.5)
    d23_point = (0, 0.5, 0.5)

    region12 = leukdev.Region.from_point(d12_point, stack)
    region13 = leukdev.Region.from_point(d13_point, stack)
    region23 = leukdev.Region.from_point(d23_point, stack)


    d123_point = (1/3, 1/3, 1/3)
    region123 = leukdev.Region.from_point(d123_point, stack)
    
    fig = plt.figure(constrained_layout = True, figsize = (14, 6))
    gs = GridSpec(3, 7, figure = fig)

    colors = _color_list()


    ax1 = fig.add_subplot(gs[:, :-4])
    ax1.set_aspect('equal')
    ax2 = fig.add_subplot(gs[0, 3:5])
    ax3 = fig.add_subplot(gs[1, 3:5])
    ax4 = fig.add_subplot(gs[2, 3:5])

    ax5 = fig.add_subplot(gs[0, 5:7])
    ax6 = fig.add_subplot(gs[1, 5:7])
    ax7 = fig.add_subplot(gs[2, 5:7])

    def plot_tern_bg(zone, ax, color):
        ternx, terny = leukdev.ternary_tools.terncoords(zone.device.dataTable[["c1", "c2", "c3"]])
        ternx_sel, terny_sel = leukdev.ternary_tools.terncoords(zone.cells[["c1", "c2", "c3"]])
        ax.scatter(ternx, terny, c = color) 

    def plot_zone_tern(zone, ax, color):
        ternx, terny = leukdev.ternary_tools.terncoords(zone.device.dataTable[["c1", "c2", "c3"]])
        ternx_sel, terny_sel = leukdev.ternary_tools.terncoords(zone.cells[["c1", "c2", "c3"]])
        #ax.scatter(ternx, terny, c = '0.8')
        ax.scatter(ternx_sel, terny_sel, c = color)

    plot_tern_bg(stack.singledrugregions[0].zones[0], ax1, "0.8")

    plot_zone_tern(stack.singledrugregions[0].zones[0], ax1, colors[0])
    plot_zone_tern(stack.singledrugregions[1].zones[0], ax1, colors[1])
    plot_zone_tern(stack.singledrugregions[2].zones[0], ax1, colors[2])

    plot_zone_tern(region12.zones[0], ax1, colors[3])
    plot_zone_tern(region13.zones[0], ax1, colors[4])
    plot_zone_tern(region23.zones[0], ax1, colors[5])

    plot_zone_tern(region123.zones[0], ax1, colors[6])

    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.spines.clear()

    # plot single drug fits
    def plot_single_drug(drug_idx, ax, color):
        y = stack.singledrugfitsGR[drug_idx-1].y
        x = stack.singledrugfitsGR[drug_idx-1].x
        hillfit = stack.singledrugfitsGR[drug_idx-1]
        xplot = np.linspace(0, max(x), 100)
        yplot = hillfit.E(xplot)
        ax.scatter(x, y, c = color)
        ax.plot(xplot, yplot, c = color)
        ax.set_ylim(-1, 1)
        ax.text(max(x)-max(x)/4, 0.6+0.05, stack.drug_names[drug_idx-1])
        ax.text(max(x)-max(x)/4, 0.3+0.05, f"GR50: {hillfit.C:.2f}")
        ax.text(max(x)-max(x)/4, 0.0+0.05, f"GRmax: {hillfit.Emax:.2f}")
        #ax.text(max(x)-max(x)/4, 0.5-0.10, f"hill: {hillfit.h:.2f}")
        ax.set_xlabel("Concentration")
        ax.set_ylabel("Viability")

    plot_single_drug(1, ax2, colors[0])
    plot_single_drug(2, ax3, colors[1])
    plot_single_drug(3, ax4, colors[2])


    # ratios and normalization
    dose1 = stack.dose[:,0]
    dose2 = stack.dose[:,1]
    dose3 = stack.dose[:,2]
    GR50_1 = stack.singledrugfitsGR[0].C
    GR50_2 = stack.singledrugfitsGR[1].C
    GR50_3 = stack.singledrugfitsGR[2].C

    print(dose1)
    print(dose2)
    print(dose3)
    print(f"GR50_1: = {GR50_1:.2f}, GR50_2: = {GR50_2:.2f}, GR50_3: = {GR50_3:.2f}")

    dose1_norm = (dose1/GR50_1)
    dose2_norm = (dose2/GR50_2)
    dose3_norm = (dose3/GR50_3)

    ratio12 = dose1_norm[1]/dose2_norm[1]
    print(f"Ratio 1:2 = {ratio12:.2f}")
    ratio13 = dose1_norm[1]/dose3_norm[1]
    print(f"Ratio 1:3 = {ratio13:.2f}")
    ratio23 = dose2_norm[1]/dose3_norm[1]
    print(f"Ratio 2:3 = {ratio23:.2f}")

    def plot_pair_combo(region, drug_idx1, drug_idx2, ax, color):
        y = region.get_viability()
        x1 = region.stack.dose[:, drug_idx1-1] / 2
        x2 = region.stack.dose[:, drug_idx2-1] / 2
        GR50 = [GR50_1, GR50_2, GR50_3]
        x = (x1/GR50[drug_idx1-1]) + (x2/GR50[drug_idx2-1])
        hillfit = Hill.fit(x, y)
        xplot = np.linspace(0, max(x), 100)
        yplot = hillfit.E(xplot)
        ax.scatter(x, y, c = color)
        ax.plot(xplot, yplot, c = color)
        ax.set_ylim(-1, 1)
        ax.text(max(x)-max(x)/4, 0.6+0.05, region.stack.drug_names[drug_idx1-1][0:4] + "+" + region.stack.drug_names[drug_idx2-1][0:4])
        ax.text(max(x)-max(x)/4, 0.3+0.05, f"GR50: {hillfit.C:.2f}")
        ax.text(max(x)-max(x)/4, 0.0+0.05, f"GRmax: {hillfit.Emax:.2f}")
        #ax.text(max(x)-max(x)/4, 0.5-0.1, f"hill: {hillfit.h:.2f}")
        ax.set_xlabel("Normalized GR50")
        ax.set_ylabel("GR Inhibition")

    plot_pair_combo(region12, 1, 2, ax5, "y")
    plot_pair_combo(region13, 1, 3, ax6, "c")
    plot_pair_combo(region23, 2, 3, ax7, "m")
    plt.show()

    fig2 = plt.figure(constrained_layout = True, figsize = (4, 3))
    ax = fig2.add_subplot(111)
    y = region123.get_viability()
    x1 = region123.stack.dose[:, 0]/3
    x2 = region123.stack.dose[:, 1]/3
    x3 = region123.stack.dose[:, 2]/3
    x = (x1/GR50_1) + (x2/GR50_2) + (x3/GR50_3)
    hillfit = Hill.fit(x, y)
    xplot = np.linspace(0, max(x), 100)
    yplot = hillfit.E(xplot)
    ax.scatter(x, y, c = colors[7])
    ax.plot(xplot, yplot, c = colors[7])
    ax.set_ylim(-1, 1)

    text = region123.stack.drug_names[0][0] + "+" + region123.stack.drug_names[1][0] + "+" + region123.stack.drug_names[2][0] + "\n" +\
        f"GR50: {hillfit.C:.2f}\n" +\
        f"GRmax: {hillfit.Emax:.2f}\n"
    ax.text(max(x)-max(x)/4, 0.5, text)

    ax.set_xlabel("Normalized GR50")
    ax.set_ylabel("GR Inhibition")
    plt.show()



def emax_plot(stack):
    # set regions
    d12_point = (0.5, 0.5, 0)
    d13_point = (0.5, 0, 0.5)
    d23_point = (0, 0.5, 0.5)

    region12 = leukdev.Region.from_point(d12_point, stack)
    region13 = leukdev.Region.from_point(d13_point, stack)
    region23 = leukdev.Region.from_point(d23_point, stack)


    d123_point = (1/3, 1/3, 1/3)
    region123 = leukdev.Region.from_point(d123_point, stack)
    
    
    # beta calculations
    E0 = 1
    E1 = stack.singledrugfits[0].Emax
    E2 = stack.singledrugfits[1].Emax
    E3 = stack.singledrugfits[2].Emax
    E12 = region12.fit_hill().Emax
    E13 = region13.fit_hill().Emax
    E23 = region23.fit_hill().Emax
    E123 = region123.fit_hill().Emax

    def beta(E1, E2, E3):
        E0 = 1
        beta = (min(E1, E2) - E3) / (E0 - min(E1, E2))
        return beta

    def beta3(E1, E2, E3, E123):
        E0 = 1
        strongest = min(E1, E2, E3)
        beta = (strongest - E123) / (E0 - strongest)
        return beta

    print(f"b12 = {beta(E1, E2, E12):.2f}")
    print(f"b13 = {beta(E1, E3, E13):.2f}")
    print(f"b13 = {beta(E2, E3, E23):.2f}")
    b123 = beta3(E1, E2, E3, E123)
    bpairs = beta(E12, E13, E23)
    print(f"b123 = {beta3(E1, E2, E3, E123):.2f}")
    print(f"bpairs = {beta(E12, E13, E23):.2f}")
    print(f"emergent = {b123+bpairs:.2f}")

    fig = plt.figure()
    ax = fig.add_subplot(111)

    d1_name = stack.drug_names[0][0]
    d2_name = stack.drug_names[1][0]
    d3_name = stack.drug_names[2][0]

    labels = [
        stack.drug_names[0][:4],
        stack.drug_names[1][:4],
        stack.drug_names[2][:4],
        d1_name + "+" + d2_name,
        d1_name + "+" + d3_name,
        d2_name + "+" + d3_name,
        d1_name + "+" + d2_name + "+" + d3_name
    ]
    #labels = ["E1", "E2", "E3", "E12", "E13", "E23", "E123"]
    data = [E1, E2, E3, E12, E13, E23, E123]
    colormap = ["tab:red", "tab:green", "tab:blue", "y", "m", "c"]
    ax.bar(labels, data, color = colormap)
    ax.set_ylabel("Max Effect (Emax)")
    ax.axhline(min(E1,E2,E3), color = "0.2", linestyle = "dotted", label = "Single Drug Max")
    ax.axhline(min(E1,E2,E3,E12,E13,E23), color = "0.2", linestyle = "dashed", label = "Overall Max")
    ax.legend()
    plt.xticks(rotation = 60)
    plt.show()

def beta_plot(stack):
        # set regions
    d12_point = (0.5, 0.5, 0)
    d13_point = (0.5, 0, 0.5)
    d23_point = (0, 0.5, 0.5)

    region12 = leukdev.Region.from_point(d12_point, stack)
    region13 = leukdev.Region.from_point(d13_point, stack)
    region23 = leukdev.Region.from_point(d23_point, stack)


    d123_point = (1/3, 1/3, 1/3)
    region123 = leukdev.Region.from_point(d123_point, stack)
    
    
    # beta calculations
    E0 = 1
    E1 = stack.singledrugfits[0].Emax
    E2 = stack.singledrugfits[1].Emax
    E3 = stack.singledrugfits[2].Emax
    E12 = region12.fit_hill().Emax
    E13 = region13.fit_hill().Emax
    E23 = region23.fit_hill().Emax
    E123 = region123.fit_hill().Emax

    def beta(E1, E2, E3):
        E0 = 1
        beta = (min(E1, E2) - E3) / (E0 - min(E1, E2))
        return beta

    def beta3(E1, E2, E3, E123):
        E0 = 1
        strongest = min(E1, E2, E3)
        beta = (strongest - E123) / (E0 - strongest)
        return beta

    # beta plot
    fig = plt.figure()
    ax = fig.add_subplot(111)

    d1_name = stack.drug_names[0][0]
    d2_name = stack.drug_names[1][0]
    d3_name = stack.drug_names[2][0]

    labels = [
        d1_name + "+" + d2_name,
        d1_name + "+" + d3_name,
        d2_name + "+" + d3_name,
        d1_name + "+" + d2_name + "+" + d3_name
    ]
    print(labels)
    #labels = ["E1", "E2", "E3", "E12", "E13", "E23", "E123"]
    data = np.array([
        beta(E1, E2, E12),
        beta(E1, E3, E13),
        beta(E2, E3, E23),
        beta3(E1, E2, E3, E123)
    ])
    #data = data * -1
    print(data)
    #data = [E1, E2, E3, E12, E13, E23, E123]
    colormap = [_color("combo12"), _color("combo13"), _color("combo23"), _color("combo123")]
    ax.bar(labels, data, color = colormap)
    ax.set_ylabel("Synergistic Efficacy (beta)")
    plt.xticks(rotation = 60)
    ax.axhline(0, color = "0.2", label = "Single Drug Max")
    plt.show()

def FIC_plot(stack):
        # set regions
    d12_point = (0.5, 0.5, 0)
    d13_point = (0.5, 0, 0.5)
    d23_point = (0, 0.5, 0.5)

    region12 = leukdev.Region.from_point(d12_point, stack)
    region13 = leukdev.Region.from_point(d13_point, stack)
    region23 = leukdev.Region.from_point(d23_point, stack)


    d123_point = (1/3, 1/3, 1/3)
    region123 = leukdev.Region.from_point(d123_point, stack)
    
    
    C12 = region12.fit_hill().C
    C13 = region13.fit_hill().C
    C23 = region23.fit_hill().C
    C123 = region123.fit_hill().C
    A = C123
    B = (C12 * C13 * C23) ** (1/3)
    FIC3 = A / B

    # FIC plot
    fig = plt.figure()
    ax = fig.add_subplot(111)

    d1_name = stack.drug_names[0][0]
    d2_name = stack.drug_names[1][0]
    d3_name = stack.drug_names[2][0]

    labels = [
        d1_name + "+" + d2_name,
        d1_name + "+" + d3_name,
        d2_name + "+" + d3_name,
        d1_name + "+" + d2_name + "+" + d3_name,
        "pairwise estimate",
        "emergent 3-way"
    ]
    #labels = ["E1", "E2", "E3", "E12", "E13", "E23", "E123"]
    data = np.array([
        C12,
        C13,
        C23,
        C123,
        B,
        FIC3
    ])
    data = np.log2(data)
    #data = [E1, E2, E3, E12, E13, E23, E123]
    colors = ["y", "m", "c", "0.2", "0.4", "0.4"]
    ax.barh(labels, data, color = colors)
    ax.set_xlabel("Synergistic Potency (log2[FIC])")
    ax.set_xlim(-1, 1.5)
    plt.xticks(rotation = 60)
    ax.axvline(0, color = "0.2", label = "Single Drug Max")
    plt.show()

def ratiometric_pair_plot(stack, drug1_idx, drug2_idx):
    # Ratiometric EC50 and Emax plots
    # start between drugs 1 and 2

    #drug1_idx = 2
    #drug2_idx = 3
    colors = ["#ff1f5b", "#00cd6c", "#009ade", "#af58ba", "#ffc61e", "#f28522", "#a0b1ba", "#a6761d"]
    EC50_1 = stack.singledrugfits[0].C
    EC50_2 = stack.singledrugfits[1].C
    EC50_3 = stack.singledrugfits[2].C

    ratio1 = np.arange(0.25, 0.8, 0.05)
    print(ratio1)

    EC50 = []
    Emax = []
    for r1 in ratio1:
        point = [0, 0, 0]
        point[drug1_idx-1] = r1
        point[drug2_idx-1] = 1-r1
        point = tuple(point)
        bounds = leukdev.Zone.generate_bounds(point)
        
        region = leukdev.Region.from_bounds(bounds, stack)
        y = region.get_viability()
        y= y/y[-1] #normalize to control device viability
        d1 = stack.dose[:,drug1_idx-1]
        d2 = stack.dose[:,drug2_idx-1]
        d1norm = d1 * point[0]
        d1norm = d1norm / EC50_1
        d2norm = d2 * point[1]
        d2norm = d2norm / EC50_2
        x = d1norm + d2norm

        hillfit = Hill.fit(x,y)
        EC50.append(hillfit.C)
        Emax.append(hillfit.Emax)
        
    fig = plt.figure(figsize = (4,4), dpi=150)
    ax = fig.add_subplot(111)
    newx = np.array(ratio1)
    newx = np.insert(newx, 0, 0)
    newx = np.append(newx, 1)
    newy = np.array(EC50)
    newy = np.insert(newy, 0, 1)
    newy = np.append(newy, 1)
    ax.set_xlabel("Proportion of Drug A to B")
    ax.set_ylabel("Normalized EC50")
    ax.plot(newx, newy, 'o-', c = "tab:red", label = "Combo EC50")
    ticks = [0, 0.25, 1/3, 0.5, 2/3, 0.75, 1]
    tick_labels = ["0% A\n100% B", "1:3", "1:2", "1:1", "2:1", "3:1", "100% A\n0% B"]
    ax.set_xticks(ticks = ticks, labels = tick_labels)

    ax.set_xlim(0.25, 0.75)

    ax2 = ax.twinx()
    ax2.plot(ratio1, Emax, 'o-', c= "tab:blue", label = "Combo Emax")
    ax2.set_ylabel("Efficacy (Emax)")
    ax2.set_ylim(0, 1)
    #fig.legend(loc = "center")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1+h2, l1+l2, loc=0)
    plt.show()
