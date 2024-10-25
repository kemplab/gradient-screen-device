import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import math
import statistics as stats
import synergy
import os
import re
from .drugfit import Hill
import configparser
from . import plot
from . import ternary_tools
from . import devicesim



class Zone():
    """A class used to represent a ternary zone for a device.

    ...
    Attributes
    ----------
    bounds : tuple
        tuple of lower and upper bounds lists for ternary coordinates ([a0,b0,c0],[a1,b1,c1])
    cells : pandas DataFrame
        table of cell locations, concentrations, and live/dead state
    device : Device
        the associated Device object for the Zone
    viability : float
        viability (live/total) of the cells in the Zone
    num_cells : int
        total number of detected cells in the zone

    Methods
    -------
    generate_bounds(tern_point, radius)
        takes a ternary point and turns it into a set of hexagonal bounds with given radius
    single_drug_bounds(drug1_idx, drug2_idx)
        ???
    from_bounds(zone_bounds, device)
        returns a Zone from Device based on bounds
    from_point(tern_point, device, radius)
        returns a Zone from Device centered at tern_point with radius
    plot()
        plots the Zone superimposed on Device in spatial and ternary coordinates
    get_point()
        estimate the average concentrations in zone, output as ternary point
    """

    def __init__(self, bounds, selected_cells, device):
        self.bounds = bounds
        self.cells = selected_cells
        self.device = device
        self.viability = selected_cells["live"].mean()
        self.num_cells = len(selected_cells)

    @classmethod
    def generate_bounds(cls, tern_point, radius = 0.05):
        a_min = tern_point[0] - radius
        a_max = tern_point[0] + radius
        b_min = tern_point[1] - radius
        b_max = tern_point[1] + radius
        c_min = tern_point[2] - radius
        c_max = tern_point[2] + radius
        zone_bounds = ([a_min, b_min, c_min], [a_max, b_max, c_max])
        return zone_bounds
    
    @classmethod
    def single_drug_bounds(cls, drug1_idx, drug2_idx):
        if drug1_idx == drug2_idx:
            raise ValueError("Same value was used for indices for drug 1 and 2")
        
        threshold_min = 0.7
        threshold_drug2 = 0.4
        threshold_drug3 = 0.1
        zone_bounds = ([0, 0, 0], [1, 1, 1])
        zone_bounds[0][drug1_idx-1] = threshold_min
        zone_bounds[1][drug2_idx-1] = threshold_drug2

        drug3_idx = [1, 2, 3]
        drug3_idx.remove(drug1_idx)
        drug3_idx.remove(drug2_idx)
        drug3_idx = drug3_idx[0]
        zone_bounds[1][drug3_idx-1] = threshold_drug3
        
        return zone_bounds

    @classmethod
    def from_bounds(cls, zone_bounds, device):
        # format: ([a_min, b_min, c_min], [a_max, b_max, c_max])
        a_min = zone_bounds[0][0]
        b_min = zone_bounds[0][1]
        c_min = zone_bounds[0][2]
        a_max = zone_bounds[1][0]
        b_max = zone_bounds[1][1]
        c_max = zone_bounds[1][2]

        selected = device.dataTable[
            (device.dataTable["c1"] > a_min) & (device.dataTable["c1"] < a_max) &
            (device.dataTable["c2"] > b_min) & (device.dataTable["c2"] < b_max) &
            (device.dataTable["c3"] > c_min) & (device.dataTable["c3"] < c_max)
            ]
        return cls(zone_bounds, selected, device)

    @classmethod
    def from_point(cls, tern_point, device, radius = 0.05):
        zone_bounds = cls.generate_bounds(tern_point, radius)
        return cls.from_bounds(zone_bounds, device)

    def plot(self):
        fig, ax = plt.subplots(1,2)
        ax[0].set_aspect('equal')
        ax[1].set_aspect('equal')

        ax[0].scatter(self.device.dataTable["x"], self.device.dataTable["y"], c = '0.8')
        ax[0].scatter(self.cells["x"], self.cells["y"])

        ternx, terny = ternary_tools.terncoords(self.device.dataTable[["c1", "c2", "c3"]])
        ternx_sel, terny_sel = ternary_tools.terncoords(self.cells[["c1", "c2", "c3"]])
        ax[1].scatter(ternx, terny, c = '0.8')
        ax[1].scatter(ternx_sel, terny_sel)

        plt.show()

    def get_point(self):
        '''Estimate the average concentrations in zone, output as ternary point.'''
        a = self.cells["c1"].mean()
        b = self.cells["c2"].mean()
        c = self.cells["c3"].mean()
        point = (a, b, c)
        return point

class Device():
    """A class used to represent the data for a Device.

    ...
    Attributes
    ----------
    dataTable : pandas DataFrame
    drug_names : list
        list of the three drug names
    input_dose : list
        
    dose_units : str

    Methods
    -------
    generate_bounds(tern_point, radius)
        takes a ternary point and turns it into a set of hexagonal bounds with given radius
    single_drug_bounds(drug1_idx, drug2_idx)
        ???
    from_bounds(zone_bounds, device)
        returns a Zone from Device based on bounds
    from_point(tern_point, device, radius)
        returns a Zone from Device centered at tern_point with radius
    plot()
        plots the Zone superimposed on Device in spatial and ternary coordinates
    get_point()
        estimate the average concentrations in zone, output as ternary point
    """
    
    def __init__(self, dataTable, drug_names, input_dose, dose_units):
        '''Create a drug Device object with a given dataTable of cell information,
        and maximum starting drug concentrations.
        DataTable just needs x,y,c1,c2,c3,live'''
        self.dataTable = dataTable # pandas dataframe
        self.drug_names = drug_names
        self.input_dose = input_dose
        self.dose_units = dose_units

    @classmethod
    def from_excel(cls, filename, drug_names, input_doses, dose_units):
        df = pd.read_excel(filename, header=None)
        df.rename({0: "x", 1: "y", 2: "c1", 3: "c2", 4: "c3", 5: "live"}, axis=1, inplace=True)
        return cls(df, drug_names, input_doses, dose_units)

    def print_summary(self):
        print(f"""Device: {self.drug_names[0][0:4]}[{self.input_dose[0]} 
            {self.dose_units[0]}]-{self.drug_names[1][0:4]}[{self.input_dose[1]} 
            {self.dose_units[1]}]-{self.drug_names[2][0:4]}[{self.input_dose[2]} 
            {self.dose_units[2]}]""")


    def get_zone_ratio(self, target_ratio, drug1_idx, drug2_idx):
        '''ratio is the x:1 ratio of drug c1:c2
        Uses a linear search algorithm to find the zone with target ratio
        Drug 1 and drug 2 are just the indexes of which drug 1, 2, or 3
        '''
        df = self.dataTable.copy()
        zone_threshold = 0.68
        combo_threshold = 0.45
        drug3_threshold = 0.27
        window_size = 0.08

        max_attempts = 100
        attempts = 0

        c1_start = 0.3
        error = 100 # large value guaranteed to trigger
        #target_ratio = 1.4 # starting
        error_thresh = 0.02 # if within 1%

        while abs(error) > error_thresh:
        
            if attempts >= max_attempts:
                print("Warning: max attempts for ratio adjustment reached. Synergy calculations may be off.")
                break
            
            if (drug1_idx == 1) & (drug2_idx == 2):
                df_filt = df[(df["c1"] > c1_start) & (df["c2"] > (1-c1_start)-window_size) & (df["c3"] < 0.10)]
            elif (drug1_idx == 1) & (drug2_idx == 3):
                df_filt = df[(df["c1"] > c1_start) & (df["c3"] > (1-c1_start)-window_size) & (df["c2"] < 0.10)]
            elif (drug1_idx == 2) & (drug2_idx == 3):
                df_filt = df[(df["c2"] > c1_start) & (df["c3"] > (1-c1_start)-window_size) & (df["c1"] < 0.10)]
            else:
                print("Some error in assigning drug indexes")
            
            #print(f"{len(df_filt)} cells counted in region")
            viability_max = df["live"].mean()
            c1_avg = df_filt["c1"].mean()
            c2_avg = df_filt["c2"].mean()
            c3_avg = df_filt["c3"].mean()
            c_avg = [c1_avg, c2_avg, c3_avg] # make a list

            ratio_12 = c_avg[drug1_idx-1] / c_avg[drug2_idx-1]
            error = (ratio_12 - target_ratio)/target_ratio
            #print(f"Error: {error:.2f}")

            c1_start -= np.sign(error) * 0.01
            #print(f"\nNew starting value = {c1_start:.2f}")

            attempts += 1 # iterate attempts

        #print(f"Viability of section: {viability_max:.2f}")
        #print(f"Drug ratio C1/C2: {ratio_12:.2f}")
        return df_filt["live"].mean(), c_avg[drug1_idx-1], c_avg[drug2_idx-1]


class Region():
    def __init__(self, zones, stack):
        self.zones = zones
        self.stack = stack
        self.bounds = self.zones[0].bounds

    @classmethod
    def from_bounds(cls, zone_bounds, stack):
        zones = []
        for i in range(stack.get_num_devices()):
            zones.append(Zone.from_bounds(zone_bounds, stack.devices[i]))
        return cls(zones, stack)
    
    @classmethod
    def from_point(cls, point, stack):
        zones = []
        for i in range(stack.get_num_devices()):
            zones.append(Zone.from_point(point, stack.devices[i]))
        region = cls(zones, stack)
        region.point = point
        return region

    def get_viability(self):
        viability = []
        for zone in self.zones:
            viability.append(zone.viability)
        return viability

    def get_GR(self):
        '''Calculates growth rate inhibition from a Region.'''

        # need to have cell counts of equivalent region to do comparison
        # make a fake device and region based on original
        conc_map = devicesim.load_map()
        conc_map["live"] = np.ones((len(conc_map), 1), dtype = int)
        ctrl_device = Device(conc_map, ["A", "B", "C"], [0, 0, 0], "a.u.")
        ctrl_zone = Zone.from_bounds(self.bounds, ctrl_device)
        proportion = ctrl_zone.num_cells / len(conc_map)

        # initialize variables
        GR = []
        #num_cells = [] # not currently used
        #cell_concentration = 2.75e6 # need to change this to get from the .ini
        #device_volume = 4 # uL

        # normalize number of starting cells to region proportion from simulated device        
        #x_0_raw = cell_concentration / 1000 * device_volume # starting cells
        #x_0 = x_0_raw * proportion

        # identify control device and get cell counts
        idx_ctrl = self.stack.get_idx_ctrl()
        x_ctrl = sum(self.zones[idx_ctrl].cells["live"]) # get live cells in control device
        
        # calculate growth rate of control
        # error occurs here when x_0 is 0
        # which only occurs when proportion is 0
        # which occurs when the ctrl_zone has a region with no cells?

        if x_ctrl == 0:
            return np.nan


        #k_0 = np.log2(x_ctrl / x_0)

        # calculate growth rate inhibition across all devices relative to control growth rate
        for i in range(len(self.zones)):
            x_c = sum(self.zones[i].cells["live"])
            #k_c = np.log2(x_c / x_0)
            #GR.append(2**(k_c / k_0) - 1)
            GR.append(x_c / x_ctrl)
            #num_cells.append(x_c)
        return GR

        

    def get_point(self):
        if hasattr(self, "point"):
            return self.point
        else:
            return self.zones[0].get_point()

    def fit_hill(self):
        y = self.get_viability()
        y=y/y[-1]  #normalize to control device viability
        point = self.get_point()
        EC50_1 = self.stack.singledrugfits[0].C        
        EC50_2 = self.stack.singledrugfits[1].C
        EC50_3 = self.stack.singledrugfits[2].C
        x1 = self.stack.dose[:, 0] * point[0] / EC50_1
        x2 = self.stack.dose[:, 1] * point[1] / EC50_2
        x3 = self.stack.dose[:, 2] * point[2] / EC50_3
        x = x1 + x2 + x3
        hillfit = Hill.fit(x, y)
        return hillfit


class DeviceStack():
    '''A collection of devices with the same drugs. Need an .ini or .config within the folder
    '''
    
    def __init__(self, devices_list):
        self.devices = devices_list
        self.drug_names = self.devices[0].drug_names
        self.dose = np.empty((len(devices_list), 3))
        for i in range(len(self.devices)):
            self.dose[i,:] = self.devices[i].input_dose
        self.dose_units = self.devices[0].dose_units
        # fix this later to check consistency between ALL devices
        self.set_singledrugfits()

    @classmethod
    def from_folder(cls, folder_name):
        working_dir = os.getcwd()

        folder_path = os.path.join(working_dir, folder_name) # assuming it is a subfolder of working directory
        
        # scan folder for initial list
        file_list = []
        for file in sorted(os.listdir(folder_path)):
            if file.endswith(".xlsx") or file.endswith(".csv"):
                file_list.append(file)
            elif file.endswith(".ini") or file.endswith(".config") or file.endswith(".cfg"):
                ini_path = os.path.join(folder_path, file)

        # drug info text file
        if "ini_path" in locals() == False:
            raise Exception("No .ini file found in data folder.")

        config = configparser.ConfigParser()
        config.read(ini_path)
        
        # check number of files matches .ini
        num_devices = len(file_list)
        if config.has_section("Files"):
            if len(file_list) != len(config.options("Files")):
                raise Exception("Number of data files does not match list given in .ini files.")

        # make sure values are unique by converting to set
        if len(config.options("Files")) != len(set(config.options("Files"))):
            raise Exception("Not all files in .ini list are unique")

        drug_names = [config["Drug1"]["name"], config["Drug2"]["name"], config["Drug3"]["name"]]
        dose_units = [config["Drug1"]["units"], config["Drug2"]["units"], config["Drug3"]["units"]]
        abbreviations = [config["Drug1"]["abbreviation"], config["Drug2"]["abbreviation"], config["Drug3"]["abbreviation"]]
        dose = np.empty((num_devices, 3))
        dose[:,0] = np.array(eval(config["Drug1"]["doses"]))
        dose[:,1] = np.array(eval(config["Drug2"]["doses"]))
        dose[:,2] = np.array(eval(config["Drug3"]["doses"]))

        # new list based on .ini file list (circumvents alphabetical misordering)
        file_list_ini = []
        file_num_ini = []
        for option in config.options("Files"):
            file_num_ini.append(re.findall(r'\d+', file))
            file_list_ini.append(option)

        devices = []
        i = 0
        for file in file_list:
            devices.append(Device.from_excel(os.path.join(folder_path, file), drug_names, dose[i,:], dose_units))
            i += 1

        output = cls(devices)
        output.abbrevs = abbreviations

        return output

    def __str__(self):
        return "DeviceStack:{}-{}-{}".format\
            (self.drug_names[0][0:4], self.drug_names[1][0:4], self.drug_names[2][0:4])

    def get_num_devices(self):
        '''get the number of devices in stack'''
        return len(self.devices)

    def bliss(drug1, drug2):
        '''get bliss synergy between any 2 drugs'''

    def set_singledrugfits(self):
        '''Hill fit for all 3 single drug regions'''
        d1_bounds = ([0.68, 0, 0], [1, 0.4, 0.4])
        region1 = Region.from_bounds(d1_bounds, self)
        d2_bounds = ([0, 0.68, 0], [0.4, 1, 0.4])
        region2 = Region.from_bounds(d2_bounds, self)
        d3_bounds = ([0, 0, 0.68], [0.4, 0.4, 1])
        region3 = Region.from_bounds(d3_bounds, self)

        def singledrugfit(region, drug_idx):
            y = region.get_viability()
            y=y/y[-1] #normalize to control device viability
            x = region.stack.dose[:, drug_idx-1]
            hillfit = Hill.fit(x, y)
            return(hillfit)

        def singledrugfitGR(region, drug_idx):
            y = region.get_GR()
            x = region.stack.dose[:, drug_idx-1]
            hillfit = Hill.fit(x, y)
            return(hillfit)

        self.singledrugfits = [
            singledrugfit(region1, 1),
            singledrugfit(region2, 2),
            singledrugfit(region3, 3)
        ]

        self.singledrugfitsGR = [
            singledrugfitGR(region1, 1),
            singledrugfitGR(region2, 2),
            singledrugfitGR(region3, 3)
        ]

        self.singledrugregions = [region1, region2, region3]

    def normalize_dose():
        pass

    def get_region(self, zone_bounds):
        '''Get single drug fit of drug, at specific condition
        Side 1, Side 2, or Average of the two'''

        zones = []
        for i in range(self.get_num_devices()):
            zones.append(self.devices[i].get_zone(zone_bounds))

        region = Region(zones, self)
        return region

    def get_idx_ctrl(self):
        '''Returns the device index for the control device (0 concentration doses)'''
        idx_ctrl = self.dose == 0
        for row in range(idx_ctrl.shape[0]):
            true_sum = sum(idx_ctrl[row, :])
            if true_sum == 3: # if there are 3 zeros
                idx_ctrl = row
            
        if type(idx_ctrl) != int:
            raise Exception("Could not find control device where all doses equal 0.")
        return idx_ctrl

    def diamond(self, drug1_idx, drug2_idx, plot = False):
        '''Apply diamond methology to synergy scoring between 2 drugs'''
        # get EC50 of drug 1 & 2
        EC50_1, EC50_2, EC50_agg, hill_agg_1 = self.single_drug_fit(drug1_idx)
        drug1_EC50 = stats.mean([EC50_1, EC50_2, EC50_agg]) # changed from median to mean
        EC50_1, EC50_2, EC50_agg, hill_agg_2 = self.single_drug_fit(drug2_idx)
        drug2_EC50 = stats.mean([EC50_1, EC50_2, EC50_agg])

        drug1_name = self.drug_names[drug1_idx-1]
        doses1 = self.dose[:,drug1_idx-1]

        drug2_name = self.drug_names[drug2_idx-1]
        doses2 = self.dose[:,drug2_idx-1]

        target1, target2 = drug1_EC50/2, drug2_EC50/2
        #print(f"Target EC50 for drug 1: {target1:.2f}")
        #print(f"Target EC50 for drug 2: {target2:.2f}")

        # I now have to find a device in which the targets are matched up and can use a ratio
        # Translate the actual device concentrations into multiples of EC50?
        #print(f"Normalized EC50 drug1: {np.around(doses1/drug1_EC50, 2)}")
        #print(f"Normalized EC50 drug2: {np.around(doses2/drug2_EC50, 2)}")

        dose1_norm = doses1 / drug1_EC50
        dose2_norm = doses2 / drug2_EC50

        ratio = 1/np.mean(dose1_norm[np.nonzero(dose1_norm)] / dose2_norm[np.nonzero(dose2_norm)])
        print(f"Target ratio is: {ratio:.2f}")

        EC50_combo, hill_combo = self.ratio_zone_fit(ratio, drug1_idx, drug2_idx, dose1_norm, dose2_norm)
        print(f"Analysis of {drug1_name} + {drug2_name} combination using diamond:")
        #print(f"EC50 of combo is: {EC50_combo:.2f}")

        # plotting of diamond
        if plot == True:
            fig, ax = plt.subplots()
            
            
            def plot_hill(hill_object, label):
                ax.scatter((hill_object.x / hill_object.get_parameters()[3]), hill_object.y)
                #x_fit = np.linspace(0, 7)
                x_fit = np.linspace(0, hill_object.x.max(), 100)
                y_fit = hill_object.E(x_fit)
                x_fit = x_fit / hill_object.get_parameters()[3]
                ax.plot(x_fit, y_fit, label = label)

            plot_hill(hill_agg_1, label = ("Drug " + str(drug1_idx)))
            plot_hill(hill_agg_2, label = ("Drug " + str(drug2_idx)))

            ax.scatter((hill_combo.x), hill_combo.y)
            x_fit = np.linspace(0, hill_combo.x.max(), 100)
            y_fit = hill_combo.E(x_fit)
            ax.plot(x_fit, y_fit, label = "Drug " + str(drug1_idx) + "+" + str(drug2_idx))


            ax.legend()
            ax.set_xlabel('Dose (Normalized EC50)')
            ax.set_ylabel('Cell Viability')
            #ax.set_xscale('log')
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            fig.set_size_inches(4, 3)
            fig.tight_layout()
            plt.show()

            



        return EC50_combo
    
    def ratio_zone_fit(self, ratio, drug1_idx, drug2_idx, dose1_norm, dose2_norm):
        '''Fits the zone indicated by the specific ratio.
        Also needs the normalized EC50 doses of drugs 1 and 2.
        Change this to calculate the normalized EC50 doses of drug 1 and 2.'''
        

        viability = np.empty(self.get_num_devices())
        c1 = np.empty(self.get_num_devices())
        c2 = np.empty(self.get_num_devices())
        c3 = np.empty(self.get_num_devices())
        plotting = False

        for i in range(self.get_num_devices()):
            viability[i], c1[i], c2[i] = self.devices[i].get_zone_ratio(ratio, drug1_idx, drug2_idx)
            #viability[i], c1[i], c2[i], c3[i] = self.devices[i].get_zone(7)
        rescaled_EC50 = ((np.multiply(c1, dose1_norm) + np.multiply(c2, dose2_norm))/2)*2
        # do I need to divide this by 2 again? I think so: or is it multiply? Let's try multiply
        #print(np.around(rescaled_EC50,2))

        hill_combo = Hill()
        hill_combo.fit(rescaled_EC50, viability)
        EC50_combo = hill_combo.get_parameters()[3]

        return EC50_combo, hill_combo

    def custom_ratio_fit(self, ratio, drug1_idx, drug2_idx, plotting = False):
        '''Fit 2-drug combo viabilities between a specified ratio of drug 1:2'''
        # get EC50 of drug 1 & 2
        EC50_1, EC50_2, EC50_agg = self.single_drug_fit(drug1_idx)
        drug1_EC50 = stats.median([EC50_1, EC50_2, EC50_agg])
        EC50_1, EC50_2, EC50_agg = self.single_drug_fit(drug2_idx)
        drug2_EC50 = stats.median([EC50_1, EC50_2, EC50_agg])

        drug1_name = self.drug_names[drug1_idx-1]
        doses1 = self.dose[:,drug1_idx-1]

        drug2_name = self.drug_names[drug2_idx-1]
        doses2 = self.dose[:,drug2_idx-1]

        target1, target2 = drug1_EC50/2, drug2_EC50/2

        dose1_norm = doses1 / drug1_EC50
        dose2_norm = doses2 / drug2_EC50

        # code from ratio_zone_fit
        viability = np.empty(self.get_num_devices())
        c1 = np.empty(self.get_num_devices())
        c2 = np.empty(self.get_num_devices())

        for i in range(self.get_num_devices()):
            viability[i], c1[i], c2[i] = self.devices[i].get_zone_ratio(ratio, drug1_idx, drug2_idx)
            #viability[i], c1[i], c2[i], c3[i] = self.devices[i].get_zone(7)

        rescaled_EC50 = ((np.multiply(c1, dose1_norm) + np.multiply(c2, dose2_norm))/2)*2

        hill_combo = synergy.single.Hill()
        hill_combo.fit(rescaled_EC50, viability)
        EC50_combo = hill_combo.get_parameters()[3]

        if plotting == True:
            fig, ax1 = plt.subplots()
            ax1.plot(dose1_norm, viability)
            ax2 = ax1.twiny()
            ax2.plot(dose2_norm, viability, alpha = 0)
            plt.show()

    def drug3_fit(self, a, b, c):
        '''3-drug fit at a specified ratio'''
        viability = np.empty(self.get_num_devices())
        for i in range(self.get_num_devices()):
            viability[i] = self.devices[i].select_zone(c1 = a, c2 = b, c3 = c, method = "nearest", n = 100, plot = True)

        print(viability)

    def plot_ternary(self, response = "viability"):
        '''Plots ternary plots of DeviceStack'''
        if response == "viability":
            plot.ternary_plots(self)
        elif response == "GR":
            plot.ternary_plots_GR(self)
        else:
            raise Exception("Invalid value for response. Acceptable values are 'viability' or 'GR'.")

    def plot_zones(self, response = "viability"):
        '''Plots the single drug, equipotent, and 3-drug regions of a DeviceStack'''
        if response == "viability":
            plot.zone_plots(self)
        elif response == "GR":
            plot.zone_plots_GR(self)
        else:
            raise Exception("Invalid value for response. Acceptable values are 'viability' or 'GR'.")

    def emergent_3way(self):
        d12_point = (0.5, 0.5, 0)
        d13_point = (0.5, 0, 0.5)
        d23_point = (0, 0.5, 0.5)

        region12 = Region.from_point(d12_point, self)
        region13 = Region.from_point(d13_point, self)
        region23 = Region.from_point(d23_point, self)


        d123_point = (1/3, 1/3, 1/3)
        region123 = Region.from_point(d123_point, self)
        
        C12 = region12.fit_hill().C
        C13 = region13.fit_hill().C
        C23 = region23.fit_hill().C
        C123 = region123.fit_hill().C

        A = C123
        B = (C12 * C13 * C23) ** (1/3)
        print(f"FIC3 = {A:.2f}")
        print(f"lFIC2 = {B:.2f}")
        FIC3 = A / B
        print(f"Emergent 3-way synergy = {FIC3:.2f}")
        return(FIC3)

    def plot_Emax(self):
        '''Plots Emax values for equipotent and single drug ratios'''
        plot.emax_plot(self)

    def plot_beta(self):
        '''Plots beta values for equipotent and single drug ratios'''
        plot.beta_plot(self)

    def plot_FIC(self):
        '''Plots FIC values for equipotent and single drug ratios'''
        plot.FIC_plot(self)

    def plot_ratiometric_pair(self, drug1_idx, drug2_idx):
        '''Plots ratiometric plot for a pair of drugs'''
        plot.ratiometric_pair_plot(self, drug1_idx, drug2_idx)

    def get_ternary_table(self, filename = None):
        '''Generates and exports a CSV file for ternary table outputs'''
        tern_grid = ternary_tools.gen_tern_grid(11)
    
        viability = np.empty([len(tern_grid), self.get_num_devices()])
        num_cells = np.empty_like(viability)
        GR = np.empty_like(viability)

        for i in range(len(tern_grid)):
            tern_point = (tern_grid["a"][i], tern_grid["b"][i], tern_grid["c"][i])
            for x in range(self.get_num_devices()):
                zone = Zone.from_point(tern_point, self.devices[x])
                viability[i, x] = zone.cells["live"].mean()
                num_cells[i, x] = len(zone.cells)
            region = Region.from_point(tern_point, self)
            GR[i, :] = region.get_GR()

        # fit and normalize to EC50
        EC50_A = self.singledrugfits[0].C
        EC50_B = self.singledrugfits[1].C
        EC50_C = self.singledrugfits[2].C

        GR50_A = self.singledrugfitsGR[0].C
        GR50_B = self.singledrugfitsGR[1].C
        GR50_C = self.singledrugfitsGR[2].C

        EC50 = np.empty(len(tern_grid))
        Emax = np.empty_like(EC50)

        GR50 = np.empty(len(tern_grid))
        GRmax = np.empty_like(GR50)

        dose_inputs_A = self.dose[:,0].transpose()
        dose_inputs_B = self.dose[:,1].transpose()
        dose_inputs_C = self.dose[:,2].transpose()

        # calculations based on viability
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

        # same for GR
        for i in range(len(EC50)):
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

        # filter out bad EC50 and GR50
        #nan_idx = np.isnan(Emax)
        #EC50[nan_idx] = 1

        #nan_idx_GR = np.isnan(GRmax)


        # set each edge to where nan is to single_drug_emax
        drug1_grid_idx = tern_grid[(tern_grid["a"]>=0.8) & (tern_grid["b"]<0.2)].index.to_list()
        drug2_grid_idx = tern_grid[(tern_grid["b"]>=0.8) & (tern_grid["a"]<0.2)].index.to_list()
        drug3_grid_idx = tern_grid[(tern_grid["c"]>=0.8) & (tern_grid["a"]<0.2)].index.to_list()
        
        #Emax[drug1_grid_idx] = stack.singledrugfits[0].Emax
        
        Emax[drug1_grid_idx] = self.singledrugfits[0].Emax
        Emax[drug2_grid_idx] = self.singledrugfits[1].Emax
        Emax[drug3_grid_idx] = self.singledrugfits[2].Emax


        # make table
        tern_grid["Norm EC50"] = EC50
        tern_grid["Log2[Norm EC50]"] = np.log2(EC50)
        tern_grid["Emax"] = Emax
        tern_grid["GR50"] = GR50
        tern_grid["Log2[GR50]"] = np.log2(GR50)
        tern_grid["GRmax"] = GRmax

        viability_col = []
        GR_col = []
        num_cells_col = []
        for n in range(self.get_num_devices()):
            viability_col.append(f"viability_dev_{n+1}")
            GR_col.append(f"GR_dev_{n+1}")
            num_cells_col.append(f"num_cells_dev_{n+1}")
        viability_df = pd.DataFrame(viability, columns = viability_col)
        GR_df = pd.DataFrame(GR, columns = GR_col)
        num_cells_df = pd.DataFrame(num_cells, columns = num_cells_col)

        output = pd.concat([tern_grid, viability_df, GR_df, num_cells_df], axis = 1)
        
        if filename != None:
            output.to_csv(filename)

        return(output)

    def plot_cell_density(self):
        pass


def beta(E1, E2, E3):
    E0 = 1
    beta = (min(E1, E2) - E3) / (E0 - min(E1, E2))
    return beta

def beta3(E1, E2, E3, E123):
    E0 = 1
    strongest = min(E1, E2, E3)
    beta = (strongest - E123) / (E0 - strongest)
    return beta
