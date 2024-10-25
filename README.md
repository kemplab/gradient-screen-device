# drug-screening-device
## About
Analytical tools for finding synergy in combinations of treatments in a 3-drug screening microfluidic for pediatric leukemia. :microscope::pill::cocktail:
## Authors
Dan Y. Zhang 
- Please direct all requests to melissa.kemp@bme.gatech.edu or evelyn.williams@choa.org

# System Requirements
There are no specific hardware requirements for this code beyond those of a standard computer. It has been tested on Windows 11, and all python dependencies are included in the `requirements.txt` file. 

# Installation instructions
Download the necessary files, ensure Jupyter notebook is installed, as well as all packages in the `requirements.txt` file. This will take ~10 minutes. The code was developed using Python v3.10 and has been tested using v3.12.

# Main Tools
## Importing experiment results
Experiment results should be contained within a single folder. The results from analysis of each individual device are contained within an `.xlsx` file. Experiment metadata is contained within an `.ini` or `.config` file: the name of the file does not matter, as long as there is only a single copy in the experiment folder. At minimum, the metadata file should contain the following data (see `leukemiadrugscreen/example_data/JURKAT DVP` for an example):
- Names of drugs
- Doses of drugs in each device, formatted as a list
- Filenames of device results corresponding to device number

To import experiments, use the `from_folder` method from the `DeviceStack` class:

        DeviceStack.from_folder("folder_name")

## Synergy analysis 
Analysis code can be found in `leukemiadrugscreen/leukemiadrugscreen.py`. 
- Example usage of developed tools is contained in `examples.ipynb` and can be run with the example experimental data in `leukemiadrugscreen/example_data/JURKAT DVP`.
- To run, open `examples.ipynb`, import desired experiment results, and execute all cells. Total run time is ~2-3 minutes.
- Example analysis includes single drug fits and synergism in potency and in efficacy across all zones or combination ratios. Hill fits and ternary plots are generated, and a summary table with per zone data including EC50, Emax, Normalized EC50, GR50, GRmax, viability, growth rate, and total cell number is saved as a .csv file.

## Generate simulated data
Simulation code can be found in `leukemiadrugsim.py`. This utilizes `leukemiadrugscreen.py` as a module. Currently, all parameters are stored inline, but the following features are coming soon:
- Will have ability to read config files (.conf) to initialize simulation paramters.
- Alternatively, also working on ability to use a spreadsheet of multiple conditions as inputs for generating a set of simulations.

# Miscellaneous things to do:
- Folder 'Image Processing' is a work-in-progress, collection of old code used to rotate images onto grid

# Appendix A: Object Hierarchy
## Classes
- Device
    - A device consists of a dataTable, name of drugs, and dose of those drugs. This is typically imported from an Excel containing the single-cell spatial locations, concentration at those locations, and live/dead status.
- Zone
    - A subset of cells from a particular device. This contains the subseted dataframe, as well as a reference to the parent device. It also stores the tuple used to create the zone, in the form `([a_min, b_min, c_min], [a_max, b_max, c_max])`.
    - Zones can also be created directly from points on the ternary plot, which is a tuple of the ternary coordinates `(a, b, c)` where concentrations a + b + c = 1.
- DeviceStack
    - An object that collects multiple devices as a list. It also contains a list of drug names as well as a list of drug concentrations, i index as devices, j index as drugs.
- Region
    - Main data is a list of zones, and the tuple of zone input used to define the zone.
- Hill
    - A more generic object for Hill function fit. Stores the parameters (E0, Emax, h, C) as well as the original data used to create the Hill fit.

## Methods
- Device
- Zone
- DeviceStack
- Region
- Hill
