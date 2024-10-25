import configparser
import numpy as np
import os
import leukemiadrugscreen as drugscreen

def info2ini():

    print("INFO: Successfully saved device information in .\info.ini")

def read_ini(filename):
    config = configparser.ConfigParser()
    config.read(filename)
    sections = config.sections()
    print(f'Sections: {sections}')

    drug1name = config[sections[2]]["doses"]
    print(drug1name)
    print(type(drug1name))

    dose_list = np.array(eval(config["Drug1"]["doses"]))
    print(dose_list)
    print(type(dose_list))

    print(config.has_option("Drug1", "abbreviation"))


    if config.has_section("Drug1"):
        print("section drug 1 exists")
    else:
        raise Exception("Section Drug1 does not exist in .ini file")
    
def make_ini():
    # lets create that config file for next time...
    cfgfile = open("template.ini",'w')
    config = configparser.ConfigParser()
        
    config.add_section("Experiment")
    config.set('Person','HasEyes', 'True')
    config.set('Person','Age', '50')
    
    
    config.write(cfgfile)
    cfgfile.close()


def load_example():
    '''Loads example DeviceStack'''
    # load example
    example_data_path = os.path.join(os.path.dirname(drugscreen.__file__), "example_data")
    stack = drugscreen.Stack.from_folder(os.path.join(example_data_path, "DVP"))
    return stack