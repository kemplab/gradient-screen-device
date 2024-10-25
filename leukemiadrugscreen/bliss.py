import leukemiadrugscreen as leukdev
from analysis import *
import pandas as pd

stack = load_folder("rep3")
dev1 = stack.devices[0]
#filtered = dev1.select_zone(c1 = 0.5, c3 = 0.5, method = "nearest", n = 100, plot = False)

# Bliss analysis
print(stack.get_num_devices())

# create an empty data frame
df12 = pd.DataFrame(columns = ["drug1", "drug2", "d1", "d2", "g1", "g2", "g12_pred", "g12"])
df13 = pd.DataFrame(columns = ["drug1", "drug3", "d1", "d3", "g1", "g3", "g13_pred", "g13"])
df23 = pd.DataFrame(columns = ["drug2", "drug3", "d2", "d3", "g2", "g3", "g23_pred", "g23"])

print(stack.dose)
d1name = stack.drug_names[0]
d2name = stack.drug_names[1]
d3name = stack.drug_names[2]

for i_dev in range(1, stack.get_num_devices()-1):
    #i_dev = 4
    d1 = stack.dose[i_dev, 0]
    d2 = stack.dose[i_dev, 1]
    d3 = stack.dose[i_dev, 2]

    g1 = stack.devices[i_dev].get_zone(1)[0]
    g2 = stack.devices[i_dev].get_zone(3)[0]
    g3 = stack.devices[i_dev].get_zone(5)[0]

    g12 = stack.devices[i_dev-1].get_zone(7)[0]
    g13 = stack.devices[i_dev-1].get_zone(8)[0]
    g23 = stack.devices[i_dev-1].get_zone(9)[0]

    g12_pred = g1 * g2
    g13_pred = g1 * g3
    g23_pred = g2 * g3




    df12 = df12.append({"drug1": d1name, "drug2":d2name, "d1": d1, "d2": d2, "g1": g1, "g2": g2, "g12_pred": g12_pred, "g12": g12}, ignore_index = True)
    df13 = df13.append({"drug1": d1name, "drug3":d3name, "d1": d1, "d3": d2, "g1": g1, "g3": g3, "g13_pred": g13_pred, "g13": g13}, ignore_index = True)
    df23 = df23.append({"drug2": d2name, "drug3":d3name, "d2": d2, "d3": d3, "g2": g2, "g3": g3, "g23_pred": g23_pred, "g23": g23}, ignore_index = True)


# how to do appending code

print(df12)
print(df13)
print(df23)