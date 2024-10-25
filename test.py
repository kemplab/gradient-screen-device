import leukemiadrugscreen as drugscreen
import matplotlib.pyplot as plt
# loading example DeviceStack
#example = drugscreen.utilities.load_example()
#example.plot_zones()
#example.plot_emax()
#example.plot_beta()
#example.plot_FIC()
#example.plot_ratiometric_pair(2,3)
#example.get_ternary_table("test_ternary_output.csv")

from leukemiadrugscreen.drugfit import Hill
import numpy as np
x = np.linspace(0, 10, 11)
test = Hill(E0 = 1, Emax = 0, h = 2, C = 3)

print(test.get_dose_from_E(0.308))
print(test.get_AUC(20))

y = test.E(x)

print(x)
print(y)

plt.scatter(x,y)
plt.show()