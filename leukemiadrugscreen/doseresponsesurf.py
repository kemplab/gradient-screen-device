import sympy as sym
import numpy as np
from matplotlib import pyplot as plt

def greco3D(C1, C2, C3, **kwargs):
   
    if len(np.shape(C1)) > 1:
        C1 = C1.flatten()
        C2 = C2.flatten()
        C3 = C3.flatten()

   # drug parameters
    E0 = kwargs.get("E0", 1) 
    Emax = kwargs.get("Emax", 0)

    # parameters for m
    alpha_m1 = kwargs.get("alpha_m1", 1)
    alpha_m2 = kwargs.get("alpha_m2", 1)
    alpha_m3 = kwargs.get("alpha_m3", 1)
    beta_m12 = kwargs.get("beta_m12", 0)
    beta_m13 = kwargs.get("beta_m13", 0)
    beta_m23 = kwargs.get("beta_m23", 0)
    gamma_m12 = kwargs.get("gamma_m12", 0)
    gamma_m13 = kwargs.get("gamma_m13", 0)
    gamma_m23 = kwargs.get("gamma_m23", 0)
    delta_m123 = kwargs.get("delta_m123", 0)

    # parameters for CI50
    alpha_D1 = kwargs.get("alpha_D1", 0)
    alpha_D2 = kwargs.get("alpha_D2", 0)
    alpha_D3 = kwargs.get("alpha_D3", 0)
    beta_D12 = kwargs.get("beta_D12", 0)
    beta_D13 = kwargs.get("beta_D13", 0)
    beta_D23 = kwargs.get("beta_D23", 0)
    gamma_D12 = kwargs.get("gamma_D12", 0)
    gamma_D13 = kwargs.get("gamma_D13", 0)
    gamma_D23 = kwargs.get("gamma_D23", 0)
    delta_D123 = kwargs.get("delta_D123", 0)


    IC50_1 = kwargs.get("IC50_1", 1)
    IC50_2 = kwargs.get("IC50_2", 1)
    IC50_3 = kwargs.get("IC50_3", 1)

    IC50 = [IC50_1, IC50_2, IC50_3] # for drugs A, B, and C

    T = C1/IC50[0] + C2/IC50[1] + C3/IC50[2]
    XA = (C1/IC50[0])/T
    XB = (C2/IC50[1])/T
    XC = (C3/IC50[2])/T

    m = (alpha_m1 * XA) + (alpha_m2 * XB) + \
        (alpha_m3 * XC) + (beta_m12 * XA * XB) + \
        (beta_m13 * XA * XC) + (beta_m23 * XB * XC) + \
        (gamma_m12 * XA * XB * (XA - XB)) + \
        (gamma_m13 * XA * XC * (XA - XC)) + \
        (gamma_m23 * XB * XC * (XB - XC)) + \
        (delta_m123 * XA * XB * XC)

    CI50 = 10 ** ((1 - XA) * (1 - XB) * (1 - XC) * ((alpha_D1 * XA) + (alpha_D2 * XB) + \
        (alpha_D3 * XC) + (beta_D12 * XA * XB) + (beta_D13 * XA * XC) + (beta_D23 * XB * XC) + \
        (gamma_D12 * XA * XB * (XA - XB)) + (gamma_D13 * XA * XC * (XA - XC)) + \
        (gamma_D23 * XB * XC * (XB - XC)) + (delta_D123 * XA * XB * XC)))

    A = (Emax - E0)
    B = (T / (CI50)) ** (m)
    E = ((A * B) / (1 + B)) + E0

    return E

def E_test():
    A = np.linspace(0, 1, 11)
    B = np.linspace(1, 0, 11)
    C = np.linspace(0, 0, 11)

    E = greco3D(A, B, C, beta_D12 = -4)
    
    fig, ax = plt.subplots()
    ax.plot(A, E)

    ax.set_xlabel("Drug A")
    ax.set_ylabel("Drug effect")
    ax.set_ylim(0, 1)
    plt.show()


if __name__ == "__main__":
    E_test()


