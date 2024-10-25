from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt


class Hill():
    '''4-parameter Hill equation. Modeling this after the synergy.single.hill package,
    but making some modifications because of nonstandard error handling leading to
    unexpected behavior when faced with a "bad" fit.'''
    
    def __init__(self, E0, Emax, h, C):
        self.E0 = E0
        self.Emax = Emax
        self.h = h
        self.C = C

    def _model(d, E0, Emax, h, C):
        dh = np.power(d,h)
        return E0 + (Emax-E0)*dh/(np.power(C,h)+dh)

    def _fit(x, y):
        initial_guess = [0.9, 0.1, 2.1, max(x)/2]
        bounds = ((0, 0, 1.5, 0), (1, 1, 10, max(x)*2))
        popt, pcov = curve_fit(Hill._model, x, y, p0 = initial_guess, bounds = bounds)
        return popt

    def E(self, d):
        '''Generate data using model and doses'''
        y_new = Hill._model(d, self.E0, self.Emax, self.h, self.C)
        return y_new

    @classmethod
    def fit(cls, d, E):
        '''Provide hill fit parameters E0, Emax, h, C from data'''
        # check if nan
        if (np.isnan(np.sum(d)) == True or np.isnan(np.sum(E)) == True):
            print("Found NaN")
            print(d)
            print(E)
            d_nan = np.isnan(d)
            E_nan = np.isnan(E)  

        if (E[-1] - E[0]) > 0.1:
            try:
                params = cls._fit(d,E)
                E0 = params[0]
                Emax = params[1]
                h = params[2]
                C = params[3]
                hill_fit = cls(E0, Emax, h, C)
                hill_fit.x = d
                hill_fit.y = E
                return hill_fit
            except RuntimeError:
                 print("Error - curve_fit failed")
                 E0 = 1
                 Emax = 1
                 C = max(d)
                 hill_fit = cls(E0, Emax, 2.1, C)
                 hill_fit.x = d
                 hill_fit.y = E
                 return hill_fit
            
        else:
            E0 = 1
            Emax = 1
            C = max(d)
            hill_fit = cls(E0, Emax, 2.1, C)
            hill_fit.x = d
            hill_fit.y = E
            return hill_fit

    def get_parameters(self):
        return (self.E0, self.Emax, self.h, self.C)

    def print_parameters(self):
        print(f"E0={self.E0:.2f}, Emax={self.Emax:.2f}, \
            h = {self.h:.2f}, C = {self.C:.2f}")

    def plot_curve(self, raw_data = True):
        '''If raw data exists, use it'''
        if hasattr(self, "x") == True:
            plt.scatter(self.x, self.y)
            x_plot = np.linspace(0, max(self.x), 100)
        else:
            x_plot = np.linspace(0, self.C*5, 100)
        y_plot = self.E(x_plot)
        plt.plot(x_plot, y_plot)

    def get_AUC(self, EC = 50):
        '''Get the AUC of the Hill fit from E0 up to the concentration that defines the ECXX'''
        x_max = self.get_dose_from_E(EC/100)
        x = np.linspace(0, x_max)
        AUC = np.trapz(self.E(x), x)
        return AUC

    def get_dose_from_E(self, E):
        h = self.h
        E0 = self.E0
        C = self.C
        Emax = self.Emax
        dose = np.power((E-E0)*C**h/(Emax-E), (1/h))
        return dose


def main():
    x = np.array([1, 2, 3, 4, 5, 6])
    y = Hill._model(x, 1, 0, 2, 2)

    rng = np.random.default_rng()
    y_noise = 0.1 * rng.normal(size = x.size)
    ydata = y + y_noise

    print(ydata)
    test_model = Hill.fit(x, ydata)
    print(test_model.x)
    print(test_model.y)
    print(test_model.E(test_model.x))
    test_model.plot_curve()
    
    sim_model = Hill(1, 0.2, 2, 2)
    sim_model.plot_curve()
    
    
    plt.show()
