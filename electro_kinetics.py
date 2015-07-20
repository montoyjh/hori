import numpy as np

from ase.units import kB

from hori.thermo import (AdsorbateThermodynamics,
                         ProtonElectronThermodynamics,
                         GasThermodynamics)
from hori.common import hbond_dict
from hori.pathway import Pathway, calculate_rxn_deltaG

Nskipsteps = 1

class electro_kinetics:
    def __init__(self, rate_functions, path, at, pet, gt, voltage, T, solver):
        self.rate_functions = rate_functions
        self.path = path
        self.at = at
        self.pet = pet
        self.gt = gt
        self.voltage = voltage
        self.T = T
        self.solver = solver
        self.gas_names = gt.gas_dict.keys()
        self.gas_names.append('pe')
        self.ads_names = at.data.keys()
        self.NG = len(self.gas_names) # number of pas phase species (including the pe-pair)
        self.NS = len(self.ads_names) # number of surface species (including the free site)
        self.NR = len(path.steps) - 1 # number of reactions
        self._make_stoich_matrix()
        self._calculate_K_r()
        self._calculate_k_r()

    def recalculate(self):
        self._calculate_K_r()
        self._calculate_k_r()

    def solve(self, theta_0, k_r=None, K_r=None):
        if k_r == None:
            self._calculate_k_r()
            k_r = self.k_r
        if K_r == None:
            self._calculate_K_r()
            K_r = self.K_r
        self.solver.initialize(self.NG, self.NS,
                               self.fugacities, self.stoich_matrix,
                               k_r, K_r)
        theta_c, dtheta_dt_c = self.solver.solve(theta_0)
        rates_c = self.solver.get_rates(theta_c)
        return theta_c, rates_c

    def get_spec_names(self):
        spec_names = []
        for gas in self.gas_names:
            spec_names.append(gas)
        for ads in self.ads_names:
            spec_names.append(ads+'*')
        self.spec_names = spec_names
        return spec_names

    def _make_stoich_matrix(self):
        NG = self.NG
        NS = self.NS
        NR = self.NR
        stoich_matrix = np.zeros([NR, NG+NS])
        r = 0 # counter for the reaction number
        for step in self.path.steps[Nskipsteps:]:
            priorstate,nextstate,rxn = step
            for x in rxn:
                if x['adsorbed'] :
                    c = NG + self.ads_names.index(x['formula'])
                    stoich_matrix[r,c] += x['count'] 
                else:
                    if x['pe'] != 0.0:
                        c = self.gas_names.index('pe')
                        stoich_matrix[r,c] += x['count'] 
                    try:
                        c = self.gas_names.index(x['formula'])
                        stoich_matrix[r,c] += x['count']
                    except ValueError:
                        pass
            r += 1 
        self.stoich_matrix = stoich_matrix

    def set_fugacities(self, fugacity_dict):
        """ Set the fugacity of 'gas phase' species.        
        The fugacity is relative to the standard state defined for the gas thermo dynamics object.
        fugacity_dict is a dictionary
        """
        #fugacities = np.ones(self.NG)
        fugacities = np.zeros(self.NG)
        for k, v in fugacity_dict.items():
            i = self.gas_names.index(k)
            fugacities[i] = v
        self.fugacities = fugacities

    def _calculate_k_r(self):
        """ Calculate the forward rate constants.
        """
        NR = self.NR
        k_r = np.zeros(NR)
        for i in range(NR):
            rf = self.rate_functions[i]
            k_r[i] = rf(U=self.voltage, T=self.T)
        self.k_r = k_r

    def _calculate_K_r(self):
        """ Calculate the equilibrium constants.
        """
        NR = self.NR
        K_r = np.zeros(NR)
        G_r = np.zeros(NR)
        for i in range(NR):
            G = calculate_rxn_deltaG(self.path.steps[i+Nskipsteps][-1],
                                     BG=self.at.G,
                                     mu=self.gt.G,
                                     mu_pe=self.pet.G)
            G_r[i] = G
            K_r[i] = np.exp(- G / (kB * self.T))
        self.G_r = G_r
        self.K_r = K_r

    def set_temperature(self, T, recalc=True):
        """ Set the temperature of the simulation.
        """
        self.T = T
        self.at.set_temperature(self.T)
        self.gt.set_all_temperatures(self.T)
        self.pet.set_temperature(self.T)
        if recalc:
            self._calculate_K_r()
            
    def set_voltage(self, voltage, recalc=True):
        self.pet.set_voltage(voltage)
        self.voltage = voltage
        if recalc:
            self._calculate_K_r()
            self._calculate_k_r()

    def set_gas_standardstate(self, gas, fugacity, recalc=True):
        """ Change standard state for gas species.
        It is recommended to use set_fugacities to change the pressure or concentration within the simulation.
        """
        self.gt.set_gas_fugacity(gas, fugacity)
        if recalc:
            self._calculate_K_r()

    def set_all_standardstates(self, fugacity, recalc=True): 
        """ Change standard state for all gas species.
        It is recommended to use set_fugacities to change the pressure or concentration within the simulation.
        """
        self.gt.set_all_fugacities(self, fugacity)
        if recalc:
            self._calculate_K_r()
