#!/usr/bin/env python

import numpy as np
from hori import io
from ase.thermochemistry import IdealGasThermo, HarmonicThermo

def gas_chemical_potential(species,temperature,fugacity):
    """For a gas, given its temperature (K) and fugacity (Pa), returns
    its chemical potential."""
    if not gasprops.has_key(species):
        raise RuntimeError('Data for this gas (%s) is not in gasprops'
                           'dictionary. Add to hori.thermo.gasprops' %
                          species)

    d = io.gasdata(species)
    p = gasprops[species]
    electronic_energy = (d['electronic energy'] +
                         p['electronic energy correction'])
    thermo = IdealGasThermo(vib_energies=d['vibrations'],
                            electronicenergy=electronic_energy,
                            geometry=p['geometry'],
                            atoms=d['atoms'],
                            symmetrynumber=p['symmetrynumber'],
                            spin=p['spin'])
    mu = thermo.get_free_energy(temperature=temperature,pressure=fugacity)
    return mu

class ProtonElectronThermodynamics:
    """Class that loads the chemical potential data for the proton-electron
    pair in the Computational Hydrogen Electrode model, meaning that it
    uses the thermodynamic properties of H2 gas at 1 atm and adjusts them
    for temperature and voltage. In this module, the voltage is verus RHE,
    so a pH correction is not needed.

    Instantiate with the temperature and voltage, and the chemical
    potential of the proton-electron pair is automatically calculated and
    stored in self.G. Similarly, after changing the temperature or
    voltage with the .set_temperature or .set_voltage method, self.G
    is auto-recalculated.

    Temperature is in K and voltage is in V vs RHE.
    """

    def __init__(self,temperature,voltage):
        self.temperature = temperature
        self.voltage = voltage
        d = io.gasdata('H2') # hydrogen electrode
        d['electronic energy correction'] = \
                    gasprops['H2']['electronic energy correction']
        d['geometry'] = gasprops['H2']['geometry']
        d['symmetrynumber'] = gasprops['H2']['symmetrynumber']
        d['spin'] = gasprops['H2']['spin']
        self.referencepressure = 101325. # Pa, RHE definition
        self._previous = {'temperature':None,'voltage':None}
        thermo = IdealGasThermo(vib_energies=d['vibrations'],
                                electronicenergy=(d['electronic energy'] +
                                       d['electronic energy correction']),
                                geometry=d['geometry'],
                                atoms=d['atoms'],
                                symmetrynumber=d['symmetrynumber'],
                                spin=d['spin'])
        self._thermo = thermo
        self._calculate_G()

    def set_temperature(self,temperature):
        """Change the temperature (K)."""
        self._previous['temperature'] = self.temperature
        self.temperature = temperature
        self._calculate_G()

    def set_voltage(self,voltage):
        """Change the potential (V vs RHE)."""
        self._previous['voltage'] = self.voltage
        self.voltage = voltage
        self._calculate_G()

    def _calculate_G(self):
        """Calculates the free energy of the proton-electron pair. Only
        recalculates those portions that would change."""
        if self.temperature != self._previous['temperature']:
            self._G_H2 = self._thermo.get_free_energy(self.temperature, 
                                        self.referencepressure, 
                                        verbose=False) 
        if self.voltage != self._previous['voltage']:
            self._eU = 1. * self.voltage
        self.G = 0.5 * self._G_H2 - self._eU

class AdsorbateThermodynamics:
    """Class that loads adsorbate free energy data for a list of adsorbates
    and has methods to report back thermodynamic properties of the list as
    a function of the temperature and the electronic energy (without
    needing to re-access the data stored on disk). Initialize with
    ads_list, which is a list in one of two formats. It is either a list of
    the adsorbate names, e.g.,

        ['CHO', 'COOH', 'CO']

    or a list of adsorbates and a vibration mode:

        [['CHO', 'generic'], ['COOH', 'specific'], ['CO', 'auto'],]

    where 'generic' means it will calculate vibrations using the vibrations
    in the generic folder, 'specific' means it will use the specific
    vibrations for that adsorbate on the specified surface (and exit with
    an error if it is missing) and 'auto' means it will attempt to use
    'specific' vibrations, but will use 'generic' in each instance that the
    'specific' vibration is not available. If the first input form is used,
    then the vibrations mode is set with the 'vibrationsmode' keyword,
    which defaults to 'auto'. For brevity the numbers 0, 1, and 2 can be
    used for 'generic', 'specific', and 'auto', respectively.

    The adsorbate name should be something like 'COOH'. Use '' to indicate
    the clean slab. The other inputs are the temperature (in K) and a
    dictionary of hydrogen-bond stabilizations by the water overlayer, of
    the form hbond_dict[adsorbate] = float (in eV).

    When this class has all needed info it auto-calculates. When changing 
    an electronic energy or the temperature, the class auto-recalculates, 
    so a free energy can just be pulled with self.G[ads].

    If other properties are needed (mainly the electronic energy), they can
    be accessed in self.data[ads].
    """

    def __init__(self, ads_list, temperature, hbond_dict=None, 
                 vibrationsmode=None, surface=None):
        self.temperature = temperature
        self.surface = surface

        data = {}
        if hbond_dict == None:
            print('No hbond_dict specified.')

        if isinstance(ads_list[0], (list,tuple)):
            if vibrationsmode:
                raise RuntimeError('vibrationsmode doubly specified.')
            for ads in ads_list:
                if ads[1] is 'generic' or ads[1] == 0:
                    data[ads[0]] = {'vibrequest':'generic'}
                elif ads[1] is 'specific' or ads[1] is 1:
                    data[ads[0]] = {'vibrequest':'specific'}
                elif ads[1] is 'auto' or ads[1] is 2:
                    data[ads[0]] = {'vibrequest':'auto'}
        else:
            if vibrationsmode is None:
                vibrationsmode = 'auto'
            for ads in ads_list:
                data[ads] = {'vibrequest':vibrationsmode}

        for ads in data.keys():
            if hbond_dict is None:
                data[ads]['hbond'] = 0.
            else:
                data[ads]['hbond'] = hbond_dict[ads]
            data[ads]['E'] = np.nan
            data[ads]['EtoG'] = np.nan
            data[ads]['G'] = np.nan
            data[ads]['vibs'] = None

        self._changed = {'vibs': True,
                         'E': True,
                         'hbonds': True,
                         'T': True} # keep track of what needs calculation
        self.data = data
        self.G = {} # direct user access to free energies
        if self.surface is not None:
            self._populate_elec_energies()
            self._populate_vibrations()
            self.recalculate()

    def recalculate(self, force=False):
        """Calculates free energy information only when things have
        changed and the data (i.e., vibrations) exist. If force is True,
        then it will recalculate everything regardless of changes."""
        if sum(self._changed.values()) is False and force is False:
            return
        if self._changed['vibs'] or self._changed['T']:
            for ads in [ads for ads in self.data.keys() if
                        self.data[ads]['vibs'] is not None]:
                self._calculate_Gcorr(adsorbate=ads)
        for ads in self.data.keys():
            self._calculate_G(ads)
        for key in self._changed.keys():
            self._changed[key] = False

    def _populate_elec_energies(self):
        """If the surface is specified, reads in electronic energies from
        the data files. Note that the io function will put a numpy nan in
        this spot if the data does not exist."""
        data = self.data
        for ads in self.data.keys():
            if self.surface is not None:
                data[ads]['E'] = io.electronicenergy(surface=self.surface,
                                                     adsorbate=ads)
        self._changed['E'] = True

    def _populate_vibrations(self, adsorbate='all',
                             vibrations_mode = 'as specified'):
        """Adds generic or specific vibrations to each adsorbate if 'all'
        is selected, otherwise adds to the individual adsorbate listed.
        Will used the type of vibration specified for the adsorbate in
        self.data[ads]['vibrequest'] unless otherwise specified. If there
        are specific vibration requests and no surface has been specified,
        this function just populates *all* vibrations with None."""
        if adsorbate == 'all':
            if (self.surface is None and 'specific' in
                [self.data[key]['vibrequest'] for key in self.data.keys()]):
                return
            for ads in self.data.keys():
                mode = self.data[ads]['vibrequest']
                self._populate_vibrations(adsorbate=ads, vibrations_mode=mode)
        elif adsorbate == '':
            self.data[adsorbate]['vibs'] = []
        else: # specific adsorbates
            if vibrations_mode == 'as specified':
                vibrations_mode = self.data[adsorbate]['vibrequest']
            if vibrations_mode == 'generic':
                vibs = io.genericvibrations(adsorbate)
            elif vibrations_mode == 'specific':
                vibs = io.specificvibrations(adsorbate, self.surface)
            elif vibrations_mode == 'ignoreimaginary':
                vibs = io.ignoreimaginaryvibrations(adsorbate)
            elif vibrations_mode == 'auto':
                try:
                    vibs = io.specificvibrations(adsorbate, self.surface)
                except RuntimeError:
                    print('No specific adsorptions for %s on %s. Using '
                          'generic' % (adsorbate, self.surface))
                    vibs = io.genericvibrations(adsorbate)
            self.data[adsorbate]['vibs'] = vibs
        self._changed['vibs'] = True

    def _calculate_Gcorr(self, adsorbate, verbose=False):
        """Re-calculates the G_correction for the specified adsorbate. Note
        since only harmonic degrees of freedom are considered and the clean
        slab is considered to be static (i.e., its vibrations are
        unaffected by the presence of the adsorbate), then its G correction
        is 0."""
        if adsorbate == '':
            self.data[adsorbate]['EtoG'] = 0.
        else:
            vibs = self.data[adsorbate]['vibs']
            thermo = HarmonicThermo(vib_energies=vibs,
                                    electronicenergy=None)
            self.data[adsorbate]['EtoG'] = thermo.get_free_energy(
                temperature=self.temperature, verbose=verbose)

    def _calculate_G(self, adsorbate, verbose=False):
        """Recalculates the free energy for an adsorbate from
        pre-calculated parameters."""
        d = self.data[adsorbate]
        d['G'] = (d['E'] + d['EtoG'] + d['hbond']) 
        self.G[adsorbate] = d['G']
        if verbose:
            fmt = '%20s%13.3f eV'
            print('='*36)
            print('Contributions to G for %s' % adsorbate)
            print('-'*36)
            print(fmt % ('Electronic energy',self.data[adsorbate]['E']))
            print(fmt % ('Hydrogen bonding',self.data[adsorbate]['hbond']))
            print(fmt % ('E -> G', self.data[adsorbate]['EtoG']))
            print('-'*36)
            print(fmt % ('G', self.G[adsorbate]))
            print('='*36)

    def _scale_E(self, adsorbate, reference, fxn):
        """Provides a scaled value of the electronic energy. adsorbate is
        the name of an adsorbate (e.g., 'COOH') and reference is the name
        of the reference (e.g., 'CH3'). fxn is the function that links the
        two, as in:
            E[adsorbate] = fxn(E[reference])
        it is assumed that fxn uses referenced values for all energies,
        that is E = E_elec[ads*] - E_elec[*] - E_ref(ads)
        where E_ref is the atomic reference energies of the adsorbate.
        Returns the adsorbate energy with normal reference, 
        e.g., -55232.24 eV.
        Note that the electronic energy of the reference should be already
        set before this function is called.
        """
        #FIXME: this function is so far only used within an external
        # script, I may want to change the approach and delete it, for
        # example if I change the reference to pull out the atomic
        # normalizing (I need to leave in the clean slab normalizing).
        # FIXME: I will need to update these to the new data form next
        # time I use it.
        x = (self._elec_energy[reference] 
             - io.electronicenergy(self._surface)
             - calculate_reference_energ(referenc))
        y = fxn(x)
        self._elec_energy[adsorbate] = y + \
             io.electronicenergy(self._surface) + \
             calculate_reference_energy(adsorbate)
        # Make sure to add in a self._changed['E'] and a
        # self.recalculate()

    def set_temperature(self,temperature):
        """Sets the temperature (K) for adsorbate free energy calcs."""
        self.temperature = temperature
        self._changed['T'] = True
        self.recalculate()

    def set_electronic_energy(self, adsorbate, elec_energy):
        """Sets the electronic energy (eV) for the specified adsorbate and
        re-calculates G based only on that change (doesn't trigger a global
        recalculate)."""
        self.data[adsorbate]['E'] = elec_energy
        self._calculate_G(adsorbate)

    def set_surface(self, surface):
        """Looks up the electronic energies, and, if specified in
        vibrequest, vibrations, for each adsorbate on the 
        surface specified using the io functions of the hori module, 
        assigns them, and re-calculates G."""
        for ads in self.data.keys():
            self.data[ads]['E'] = io.electronicenergy(surface=surface,
                                                      adsorbate=ads)
        self.surface = surface
        self._populate_vibrations()
        self._changed['E'] = True
        self._changed['vibs'] = True
        self.recalculate()

    def set_energies_by_metal(self, metal):
        """Legacy function -- use set_surface."""
        self.set_surface(metal)

    def print_thermodynamic_contributions(self, ads_list='all'):
        """Prints the contributions to U,S,G for each adsorbate in ads_list
        at the current conditions.  If gas_list is 'all', then prints for 
        all gases. Note that this forces a re-calculation -- if any
        parameters have changed, then what is reported will be for those
        new values.
        """
        print('='*75)
        print('Thermodynamic contributions for adsorbates:')
        print('-'*75)
        if ads_list == 'all':
            ads_list = self.data.keys()
        for ads in ads_list:
            print('*'*70)
            print(ads)
            print('*'*70)
            if ads == '':
                print('No G-correction for clean slab.')
                print('Electronic energy: %.3f eV' % self.data['']['E'])
                print('Free energy:       %.3f eV' % self.G[''])
            else:
                self._calculate_Gcorr(ads,verbose=True)
                self._calculate_G(ads,verbose=True)
            print('\n\n')
        print('='*75)




class GasThermodynamics:
    """Class that loads gas chemical potential data for a list of gases,
    and has methods to report back thermodynamic properties of the list
    of gases as a function of temperature and fugacity (without needing
    to re-access the data stored on disk.
    
    Initialize with gas_list, which is a list of gases and their
    corresponding temperature and pressure (fugacity), of the form:

        [ [gas_name, temperature, fugacity], ... ]
    
    After instantiating this class or changing a temperature or pressure,
    this class auto-recalculates, so a free energy can just be pulled with
    self.G[gas]
    """

    def __init__(self,gas_list):
        gas_dict = {}
        for gas,temperature,fugacity in gas_list:
            gas_dict[gas] = io.gasdata(gas)
            if not gasprops.has_key(gas):
                raise RuntimeError('Data for this gas (%s) is not in '
                                   'gasprops dictionary. Add to '
                                   'hori.thermo.gasprops' % gas)
            gas_dict[gas]['geometry'] = gasprops[gas]['geometry']
            gas_dict[gas]['symmetrynumber'] = \
                    gasprops[gas]['symmetrynumber']
            gas_dict[gas]['spin'] = gasprops[gas]['spin']
            gas_dict[gas]['electronic energy correction'] = \
                    gasprops[gas]['electronic energy correction']
            gas_dict[gas]['temperature'] = temperature
            gas_dict[gas]['fugacity'] = fugacity
        self.gas_dict = gas_dict
        self.G = {}
        self._calculate_all_Gs()

    def _calculate_one_G(self,gas,verbose=False):
        d = self.gas_dict[gas]
        thermo = IdealGasThermo(vib_energies=d['vibrations'],
                                electronicenergy=(d['electronic energy'] +
                                       d['electronic energy correction']),
                                geometry=d['geometry'],
                                atoms=d['atoms'],
                                symmetrynumber=d['symmetrynumber'],
                                spin=d['spin'])
        self.G[gas] = thermo.get_free_energy(d['temperature'], 
                                             d['fugacity'],verbose=verbose)

    def _calculate_all_Gs(self):
        for gas in self.gas_dict.keys():
            self._calculate_one_G(gas)

    def set_gas_temperature(self,gas,temperature):
        """Change the temperature of an individual gas."""
        self.gas_dict[gas]['temperature'] = temperature
        self._calculate_one_G(gas)

    def set_gas_fugacity(self,gas,fugacity):
        """Change the fugacity of an individual gas."""
        self.gas_dict[gas]['fugacity'] = fugacity
        self._calculate_one_G(gas)

    def set_all_temperatures(self,temperature):
        """Set the temperature of all gases to the same value."""
        for gas in self.gas_dict.keys():
            self.gas_dict[gas]['temperature'] = temperature
        self._calculate_all_Gs()

    def set_all_fugacities(self,fugacity):
        """Set the fugacity of all gases to the same value."""
        for gas in self.gas_dict.keys():
            self.gas_dict[gas]['fugacity'] = fugacity 
        self._calculate_all_Gs()

    def print_thermodynamic_contributions(self,gas_list='all'):
        """Prints the contributions to H,S,G for each gas in gas_list at
        the current conditions.  If gas_list is 'all', then prints for all 
        gases.
        """
        if gas_list == 'all':
            gas_list = self.gas_dict.keys()
        for gas in gas_list:
            print('*'*60)
            print(gas)
            print('*'*60)
            print('Electronic energy contains %+.3f eV correction.' %
                  self.gas_dict[gas]['electronic energy correction'])
            self._calculate_one_G(gas=gas,verbose=True)
            print('\n\n')
                            


# Immutable gas properties follow.
gasprops = {
    'CO2':{'geometry':'linear',
           'symmetrynumber':2,
           'spin':0,
           'electronic energy correction':+0.45},
    'CH4':{'geometry':'nonlinear',
           'symmetrynumber':12,
           'spin':0,
           'electronic energy correction':+0.00},
    'H2': {'geometry':'linear',
           'symmetrynumber':2,
           'spin':0,
           'electronic energy correction':+0.00},
    'HCOOH' : {'geometry':'nonlinear',
               'symmetrynumber':1.,
               'spin':0.,
               'electronic energy correction':+0.45},
    'CO' : {'geometry':'linear',
            'symmetrynumber':1.,
            'spin':0.,
            'electronic energy correction':+0.00},
    'H2O' : {'geometry':'nonlinear',
             'symmetrynumber':2.,
             'spin':0.,
             'electronic energy correction':+0.00},
    'CH2O' : {'geometry' : 'nonlinear',
              'symmetrynumber' : 2.,
              'spin' : 0.,
              'electronic energy correction' : +0.00},
    'C2H4' : {'geometry' : 'nonlinear',
              'symmetrynumber' : 4.,
              'spin' : 0.,
              'electronic energy correction' : +0.00},
    'C3H7COOH' : {'geometry' : 'nonlinear',
                  'symmetrynumber' : 1.,
                  'spin' : 0.,
                  'electronic energy correction' : +0.00},
    'C3H8' : {'geometry' : 'nonlinear',
              'symmetrynumber' : 2.,
              'spin' : 0.,
              'electronic energy correction' : +0.00},
    'CH3OH' : {'geometry' : 'nonlinear',
               'symmetrynumber' : 1.,
               'spin' : 0.,
               'electronic energy correction' : 0.},
    'OCHCHO': {'geometry': 'nonlinear',
               'symmetrynumber' : 2.,
               'spin' :0.,
               'electronic energy correction' : +0.00},
    'CH3CH2OH': {'geometry': 'nonlinear',
                 'symmetrynumber' : 1.,
                 'spin' :0.,
                 'electronic energy correction' : +0.00},
    'OHCH2CH2OH': {'geometry' : 'nonlinear',
                   'symmetrynumber' : 1., 
                   'spin' : 0.,
                   'electronic energy correction' :0.00},
    'OCHCH2OH': {'geometry' : 'nonlinear',
                 'symmetrynumber' : 1.,
                 'spin' : 0.,
                 'electronic energy correction' : 0.00}

}
