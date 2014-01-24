import numpy as np
#import scipy
from scipy.optimize import fsolve
from scipy.integrate import odeint
from scipy.integrate import ode 

def constraint_matrices(NS, constraints):
    """ Create matrices to map between the complete system and system constrained by
    site conservation.
    In the 'full' system all sites including free sites are treated explicitly.
    In the constrained systems some coverages are eliminated by site conservation.
    NS is the number of surface species in the full system.
    'constraints' is a list of constraints of the type
    ((ic, (i1, i2, ... iN), (s1, s2, ... sN), sum_theta), ...)
    where the constraint
    sum_theta = s1 * theta_i1 + s1 * theta_i2 + ... s1 * theta_iN
    is replaceing the ode for the ic'th coverage.
    ic should be included in (i1, i2, ... iN).
    The transformation is as follows:
    y_full = Msc * y_constrained + b
    y_constrained = Mcs * y_full
    """
    NC = NS - len(constraints) # DOF of constrained system
    Msc = np.zeros((NS, NC))   #
    Mcs = np.zeros((NC, NS))   #
    b = np.zeros(NS)           #
    # figure out which coverages are constrained
    constrained_i = [] # list with constrained coverages in the full system
    for c in constraints:
        if c[0] not in constrained_i:
            constrained_i.append(c[0])
        else:
            print 'Warning: multiple constraints specified for the same coverage.'
    # dict mapping an independent coverage to the constrained system
    # (y1, y2, y3, y4) -> (y1, y2, y4)
    independent_i = {} 
    ii = 0 
    for i in range(NS):
        if i in  constrained_i:
            pass
        else:
            independent_i[i] = ii
            ii += 1
    # set up b
    for c in constraints:
        b[c[0]] = c[3]
    # set up Msc
    ic_count = 0 # count how many 'independent' coverages have been encountered
    for i in range(NS): # loop over rows of MSc
        constrained = False
        for c in constraints:
            if c[0] == i:
                constrained = True
                for j in range(len(c[1])):
                    if c[1][j] in constrained_i: # constrained coordinate, skip
                        pass
                    else:
                        ni = independent_i[c[1][j]] # index in the reduced system
                        Msc[i, ni] = - c[2][j]
        if not constrained:
            Msc[i, ic_count] = 1
            Mcs[ic_count, i] = 1
            ic_count += 1            
    return Mcs, Msc, b 

def pack(transformation, y):
    """ Eliminate dependent coverages using the constraint defined by 'transformation'.
    """
    return np.dot(transformation[0], y)

def unpack(transformation, y):
    """ Use the constraints defined by 'transformation' to get the complete set of coverages
    from the site conservation
    """
    return np.dot(transformation[1], y) + transformation[2]

class kinetic_solver:
    """ Base class for the kinetic solver.
    """
    def __init__(self):
        assert 0**0 == 1 # 
        pass

    def initialize(self, NG, NS, fugacities, stoich_mat, k_r, K_r):
        """ Initialize the kinetic solver.
        NG - number of non-adsorbed species
        NS - number of surface species 
        stoich_mat - the stoichiometric matrix
        fugacities - a list of fugacities (must match the order in stoich_mat).
        k_r - forward rate constants
        K_r - equilibrium constant 
        """
        assert len(fugacities) == NG
        assert stoich_mat.shape[1] == NG + NS
        self.NG = NG
        self.NS = NS
        self.NR = stoich_mat.shape[0]
        self.stoich_mat = stoich_mat
        self.fugacities = fugacities
        self.k_r = k_r
        self.K_r = K_r
        self.krev_r = k_r/K_r
        self._calc_kf_r()  # rate constants multiplied by fugacities.
        
    def _calc_kf_r(self):
        """ Calculate rate constants multiplied by fugacities.
        These are constants in the simulation.
        """
        kf_r =  np.zeros(self.NR)
        krevf_r =  np.zeros(self.NR)
        for r in range(self.NR):
            kf = self.k_r[r]
            krevf = self.krev_r[r]
            for c in range(self.NG):
                alpha = self.stoich_mat[r,c]
                f = self.fugacities[c]
                if alpha < -0.:
                    kf *= f**(-alpha)
                elif alpha > 0.:
                    krevf *= f**(alpha)
            kf_r[r] = kf
            krevf_r[r] = krevf
        self.kf_r = kf_r
        self.krevf_r = krevf_r

    def get_odes_c(self, coverages, t=0, verbose=0):
        """ Calculate d/dt on the coverages of surface species as function of the coverage.
        """
        if verbose > 0: print 'get_odes_c; c,t:',  coverages, t
        dcov_dt_c = np.zeros(self.NS)
        for c in range(self.NS):  # calc d/dt for all surface species
            for r in range(self.NR):    # loop over reactions
                r_forward  = self.kf_r[r]
                r_reverse = self.krevf_r[r]
                for c2 in range(self.NS):
                    alpha = self.stoich_mat[r,c2+self.NG]
                    theta = coverages[c2]
                    if alpha < 0.:
                        r_forward *= theta**(-alpha)
                    elif alpha > 0.:
                        r_reverse *= theta**(alpha)
                netrate_r = r_forward - r_reverse
                dcov_dt_c[c] += self.stoich_mat[r,c+self.NG] * netrate_r
        return dcov_dt_c

    def get_odes_tc(self, t, coverages):
        """ Calculate d/dt on the coverages of surface species as function of the coverage.
        """
        return self.get_odes_c(coverages, t)

    def get_jacobian_ci(self, coverages, t=0, verbose=0):
        """ Calculate the Jacobian d f_c / d theta_i
        The current implementation assumes the (effective) rate constants are independent
        of coverages.
        Has not been tested extensively.
        Could probably be optimized for speed.
        """
        if verbose > 0: print 'get_jacobian_ci; c,t:',  coverages, t
        J_ci = np.zeros((self.NS, self.NS))
        for i in range(self.NS):  # surface species
            for c in range(self.NS):  # surface species
                for r in range(self.NR):    # loop over reactions
                    J_forward = self.kf_r[r]
                    J_reverse = self.krevf_r[r]
                    for c2 in range(self.NS): # loop over theta_c2
                        theta = coverages[c2]
                        alpha = self.stoich_mat[r,c2+self.NG]
                        # the theta we differentiate wrt to
                        if i==c2:
                            #if alpha ==0:
                            #    pass
                            if alpha < 0 :
                                J_forward *= -alpha
                                alpha = alpha + 1
                                #J_forward *= theta**(-alpha)
                                J_reverse = 0.
                            elif alpha > 0 :
                                J_reverse *= alpha
                                alpha = alpha -1
                                #J_reverse *= theta**(alpha)
                                J_forward = 0.
                        # the other coverages
                        if alpha < 0.:
                            J_forward *= theta**(-alpha)
                        elif alpha > 0.:
                            J_reverse *= theta**(alpha)
                    # contribution corresponding to \partial fc / \partial thetai not implemented
                    J_ci[c,i] += self.stoich_mat[r,c+self.NG] * ( J_forward - J_reverse )
        return J_ci

    def get_rates(self, coverages):
        """ Calculate the rates for all gas phase and adsorbed species as function of coverages.
        The result is storedin self.d_dt_c
        """
        d_dt_c = np.zeros(self.NG+self.NS)
        for c in range(self.NG+self.NS):  # calc d/dt for all surface species
            for r in range(self.NR):    # loop over reactions
                r_forward  = self.kf_r[r]
                r_reverse = self.krevf_r[r]
                for c2 in range(self.NS):
                    alpha = self.stoich_mat[r,c2+self.NG]
                    theta = coverages[c2]
                    if alpha < 0.:
                        r_forward *= theta**(-alpha)
                    elif alpha > 0.:
                        r_reverse *= theta**(alpha)
                netrate_r = r_forward - r_reverse
                d_dt_c[c] += self.stoich_mat[r,c] * netrate_r
        self.d_dt_c = d_dt_c
        return d_dt_c

class integrate_to_time_t(kinetic_solver):
    """ Integrate the ODE's to through time t.
    No site conservation is explicitly invoked, so the initial guess must be chosen carefully.
    The ODE's should conserve the number of sites, yet this should probably be verified for long integration times.
    It may be a good idea to define a new class which explicitly takes site conservation into account.
    Based on scipy's odeint.
    """
    def solve(self, thetas0, ts=np.logspace(-15,10,30,base=10), rtol=1e-3, useJac=False, full_output=False):
        if useJac:
            Dfun = self.get_jacobian_ci
        else:
            Dfun = None
        theta_ts, infodict= odeint(self.get_odes_c, thetas0, ts, Dfun=Dfun,
                                   #ixpr=True,
                                   #printmessg=True,
                                   full_output=True
                                   )
        for i in range(self.NS):
            if abs(1 - theta_ts[-2,i] / theta_ts[-1,i]) > rtol :
                print "WARNING: ODE's have not achieved steady state."
        theta_c = theta_ts[-1,:]
        rates_c = self.get_rates(theta_c)
        #return theta_c, rates_c
        if full_output:
            return theta_c, rates_c, infodict
        else:
            return theta_c, rates_c,

class integrate_to_time_t_constrained(kinetic_solver):
    """ Integrate the ODE's to through time t. Use constraints to reduce the number of odes.
    Based on scipy's odeint.
    """
    def solve(self, thetas0_full, constraints, ts=np.logspace(-15,10,30,base=10), rtol=1e-3,
              full_output=False):
        transform = constraint_matrices(self.NS, constraints)
        thetas0 = pack(transform, thetas0_full)
        def func(c, t):
            cfull = unpack(transform, c)
            dcov_dt_c_full = self.get_odes_c(cfull, t)
            return pack(transform, dcov_dt_c_full)
        theta_ts, infodict = odeint(func, thetas0, ts, printmessg=False, full_output=True,
                                    rtol=1e-3, atol=1e-3,
                                    h0=1e-9, # first step length
                                    #ixpr=True,
                                    )
        for i in range(self.NS-len(constraints)):
            if abs(1 - theta_ts[-2,i] / theta_ts[-1,i]) > rtol :
                print "WARNING: ODE's have not achieved steady state."
        theta_c = theta_ts[-1,:]
        theta_c_full = unpack(transform, theta_c)
        rates_c_full = self.get_rates(theta_c_full)
        if full_output:
            return theta_c_full, rates_c_full, infodict
        else:
            return theta_c_full, rates_c_full

class integrate_to_time_t_vode(kinetic_solver):
    """ Integrate the ODE's to through time t.
    using scipy's vode integrator.
    Scipy's odeint seems to work better than this vode implementation.
    """
    def solve(self, thetas0, ts=np.logspace(-15,10,30,base=10), rtol=1e-3):
        t0 = 0.
        dts = np.diff(ts)
        theta_ts = np.zeros((len(ts), self.NS))
        r = ode(self.get_odes_tc) # takes t or theta as first argument?
        r.set_integrator('vode', method='bdf', order=3, nsteps=5000)
        r.set_initial_value(thetas0, t0)
        t_i = 0
        while r.successful() and r.t < ts[-1]:
            theta_ts[t_i, :] = r.integrate(r.t+dts[t_i])
            t_i += 1
        for i in range(self.NS):
            if abs(1 - theta_ts[-2,i] / theta_ts[-1,i]) > rtol :
                print "WARNING: ODE's have not achieved steady state."
        theta_c = theta_ts[-1,:]
        rates_c = self.get_rates(theta_c)
        return theta_c, rates_c

class steady_state_solver(kinetic_solver):
    """ Class for finding the steady state solution. 
    The implementation replaces one of the rate equations with the site conservation constraint:
    sum_i theta_i = 1.
    It may be desirable to implement more flexible constrains (similar to the integrators).
    """
    #import scipy.optimize as opt
    def func(self, covs):
        """ Solve N-1 coupled ODE's  with the constraint of site conservation.
        """
        y = self.get_odes_c(covs)
        y[0] = covs.sum() - 1. # site conservation.
        #y[-1] = covs.sum() - 1.
        return y
        
    def solve(self, theta_0):
        theta_c, infodict, iter, msg = fsolve(self.func, theta_0, full_output=True, xtol=1e-14)
        rates_c = self.get_odes_c(theta_c)
        return theta_c, rates_c
