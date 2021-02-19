import load
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import ipywidgets as widgets
from ipywidgets import interact


# Step Length
def Teng_LC():
    dat_LC = load.Teng_LC()
    # this loop looks at the acute and obtuse step data independently
    for tp in ['acute', 'obtuse']:
        ind = dat_LC.type == tp  # this creates a boolean index for selecting data
        plt.scatter(dat_LC.loc[ind, '1_sigma'], dat_LC.loc[ind, 'LC'], label=tp)  # this plots ths data separately

    # make it look nice
    plt.title('Critical Step Lengths for Calcite')
    plt.legend()
    plt.ylabel('$L_C$ (nm)')
    plt.xlabel('1 / $\sigma$')
    plt.xlim(0, 4.2)
    plt.ylim(0, 200)

def critLen(sigma_inv, gamma, sigma_c):
    """
    Calculates critical step length (LC) as a function of 1/sigma (sigma_inv, 
    the step energy (gamma) and the energy barrier to 1D nucleation (g1D).    
    """
    b = 6.4e-1  # calcite intermolecular distance along step (nm)
    c = 3.2e-1  # calcite intermolecular distance between rows (nm)
    kB = 8.617333262145e-5  # Boltzman constant eV K-1
    T = 298.15  # K
    
    # sigma_r_inv = 1 / (g1D / (kB * T))  # this converts energy barrier to 1D nuclation into inverse saturation.
    sigma_c_inv = 1 / sigma_c

    return 2 * b * c * gamma * (sigma_inv - sigma_c_inv) / (kB * T)

def _fit_critLen(gamma_ac, sigma_c_ac, gamma_ob, sigma_c_ob):

    dat_LC = load.Teng_LC()

    sigma_inv = np.linspace(0, 4.2)
    pred_LC_ac = critLen(sigma_inv, gamma_ac, sigma_c_ac)
    pred_LC_ob = critLen(sigma_inv, gamma_ob, sigma_c_ob)

    for tp in ['acute', 'obtuse']:
        ind = dat_LC.type == tp  # this creates a boolean index for selecting data
        plt.scatter(dat_LC.loc[ind, '1_sigma'], dat_LC.loc[ind, 'LC'], label=tp)  # this plots ths data separately

    plt.plot(sigma_inv, pred_LC_ac, color='C0')
    plt.plot(sigma_inv, pred_LC_ob, color='C1')

    plt.title('Critical Step Lengths for Calcite')
    plt.legend()
    plt.ylabel('$L_C$ (nm)')
    plt.xlabel('1 / $\sigma$')
    plt.xlim(0, 4.2)
    plt.ylim(0, 200)

def fit_critLen():
    interact(
        _fit_critLen,
        gamma_ac=widgets.FloatSlider(2.4, min=0, max=5, step=0.01, description='$\gamma_-$'),
        sigma_c_ac=widgets.FloatSlider(1, min=0.01, max=10, step=0.01, description='$\sigma_{c-}$'),
        gamma_ob=widgets.FloatSlider(2.4, min=0, max=5, step=0.01, description='$\gamma_+$'),
        sigma_c_ob=widgets.FloatSlider(1, min=0.01, max=10, step=0.01, description='$\sigma_{c+}$')
    )

# Step Velocity
def Teng_v():
    dat_v = load.Teng_v()

    for tp, m in zip(['acute', 'obtuse'], ['o', 's']):
        ind = dat_v.type == tp
        plt.scatter(dat_v.loc[ind, 'sigma'], dat_v.loc[ind, 'v'], marker=m, label=tp)
        
    plt.xlabel('$\sigma$')
    plt.ylabel('v (nm/s)')
    plt.legend()

def calc_vinf(omega, beta, activity, sigma):
    """
    Calculates the velocity of an infinitely long step as a function of
    sigma.
    """
    return omega * beta * activity * (np.exp(sigma) - 1)

def calc_LC_params():
    dat_LC = load.Teng_LC()
    params = {}  # a dictionary to store the best-fit parameters in
    for tp in ['acute', 'obtuse']:  # loop through the two step types
        ind = dat_LC.type == tp
        # opt.curve_fit gives us back two arrays: the optimised parameters, and the covariance 
        # matrix of those parameters. We store these below as 'p' and 'cov', and then save them to
        # the 'params' dictionary.
        p, cov = opt.curve_fit(critLen, dat_LC.loc[ind, '1_sigma'], dat_LC.loc[ind, 'LC'], (1., .1))
        params[tp] = p, cov
    return params

params = calc_LC_params()
dat_v = load.Teng_v()

def normalise_Teng_v():
    beta = 0.4e7  # calcite kinetic coefficient - measured range 0.3-0.5 (nm/s)
    Ksp = 10**-8.54  # equilibrium activity product (mol2 kg-2)
    b = 6.4e-1  # calcite intermolecular distance along step (nm)
    c = 3.2e-1  # calcite intermolecular distance between rows (nm)
    h = 3.1e-1  # calcite step height (nm)
    omega = c * b * h  # volume per growth unit (average) (nm3)

    # We need to treat the acute and obtuse steps differently, as they have different
    # step energies and therefore critical lengths.
    for tp in ['acute', 'obtuse']:
        ind = dat_v.type == tp
        # use our LC function from above
        dat_v.loc[ind, 'LC'] = critLen(1/dat_v.loc[ind, 'sigma'], *params[tp][0])

    # assume solution is equimolar, so activity of single species is Ksp**0.5
    dat_v.loc[:, 'Vinf'] = calc_vinf(omega, beta, Ksp**0.5, dat_v.loc[:, 'sigma'])

    # normalise measured data to theoretical maximum velocity and critical step length
    dat_v.loc[:, 'V/Vinf'] = dat_v.loc[:, 'v'] / dat_v.loc[:, 'Vinf']
    dat_v.loc[:, 'L/LC'] = dat_v.loc[:, 'L'] / dat_v.loc[:, 'LC']

normalise_Teng_v()

def Teng_v_normalised():

    for tp, m in zip(['acute', 'obtuse'], ['o', 's']):
        ind = dat_v.type == tp
        plt.scatter(dat_v.loc[ind, 'L/LC'], dat_v.loc[ind, 'V/Vinf'], label=tp, marker=m, c=dat_v.loc[ind, 'sigma'])

    plt.colorbar(label='$\sigma$')
        
    plt.xlabel('$L/L_C$')
    plt.ylabel('$v/v_{\infty}$')
    plt.legend()

# Gibbs-Thomson
def calc_v_vinf_GT(sigma, LC_L):
    return 1 - (np.exp(sigma * LC_L) - 1) / (np.exp(sigma) - 1)

def Teng_v_GT():
    L_LC = np.linspace(1, 4)
    LC_L = 1 / L_LC
    v_vinf = calc_v_vinf_GT(1, LC_L)

    for tp, m in zip(['acute', 'obtuse'], ['o', 's']):
        ind = dat_v.type == tp
        plt.scatter(dat_v.loc[ind, 'L/LC'], dat_v.loc[ind, 'V/Vinf'], label=tp, marker=m, c=dat_v.loc[ind, 'sigma'])

    plt.plot(L_LC, v_vinf, label='Gibbs-Thomson Prediction')
    plt.colorbar(label='$\sigma$')

    plt.xlabel('$L/L_C$')
    plt.ylabel('$v/v_{\infty}$')
    plt.legend()

# Terrace Width
dat_W = load.Teng_W()

def Teng_W():
    sigma = np.linspace(0.16, 1.5)

    for c, tp in zip(['C0', 'C1'], ['acute', 'obtuse']):
        ind = dat_W.type == tp
        LC = critLen(1/sigma, *params[tp][0])
        plt.scatter(dat_W.loc[ind, 'sigma'], dat_W.loc[ind, 'W'], color=c, label=tp)
        
        plt.plot(sigma, 4 * LC, color=c)
        plt.plot(sigma, 9.6 * LC, ls='dashed', color=c)
    
    for ls, label in zip(['solid', 'dashed'], ['isotropic', 'Gibbs-Thomson']):
        plt.plot([], [], color='k', ls=ls, label=label)
    
    plt.legend()
    plt.ylabel('W (nm)')
    plt.xlabel('$\sigma$')

# Mg inhibition
dat_Davis = load.Davis_v()

def _Davis_Mg_v(show_fitlines=False):
    fig, ax = plt.subplots(1,1)
    cb = ax.scatter(dat_Davis.loc[:, 'aCa'] * 1e6, dat_Davis.loc[:, 'v'], c=dat_Davis.loc[:, 'xMg'])
    fig.colorbar(cb, label='xMg')
    ax.set_xlabel(r'$\alpha Ca^{2+}$')
    ax.set_ylabel('v (nm/s)')

    if show_fitlines:
        aCa = np.linspace(60, 120)
        for xMg in dat_Davis.loc[:, 'xMg'].unique():
            ind = (dat_Davis.loc[:, 'xMg'] == xMg) & (dat_Davis.loc[:, 'aCa'] * 1e6 > 80)
            
            p = np.polyfit(dat_Davis.loc[ind, 'aCa'] * 1e6, dat_Davis.loc[ind, 'v'], 1)
            ax.plot(aCa, np.polyval(p, aCa), color=(0,0,0,0.6), ls='dashed')
    
    ax.set_xlim(50, 120)
    ax.set_ylim(0, 7)



def Davis_Mg_v():
    interact(
        _Davis_Mg_v,
        show_fitlines=widgets.Checkbox(False, description='Show Fits for single [Mg] levels')
    )