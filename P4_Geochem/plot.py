import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import ipywidgets as widgets
from ipywidgets import interact

def fn_Rayleigh(f, K):
    return (1 - f**K) / (1 - f)

def _Rayleigh(K):
    f = np.concatenate([[0], np.logspace(-6, -0.0000001,100)])
    
    fig, ax = plt.subplots(1, 1)
    
    ax.plot(f, fn_Rayleigh(f, K))
    ax.invert_xaxis()
    
    ax.set_xlabel('f (Fraction of Ca Remaining)')
    ax.set_ylabel('Distribution Coefficient')
    ax.set_ylim(0, 2)
    ax.set_xlim(1, 0)
    
def Rayleigh():
    interact(_Rayleigh, K=widgets.FloatSlider(min=0.01, max=2, value=1, step=0.01, decription='K'))
    
dat = pd.DataFrame.from_dict({
    'Element': ['Mg', 'Sr', 'Ba', 'Na', 'K', 'Li', 'Mn'],
    'Inorganic': [0.02, 0.11397333218331, 0.0296982050896016, 0.0000425, 0.000055, 0.003505, 17.5],
    'Foraminifera': [0.00126505111700114, 0.171333333333333, 0.176772486772487, 0.000131497036626444, 0.000201371204701273, 0.01000748, 30.3],
    # 'Kf': [0.03, 0.24, 0.0562341325190349, 7.40E-03, 7.40E-03, 1.50E-03, 5],
    # 'Keq': [0.005, 0.03, 0.00316227766016838, 6.30E-05, 6.30E-05, 5.00E-04, 50]
    'Kf': [0.06, 0.24, 0.08, 7.40E-03, 7.40E-03, 1.50E-03, 10],
    'Keq': [0.01, 0.06, 0.015, 1.80E-05, 2.00E-05, 5.00E-03, 60]
})

def _foraminifera_Rayleigh(log10f):
    fig, ax = plt.subplots(1, 1)
    
    for i, r in dat.iterrows():
        ax.scatter(r['Inorganic'], r['Foraminifera'] / r['Inorganic'], label=r['Element'])

    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.set_xlim(ax.get_xlim())
    ax.set_ylim(ax.get_ylim())
    
    ax.axhline(1, color='black', linestyle='--', alpha=0.5)
    ax.axvline(1, color='black', linestyle='--', alpha=0.5)
    
    K = np.logspace(-5, 2, 100)
    ax.plot(K, fn_Rayleigh(10**log10f, K) / K)
    
    ax.set_xlabel('Inorganic K')
    ax.set_ylabel('Foraminifera D / Inorganic K')
    
def foraminifera_Rayleigh():
    interact(_foraminifera_Rayleigh, log10f=widgets.FloatSlider(min=-6, max=-0.001, value=-0.001, step=0.001, decription='log10(f)'))
    
def fn_SKM(Rp, Kf, Keq, Rb):
    return Kf / (1 + Rb / (Rp + Rb) * (Kf / Keq - 1))

def _SKM(Kf, Keq, log10Rb):
    Rp = np.logspace(-9, -3, 100)
    
    fig, ax = plt.subplots(1, 1)
    
    ax.plot(Rp, fn_SKM(Rp, Kf, Keq, 10**log10Rb))
    
    ax.set_xscale('log')
    
    ax.set_xlim(10**-9, 10**-3)
    ax.set_ylim(0.01, 2)
    
    ax.axhline(Kf, color='black', linestyle='--', alpha=0.5, lw=1, zorder=-1)
    ax.text(1.2e-3, Kf, '$K_f$', va='center')
    ax.axhline(Keq, color='black', linestyle='--', alpha=0.5, lw=1, zorder=-1)
    ax.text(1.2e-3, Keq, '$K_{eq}$', va='center')

    ax.axvline(10**log10Rb, color='black', linestyle='--', alpha=0.5, lw=1, zorder=-1)
    ax.text(10**log10Rb, 2.05, '$R_b$', ha='center')
    
    ax.set_xlabel('Precipitation Rate ($mol/m^2/s$)')
    ax.set_ylabel('Partitioning Coefficient')
    
def SKM():
    interact(_SKM, 
             Kf=widgets.FloatSlider(min=0.01, max=2, value=1.5, step=0.01, decription='Kf'), 
             Keq=widgets.FloatSlider(min=0.01, max=2, value=0.5, step=0.01, decription='Keq'), 
             log10Rb=widgets.FloatSlider(min=-8, max=-5, value=-6.2, step=0.01, decription='$log_{10}R_b$'))
    
def _foraminifera_SKM(log10Rp, log10Rb):
    fig, ax = plt.subplots(1, 1)
    
    inorg = fn_SKM(Kf=dat['Kf'], Keq=dat['Keq'], Rp=10**log10Rp, Rb=10**log10Rb)
    
    for i, r in dat.iterrows():
        ax.scatter(inorg[i], r['Foraminifera'] / inorg[i], label=r['Element'])

    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    # ax.set_xlim(ax.get_xlim())
    # ax.set_ylim(ax.get_ylim())
    
    ax.axhline(1, color='black', linestyle='--', alpha=0.5)
    ax.axvline(1, color='black', linestyle='--', alpha=0.5)
    
    # K = np.logspace(-5, 2, 100)
    # ax.plot(K, fn_Rayleigh(10**log10f, K) / K)
    
    ax.set_xlabel('Inorganic K')
    ax.set_ylabel('Foraminifera D / Inorganic K')
    
def foraminifera_SKM():
    interact(_foraminifera_SKM,
            log10Rp=widgets.FloatSlider(min=-9, max=-3, value=-6, step=0.01, decription='$log_{10}R_p$'),
            log10Rb=widgets.FloatSlider(min=-8, max=-5, value=-6.6, step=0.01, decription='$log_{10}R_b$')
            )
    
def _foraminifera_SKM_Rayleigh(log10Rp, log10Rb, log10f):
    fig, ax = plt.subplots(1, 1)
    
    inorg = fn_SKM(Kf=dat['Kf'], Keq=dat['Keq'], Rp=10**log10Rp, Rb=10**log10Rb)
    
    for i, r in dat.iterrows():
        ax.scatter(inorg[i], r['Foraminifera'] / inorg[i], label=r['Element'])

    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    ax.set_xlim(ax.get_xlim())
    # ax.set_ylim(ax.get_ylim())
    
    ax.axhline(1, color='black', linestyle='--', alpha=0.5)
    ax.axvline(1, color='black', linestyle='--', alpha=0.5)
    
    K = np.logspace(-5, 2, 100)
    ax.plot(K, fn_Rayleigh(10**log10f, K) / K)
    
    ax.set_xlabel('Inorganic K')
    ax.set_ylabel('Foraminifera D / Inorganic K')
    
def foraminifera_SKM_Rayleigh():
    interact(_foraminifera_SKM_Rayleigh,
            log10Rp=widgets.FloatSlider(min=-9, max=-3, value=-6, step=0.01, decription='$log_{10}R_p$'),
            log10Rb=widgets.FloatSlider(min=-8, max=-5, value=-6.6, step=0.01, decription='$log_{10}R_b$'),
            log10f=widgets.FloatSlider(min=-6, max=-0.001, value=-0.001, step=0.001, decription='log10(f)')
            )