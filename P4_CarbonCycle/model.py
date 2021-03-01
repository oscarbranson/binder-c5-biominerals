'''
Simple ocean carbonate model from Zeebe and Westbroek (2003; doi:10.1029/2003GC000538)
'''
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from ipywidgets import interact, interact_manual

# global variables
r = 4.3
c_n = 0.62  # neritic aragonite saturation
c_s = 1  # calcite saturation at z = 1 km
c_b = 1.73  # pelagic biogenic saturation 
c_d = 2.83  # calcite saturation at z = 6.5 km
c_c = 4.0  # critical inorganic supersaturation
Delta = c_d - c_s

# fluxes
F_in = 1  # riverine flux
F_n0 = 1 * F_in  # Biogenic neritic flux parameter
F_pm = r * F_in  # Biogenic pelagic flux, maximum

# constants
gamma = 500  # speed of abiotic precipitation above saturation limit

# converters
def c_to_CO3(c):
    return c * 52

def CO3_to_c(CO3):
    return CO3 / 52

# Fluxe Expressions
def strangelove_F_out(c, c_c=4, gamma=500):
    F_out = gamma * (c - c_c)**3
    F_out[c < c_c] = 0
    
    return F_out    

def neritic_F_out(c, F_n0, c_n, c_s):
    F_out = F_n0 * (c - c_n)**2 / c_s**2
    F_out[c < c_n] = 0
    return F_out

def cretan_F_out(c, F_pm, c_b):
    F_out = F_pm * c / c_b
    F_out[c >= c_b] = F_pm
    return F_out

def cretan_diss(c, F_pm, c_b, c_s, c_d, Delta):
    diss = np.zeros(c.size)
    ind = c < c_s
    diss[ind] = - F_pm * c[ind] / c_b
    ind = (c_s <= c) & (c < c_b)
    diss[ind] = (-F_pm * c[ind] / c_b) * (c_d - c[ind]) / Delta
    ind = (c_b <= c) & (c <= c_d)
    diss[ind] = -F_pm * (c_d - c[ind]) / Delta
    ind = c > c_d
    diss[ind] = 0
    
    return diss

# Steady State
def strangelove_c(F_in=1, c_c=4.0, gamma=500):
    return c_c + (F_in / gamma)**(1/3)

def neritic_c(F_in, F_n0, c_n, c_s):
    return c_n + c_s * np.sqrt(F_in / F_n0)

def cretan_c(r, f, c_s, c_b, c_d):
    Delta = c_d - c_s
    fprime = f / (1 - f)
    rinv = 1 / r
    c_D2 = (c_d - Delta * (1 + fprime)) / 2 + np.sqrt(((c_d - Delta * (1 + fprime)) / 2)**2 + rinv * c_b * Delta)
    c_D3 = c_d - Delta * (1 + fprime - r**-1)
    
    c = c_D2
    c[c > c_b] = c_D3[c > c_b]
    c[c < c_s] = np.nan  # Domain 1 behaviour
    
    return c

# Interactive plots
def _plot_strangelove_F(c_c=4, gamma=500):
    c = np.linspace(0, 5)
    plt.plot(c, strangelove_F_out(c, c_c, gamma), color='C0')
    plt.xlabel('Deep Ocean $c$')
    plt.ylabel('$F_{out}$')
    plt.axvline(c_c, ls='dashed', color=(0,0,0,0.5))
    plt.xlim(0,5)
    
def plot_strangelove_F():
    interact(_plot_strangelove_F,
             c_c=widgets.FloatSlider(value=c_c, min=0., max=5, step=0.01, description='$c_c$'),
             gamma=widgets.FloatSlider(value=gamma, min=100., max=1000, step=50, description='$\gamma$'))

def _plot_neritian_F(F_n0, c_n):
    c = np.linspace(0, 5, 500)
    neritic = neritic_F_out(c, F_n0, c_n, c_s)
    plt.plot(c, neritic, label='Neritian', color='C1')
    plt.plot(c, strangelove_F_out(c), alpha=0.5, label='Stangelove', color='C0')
    plt.ylim(-0.05 * neritic.max(), neritic.max())
    plt.xlabel('Deep Ocean $c$')
    plt.ylabel('$F_{out}$')
    plt.axhline(F_n0, ls='dashed', color=(0,0,0,0.5))
    plt.axvline(c_s + c_n, ls='dashed', color=(0,0,0,0.5))
    plt.xlim(0,5)
    plt.legend()

def plot_neritian_F():
    interact(
        _plot_neritian_F,
        F_n0=widgets.FloatSlider(1., min=0., max=5., step=0.1, description='$F_n^0$'),
        c_n=widgets.FloatSlider(c_n, min=0, max=2, step=0.01, description='$c_n$')
        )

def _plot_cretan_F(F_pm, c_b):
    c = np.linspace(0, 5, 500)
    cretan = cretan_F_out(c, F_pm, c_b)
    cretan_d = cretan_diss(c, F_pm, c_b, c_s, c_d, c_d - c_s)

    plt.plot(c, cretan + cretan_d, label='Pelagic', color='C2', lw=3)
    plt.plot(c, cretan, label='Pelagic Dep.', color='C2')
    plt.plot(c, cretan_d, label='Pelagic Diss.', color='C2', ls='dashed')

    plt.plot(c, neritic_F_out(c, F_n0, c_n, c_s), label='Neritian', color='C1', alpha=0.5)
    plt.plot(c, strangelove_F_out(c), alpha=0.5, label='Stangelove', color='C0')
    
    pad = 0.05 * (cretan.max() - cretan_d.min())
    plt.ylim(cretan_d.min() - pad, cretan.max() + pad)
    plt.xlabel('Deep Ocean $c$')
    plt.ylabel('$F_{out}$')
    # plt.axhline(F_n0, ls='dashed', color=(0,0,0,0.5))
    # plt.axvline(c_s + c_n, ls='dashed', color=(0,0,0,0.5))
    plt.xlim(0,5)
    plt.legend(loc='lower right')

def plot_cretan_F():
    interact(
        _plot_cretan_F,
        F_pm=widgets.FloatSlider(F_pm, min=0.5 * F_pm, max=2 * F_pm, step=0.1, description='$F_p^m$'),
        c_b=widgets.FloatSlider(c_b, min=1, max=4, step=0.01, description='$c_b$')
        )

def _plot_ss(f, c_c, c_b, c_s, c_n, show_modern):
    fig, ax = plt.subplots(1,1)
    
    F_in = 1
    sF_b = np.linspace(1, 10)
    F_pm = sF_b * (1 - f)
    r = F_pm / F_in

    cs = strangelove_c(F_in, c_c, gamma)
    cn = neritic_c(F_in, F_n0, c_n, c_s)
    cc = cretan_c(r, f, c_s, c_b, c_d)
    
    # strangelove
    ax.plot([0, F_in], [cs, cs], c='C0', lw=3, label='Strangelove')

    # neritic
    ax.plot([F_in, F_in], [c_n, cs], c='C1', lw=3, label='Neritian')

    # cretan
    ax.plot(sF_b, cc, c='C2', lw=3, label='Cretan')
    if any(np.isnan(cc)):
        x = np.max(sF_b[~np.isnan(cc)])
        y1 = np.nanmin(cc)
        ax.plot([x, x], [y1, c_n], color='C2', lw=3)

    ax.axhline(c_n, ls='dashed', color=(0,0,0,0.3), lw=1, zorder=-3)
    ax.axhline(c_s, ls='dashed', color=(0,0,0,0.3), lw=1, zorder=-3)
    ax.axhline(c_b, ls='dashed', color=(0,0,0,0.3), lw=1, zorder=-3)
    ax.axhline(c_d, ls='dashed', color=(0,0,0,0.3), lw=1, zorder=-3)
    ax.axhline(c_c, ls='dashed', color=(0,0,0,0.3), lw=1, zorder=-3)
    
    if show_modern:
        ax.axhspan(*CO3_to_c(np.array([71, 85])), color='C3', alpha=0.4, label='Modern')
    plt.legend()

    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)

    ax.set_xticks([0,1,2,4,6,8,10])
    ax.set_xticklabels([0, '$F_{in}$', 2, 4, 6, 8, 10])
    
    lax = ax.twinx()
    lax.set_ylim(0,5)
    lax.set_yticks([c_n, c_s, c_b, c_d, c_c])
    lax.set_yticklabels(['$c_n$', '$c_s$', '$c_b$', '$c_d$', '$c_c$'])
    
    ax.set_xlabel('Total Biogenic Flux ($F_b$)')
    ax.set_ylabel('Steady-State $c$')

def plot_ss():
    interact(_plot_ss,
            f=widgets.FloatSlider(value=0., min=0., max=0.99, step=0.01, description='f'),
            c_c=widgets.FloatSlider(value=c_c, min=c_d, max=5., step=0.01, description='$c_c$'),
            c_s=widgets.FloatSlider(value=c_s, min=c_n, max=c_b, step=0.01, description='$c_s$'),
            c_b=widgets.FloatSlider(value=c_b, min=c_s, max=c_d, step=0.01, description='$c_b$'),
            c_n=widgets.FloatSlider(value=c_n, min=0., max=c_s, step=0.01, description='$c_n$'),
            show_modern=widgets.Checkbox(False, description='Modern Ocean')
            )

