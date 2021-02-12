'''
Simple ocean carbonate model from Zeebe and Westbroek (2003; doi:10.1029/2003GC000538)
'''
import numpy as np
import matplotlib.pyplot as plt
import ipywidgets as widgets
from ipywidgets import interact, interact_manual


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
    
    return diss

# Steady State
def strangelove_c(F_in=1, c_c=4.0, gamma=500):
    return c_c + (F_in / gamma)**(1/3)

def neritic_c(F_in, F_n0, c_n, c_s):
    return c_n + c_s * np.sqrt(F_in / F_n0)

def cretan_c(F_in, F_pm, F_n0, f, c_s, c_b, c_d, c_n):
    Delta = c_d - c_s
    fprime = f / (1 - f)
    r = F_pm / F_in
    rinv = 1 / r
    c_D2 = (c_d - Delta * (1 + fprime)) / 2 + np.sqrt(((c_d - Delta * (1 + fprime)) / 2)**2 + rinv * c_b * Delta)
    c_D3 = c_d - Delta * (1 + fprime - r**-1)
    
    c = c_D2
    c[c > c_b] = c_D3[c > c_b]
    c[c < c_s] = np.nan  # Domain 1 behaviour
    
    return c

# Interactive plots
def _plot_ss(f, F_pm, F_n0, c_s, c_b, c_d, c_n, c_c, gamma):
    fig, ax = plt.subplots(1,1)
    
    F_in = 1
    sF_b = np.linspace(1, 10)
    F_pm = sF_b * (1 - f)

    cs = strangelove_c(F_in, c_c, gamma)
    cn = neritic_c(F_in, F_n0, c_n, c_s)
    cc = cretan_c(F_in, F_pm, F_n0, f, c_s, c_b, c_d, c_n)
    
    # strangelove
    ax.plot([0, F_in], [cs, cs])

    # neritic
    ax.plot([F_in, F_in], [c_n, cs])

    # cretan
    ax.plot(sF_b, cc)

    ax.axhline(c_n, ls='dashed', color=(0,0,0,0.4))
    ax.axhline(c_s, ls='dashed', color=(0,0,0,0.4))
    ax.axhline(c_b, ls='dashed', color=(0,0,0,0.4))
    ax.axhline(c_d, ls='dashed', color=(0,0,0,0.4))
    ax.axhline(c_c, ls='dashed', color=(0,0,0,0.4))
    
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    
    lax = ax.twinx()
    lax.set_ylim(0,5)
    lax.set_yticks([c_n, c_s, c_b, c_d, c_c])
    lax.set_yticklabels(['$c_n$', '$c_s$', '$c_b$', '$c_d$', '$c_c$'])
    
    ax.set_xlabel('Total Biogenic Flux ($F_b$)')
    ax.set_ylabel('Steady-State $[CO_3^{2-}]$')

def plot_ss():
    interact(plot_ss,
            f=widgets.FloatSlider(value=0., min=0., max=1., step=0.01, description='f'),
            F_pm=widgets.FloatText(value=1., min=0., description='$F_p^m$'), 
            F_n0=widgets.FloatText(value=1., min=0., description='$F_n^0$'), 
            c_c=widgets.FloatSlider(value=4., min=1., max=10., step=0.01, description='$c_c$'),
            c_d=widgets.FloatSlider(value=2.83, min=0., max=10., step=0.01, description='$c_d$'),
            c_b=widgets.FloatSlider(value=1.73, min=0., max=10., step=0.01, description='$c_b$'),
            c_s=widgets.FloatSlider(value=1., min=0., max=10., step=0.01, description='$c_s$'), 
            c_n=widgets.FloatSlider(value=0.62, min=0., max=10., step=0.01, description='$c_n$'),
            gamma=widgets.FloatText(value=500, min=0., step=10, description='$\gamma$')
            )