import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import ipywidgets as widgets
from ipywidgets import interact, interact_manual

# Rule of mixtures fibre composite

def calc_axial(Vf, Ef, Em):
    return Ef * Vf + Em * (1 - Vf)

def calc_transverse(Vf, Ef, Em):
    return 1 / ((1 - Vf) / Em + Vf / Ef)

def calc_halpintsai(Vf, Ef, Em, direction='transverse'):
    # xi values from http://www.mse.mtu.edu/~drjohn/my4150/ht/ht.html
    # if mode == 'fibre':
    if direction == 'axial':
        return calc_axial(Vf, Ef, Em)
    elif direction == 'transverse':
        xi = 2 + 40 * Vf**10

    Vm = (1 - Vf)
    return Em * (Ef + xi * (Vf * Ef + Vm * Em)) / (Vf * Em + Vm * Ef + xi * Em)

def _fibre_composite_E(Ef, Em, HT=True):
    Vf = np.linspace(0,1)

    fig, ax = plt.subplots(1, 1)

    ax.plot(Vf, calc_axial(Vf, Ef, Em), label='Axial')
    ax.plot(Vf, calc_transverse(Vf, Ef, Em), label='Transverse')
    if HT:
        ax.plot(Vf, calc_halpintsai(Vf, Ef, Em), label='Halpin-Tsai')
    
    ydiff = abs(Ef - Em)
    ax.axhline(Ef, ls='dashed', color=(0,0,0,0.5))
    ax.text(0, Ef - 0.02 * ydiff, '$E_f$', va='top', fontsize=14)
    ax.axhline(Em, ls='dashed', color=(0,0,0,0.5))
    ax.text(1, Em + 0.02 * ydiff, '$E_m$', ha='right', fontsize=14)
    
    ax.legend(title='Direction')
    ax.set_xlabel('Fractional Volume of Fibre ($V_f$)')
    ax.set_ylabel("Young's Modulus ($E$)")

def fibre_composite_E():
    # interact(_fibre_composite, Ef=widgets.FloatSlider(120, min=1, max=200, description='Fibre E'), Em=widgets.FloatSlider(10, min=1, max=200, description='Matrix E'))
    interact(_fibre_composite_E, 
             Ef=widgets.FloatText(120, description='$E_f$'),
             Em=widgets.FloatText(10, description='$E_m$'),
             HT=widgets.Checkbox(value=False, description='Halpin-Tsai'))

# Fracture strength fibre composite

def calc_fracturestrength(Vf, Ef, Em, epsilonb_fibre, epsilonb_matrix):
    Ea = calc_axial(Vf, Ef, Em)
    sigma_bf = Ea * epsilonb_fibre
    # sigma_bf = Vf 

    sigma_bm = epsilonb_matrix * Em * (1 - Vf)

    out = sigma_bf
    out[sigma_bm > sigma_bf] = sigma_bm[sigma_bm > sigma_bf]

    return out

def _fibre_composite_frac(Ef, Em, epsilonb_fibre, epsilonb_matrix):
    Vf = np.linspace(0,1)

    fig, ax = plt.subplots(1, 1)

    ax.plot(Vf, calc_fracturestrength(Vf, Ef, Em, epsilonb_fibre, epsilonb_matrix))
    ax.set_ylabel('Axial Fracture Strength ($\sigma_b$)')
    ax.set_xlabel('Fractional Volume of Fibre ($V_f$)')

    ydiff = np.diff(ax.get_ylim())
    ax.axhline(Ef * epsilonb_fibre, ls='dashed', color=(0,0,0,0.5))
    ax.text(0, Ef * epsilonb_fibre - 0.02 * ydiff, '$\sigma_{fb}$', va='top', fontsize=14)
    ax.axhline(Em * epsilonb_fibre, ls='dashed', color=(0,0,0,0.5))
    ax.text(1, Em * epsilonb_fibre + 0.02 * ydiff, '$\sigma_{fm}$', ha='right', fontsize=14)

def fibre_composite_fracture():
    interact(_fibre_composite_frac, 
             Ef=widgets.FloatText(120, description='$E_f$'),
             Em=widgets.FloatText(10, description='$E_m$'),
             epsilonb_fibre=widgets.FloatText(1.2, description='$\epsilon_{fb}$'),
             epsilonb_matrix=widgets.FloatText(1.2, description='$\epsilon_{fm}$'))



## Moment as a function of distance

def _moment_distance(d):
    x = np.linspace(0, 1, 500)
    L = 1
    P = 1
    M = (1 - d) * P * x
    M[x > d] = d * P * (L - x[x > d])

    plt.plot(x, M)
    plt.ylim(0, 0.26)
    plt.axvline(d, ls='dashed', c=(0,0,0,0.5))
    plt.axhline(max(M), ls='dashed', c=(0,0,0,0.5))
    plt.ylabel('M / P')

    plt.xticks([0, d, 1], ['S', 'P', 'S'])

def moment_distance():
    interact(_moment_distance,
             d=widgets.FloatSlider(0.2, min=0, max=1, step=0.01, description='P Location'))

# Flexural Stress

def _flexural_stress(L, P):
    rx = np.linspace(0.5, 5)
    M = L * P / 4
    I = rx**4 * np.pi / 4

    sigma = M * rx / I

    plt.plot(rx, sigma)
    plt.xlabel('Radius')
    plt.ylabel('$\sigma$')

def flexural_stress():
    interact(
        _flexural_stress,
        L=widgets.FloatSlider(1, min=0.1, max=10),
        P=widgets.FloatSlider(1, min=0.1, max=10),
    )

# Cylinder fracture

def calc_bend(R, L):
    theta = np.arcsin(0.5 * L / R)
    H = np.cos(theta) * R

    theta_range = np.linspace(-theta - np.pi/2, theta - np.pi/2)

    x = R * np.cos(theta_range)
    y = R * np.sin(theta_range)
    
    return x, y + H

def _cylinder_mechanics(E, P, r1, r2, L, frac_stress=np.inf):
    try:
        fig = plt.figure(figsize=[10, 5.5])

        gs = GridSpec(3, 2)

        ax0 = fig.add_subplot(gs[0, :])
        ax1 = fig.add_subplot(gs[1:, 0])
        ax2 = fig.add_subplot(gs[1:, 1])

        axs = [ax0, ax1, ax2]
        
        # handle units
        # lengths to m
        r1 /= 1e6
        r2 /= 1e6
        L /= 1e6  # m

        # GPa to Pa
        E *= 1e9
        frac_stress *= 1e9

        ngrid = 500
        x = y = np.linspace(-r2, r2, ngrid)
        X, Y = np.meshgrid(x, y)

        ir = np.sqrt(X**2 + Y**2)
        mask = (ir > r1) & (ir <= r2)

        I = (r2**4 - r1**4) * np.pi / 4
        M = P * L / 4

        stress = M * Y / I
        strain = stress / E

        stress_max = np.nanmax(abs(stress))

        if np.isfinite(frac_stress):
            stresslim = (-min(stress_max, frac_stress), min(stress_max, frac_stress))
        else:
            stresslim = (np.nanmin(stress), np.nanmax(stress))

        strain[~mask] = np.nan
        stress[~mask] = np.nan

        cm1 = plt.cm.RdBu
        cm2 = plt.cm.PRGn

        cm2.set_over('r')
        cm2.set_under('r')

        c1 = ax1.pcolormesh(1e6 * X, 1e6 * Y, strain, cmap=cm1)
        c2 = ax2.pcolormesh(1e6 * X, 1e6 * Y, stress * 1e-9, cmap=cm2, vmin=stresslim[0] * 1e-9, vmax=stresslim[1] * 1e-9)
        
        fig.colorbar(c1, ax=ax1)
        fig.colorbar(c2, ax=ax2)
            
        ax1.set_title('Strain')
        ax2.set_title('Stress (GPa)')
        
        for ax in axs[1:]:
            ax.set_aspect(1)
            ax.set_xlim(-r2 * 1e6 * 1.1, r2 * 1e6 * 1.1)
            ax.set_ylim(-r2 * 1e6 * 1.1, r2 * 1e6 * 1.1)
            ax.axhline(0, ls='dashed', color=(0,0,0,0.5))
            
            ecircle2 = plt.Circle((0,0), 1e6 * r2, edgecolor='k', facecolor=(0,0,0,0), lw=3)
            ecircle1 = plt.Circle((0,0), 1e6 * r1, edgecolor='k', facecolor=(0,0,0,0), lw=3)
            ax.add_artist(ecircle2)
            ax.add_artist(ecircle1)
        
        # bend graphic
        D = P * L**3 / (48 * E * I)
        try:
            R = (L**2 + 4 * D**2) / (8 * D)
        except ZeroDivisionError:
            R = 10000
        if stress_max > frac_stress:
            c = 'r'
            lw = 3
        else:
            c = 'k'
            lw = 2

        x, y = calc_bend(R, L)
        
        ax0.plot(1e6 * x, 1e6 * y, color=c, lw=lw)
        ax0.set_xlim(-1e6 * L/2, 1e6 * L/2)

        if min(y) < -0.05e-6:
            ymin = min(y * 1e6) - 0.05
        else:
            ymin = -0.05
        ax0.set_ylim(ymin, 0.02)
        ax0.axhline(0, c=(0,0,0,0.25))
        ax0.axhline(-1e6 * D, ls='dashed', c=(0,0,0,0.25))
        ax0.axvline(0, ls='dashed', c=(0,0,0,0.5))
        ax0.set_xticks([-1e6 * L/2, 0, 1e6 * L/2])
        ax0.set_xticklabels([0, int(1e6 * L / 2), int(L * 1e6)])
        ax0.set_ylabel('Deflection ($\mu m$)')
        ax0.set_xlabel('Length ($\mu m$)')

        ax0.text(.01, .05, f'Deflection: {1e6 * D:.2e}' + ' $\mu m$, ' + f'{100 * D/L:.2f} %', transform=ax0.transAxes)

        fig.tight_layout()
    
    except:
        print('Nope. Looks like one of the values is wrong.')

def cylinder_mechanics():
    interact(
        _cylinder_mechanics,
        E=widgets.FloatText(30, description="$E$ (GPa)", step=0.1),
        P=widgets.FloatText(5, description="$P$ (N)", step=0.05),
        r1=widgets.FloatText(0, description="$r_i$ ($\mu m$)", step=0.1),
        r2=widgets.FloatText(200, description="$r_o$ ($\mu m$)", step=0.1),
        L=widgets.FloatText(500, description="$L$ ($\mu m$)", step=5),
        frac_stress=widgets.FloatText(0.1, description="$\sigma_f$ (GPa)", step=0.1)
        )

# def _cylinder_mechanics(E, R, r1, r2, fig, axs, frac_stress=np.inf):
#     ax0, ax1, ax2 = axs
    
#     ngrid = 500
#     x = y = np.linspace(-r2, r2, ngrid)
#     X, Y = np.meshgrid(x, y)

#     ir = np.sqrt(X**2 + Y**2)
#     mask = (ir > r1) & (ir <= r2)

#     strain = - Y / R
#     stress = strain * E
#     stress_max = np.nanmax(abs(stress))
#     if np.isfinite(frac_stress):
#         stresslim = (-min(stress_max, frac_stress), min(stress_max, frac_stress))
#     else:
#         stresslim = (np.nanmin(stress), np.nanmax(stress))

#     strain[~mask] = np.nan
#     stress[~mask] = np.nan

#     cm1 = plt.cm.RdBu
#     cm2 = plt.cm.PRGn

#     cm2.set_over('r')
#     cm2.set_under('r')

#     c1 = ax1.pcolormesh(X, Y, strain, cmap=cm1)
#     c2 = ax2.pcolormesh(X, Y, stress, cmap=cm2, vmin=stresslim[0], vmax=stresslim[1])
    
#     fig.colorbar(c1, ax=ax1)
#     fig.colorbar(c2, ax=ax2)
        
#     ax1.set_title('Strain')
#     ax2.set_title('Stress')
    
#     for ax in axs[1:]:
#         ax.set_aspect(1)
#         ax.set_xlim(-r2 * 1.1, r2 * 1.1)
#         ax.set_ylim(-r2 * 1.1, r2 * 1.1)
#         ax.axhline(0, ls='dashed', color=(0,0,0,0.5))
        
#         ecircle2 = plt.Circle((0,0), r2, edgecolor='k', facecolor=(0,0,0,0), lw=3)
#         ecircle1 = plt.Circle((0,0), r1, edgecolor='k', facecolor=(0,0,0,0), lw=3)
#         ax.add_artist(ecircle2)
#         ax.add_artist(ecircle1)
    
#     # bend graphic
    
#     Rmax = -E * r2 / frac_stress


#     if stress_max > frac_stress:
#         c = 'r'
#         lw = 3
#     else:
#         c = 'k'
#         lw = 2

#     x, y = calc_bend(R)
    
#     ax0.plot(x, y, color=c, lw=lw)
#     ax0.set_xlim(-.5, .5)

#     if min(y) < -0.1:
#         ymin = min(y) - 0.05
#     else:
#         ymin = -0.15
#     ax0.set_ylim(ymin, 0.05)
#     ax0.axhline(0, ls='dashed', c=(0,0,0,0.5))
#     ax0.axvline(0, ls='dashed', c=(0,0,0,0.5))
#     ax0.set_xticks([])

# def _plot_cylinder_mechanics(E, R, r1, r2, frac_stress=np.inf):
#     fig = plt.figure(figsize=[10, 6])

#     gs = GridSpec(3, 2)

#     ax0 = fig.add_subplot(gs[0, :])
#     ax1 = fig.add_subplot(gs[1:, 0])
#     ax2 = fig.add_subplot(gs[1:, 1])

#     axs = [ax0, ax1, ax2]

#     _cylinder_mechanics(E, R, r1, r2, fig, axs, frac_stress)

#     fig.tight_layout()

# def cylinder_mechanics():
#     interact(
#         _plot_cylinder_mechanics,
#         E=widgets.FloatText(1.5, description="$E$", step=0.1),
#         R=widgets.FloatText(10, description="$R$", step=0.1),
#         r1=widgets.FloatText(0, description="$r_i$", step=0.1),
#         r2=widgets.FloatText(5, description="$r_o$", step=0.1),
#         frac_stress=widgets.FloatText(10, description="$\sigma_f$", step=0.1)
#         )