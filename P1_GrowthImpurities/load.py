import pandas as pd

def Teng_LC():
    """
    Load Critical Step Length (LC) data from Teng et al. (1998).

    Columns
    -------
    1_sigma : float
        1 / saturation (unitless)
    LC : float
        Critical step length (nm)
    type : string
        The type of step (acute of obtuse)
    """
    return pd.read_csv('data/Teng1998_LC.csv', comment='#')

def Teng_v():
    """
    Load Step Velocity (v) data from Teng et al. (1998).

    Columns
    -------
    type : string
        The type of step (acute of obtuse)
    sigma : float
        Saturation (unitless)
    L : float
        The length of the step (nm)
    v : float
        Speed of step propagation normal to terrace edge direction (nm/s)
    """
    return pd.read_csv('data/Teng1998_v.csv', comment='#')

def Teng_W():
    """
    Load Terrace Width (W) data from Teng et al. (1998).

    Columns
    -------
    type : string
        The type of step (acute of obtuse)
    sigma : float
        Saturation (unitless)
    W : float
        Terrace width normal to growth direction (nm)
    """
    return pd.read_csv('data/Teng1998_W.csv', comment='#')

def Davis_v():
    """
    Load Step Velocity (v) data from Davis et al. (2000).

    Columns
    -------
    aCa : float
        The activity of Ca in solution (mol/kg).
    v : float
        Speed of step propagation normal to terrace edge direction (nm/s)
    treatment : string
        Name of experiment
    xMg : float
        Mg activity as multiple of 4.96875e-7 mol/kg.
    """
    return pd.read_csv('data/Davis2000.csv', comment='#')