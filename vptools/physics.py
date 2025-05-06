#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Contains various physics models useful for analyzing experimental data.
"""
#### basic imports
import numpy as np
from scipy.optimize import root

def rigid_rotor(r, B0, n0, rs, rw=0.8):
    """
    Computes the rigid rotor model.
    
    Arguments:
        r : radial coordinate, 1D array-like [m]
        B0 : external field, float [T]
        n0 : maximum density, float [m^-3]
        rs : seperatrix radius, float [m]
        rw : conducting wall radius, float [m]
        
    Returns:
        Bz : axial magnetic flux density, 1D array like len(r), [T]
        psi : magnetic flux, 1D array like len(r), [Wb]
        ne : electron density, 1D array like len(r), [m^-3]
    """
    
    k = root(_rrk, 0.5, args=(rs,rw)).x
    u = 2 * r**2 / rs**2 - 1
    Bz = B0 * np.tanh(k * u)
    psi = (rs**2 * B0 / (4 * k)) * np.log(np.cosh(k * u) / np.cosh(k))
    ne = n0 / np.cosh(k * u)**2
    
    return Bz, psi, ne

def solovev(r, z, B0, rs, zs, rw=0.8):
    """
    Computes the solovev flux model
    
    Arguments:
        r : radial coordinate, 1D array-like [m]
        z : axial coordinate, 1D array-like [m]
        B0 : external field, float [T]
        rs : seperatrix radius, float [m]
        zs : seperatrix half-length, float [m]
        rw : conducting wall radius, float [m]
        
    Returns:
        psi : magnetic flux, 2D array like [len(z),len(r)], [Wb]
        Bz : axial magnetic flux density, 2D array like [len(z),len(r)], [T]
        Br : radial magnetic flux density, 2D array like [len(z),len(r)], [T]
    """
    
    R, Z = np.meshgrid(r, z)
    psi0 = rs**2 * B0 / (2 * (1 - 2 * rw**2 / rs**2))
    psi = psi0 * (R / rs)**2 * (1 - (R / rs)**2 - (Z / zs)**2)
    Bz = (psi0 / rs**2) * (1 - 2 * (R / rs)**2 - (Z / zs)**2)
    Br = (2 * psi0 / (rs * zs)) * (R / rs) * (Z / zs)
    
    return psi, Bz, Br

def calc_sigma_cx(E):
    """
    Computes the Janev CX cross section.
    
    Arguments:
        E : eV
    
    Returns:
        sigma_cx : cross section, in m^2
    """
    
    return 0.6937e-14 * (1 - 0.155 * np.log10(E))**2 / (1 + 0.1112e-14 * E**(3.3)) * (1/100)**2

def _rrk(k, rs, rw=0.8):
    """
    Helper function for rigid rotor.
    """
    return k * (1 - 0.5 * rs**2 / rw**2) - np.tanh(k)


if __name__=='__main__':
    
    #### local imports
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    mpl.style.use('../basic.mplstyle')
    
    dx = 0.01
    r = np.linspace(0,0.8,int(0.8/dx)+1)
    z = np.linspace(-3,3,int(6/dx)+1)
    
    B0, n0, rs, zs = 0.1, 1e19, 0.4, 1.6
    Bz_rr, psi_rr, ne_rr = rigid_rotor(r, B0, n0, rs)
    psi_sv, Bz_sv, Br_sv = solovev(r, z, B0, rs, zs)
    
    fig, ax = plt.subplots(1,1,figsize=(10.0,7.0))
    ax.plot(r, Bz_rr * 10, label=r'B$_z$ [kG]')
    ax.plot(r, psi_rr * 1e2, label=r'$\psi$ [$10^{-1}$ mWb]')
    ax.plot(r, ne_rr / 1e19, label=r'n$_e$ [$10^{-19}$ m$^{-3}$]')
    ax.set_xlabel('r [m]')
    ax.set_ylabel('Arb.')
    ax.set_title('Rigid Rotor Model')
    ax.legend(fontsize=24)
    
    fig, ax = plt.subplots(1,1,figsize=(10.0,6.0))
    p = ax.pcolor(z,r,Bz_sv.T,cmap='seismic', vmin=-B0, vmax=B0)
    ax.contour(z,r,psi_sv.T,[0],colors='k')
    ax.contour(z,r,Bz_sv.T,[0],colors='k', linestyles='--')
    ax.set_xlabel('z [m]')
    ax.set_ylabel('r [m]')
    ax.set_title('Solovev Model')
    c = plt.colorbar(p)
    c.set_label(r'B$_z$ [T]')
    handles = [Line2D([],[],lw=3,color='k',label='Seperatrix'),
               Line2D([],[],lw=3,color='k',ls='--',label='Magnetic Null')]
    ax.legend(handles=handles, loc='upper left', fontsize=24)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    