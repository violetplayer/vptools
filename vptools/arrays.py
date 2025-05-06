#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.interpolate import griddata

def get_bin_edges(x):
    """
    General utility for grabbing bin edges, given bin centers.
    
    Useful for histograms and correctly formatting pcolormesh maps.
    
    Arguments:
        x : array-like, 1D
        
    Returns:
        edg : 1D numpy array, len(x) + 1
    """
    
    dx = x[1] - x[0]
    edg = np.zeros(len(x) + 1)
    edg[0] = x[0] - dx / 2
    edg[1:] = x + dx /2
    
    return edg

def map_over_nans(x, y, A, new_coords=None, method='linear', fill_value=np.nan):
    
    """
    Interpolates across NaN values in a 2D array.
    """
    array  = np.ma.masked_invalid(A) ### create masked array
    xx, yy = np.meshgrid(x, y, indexing='ij')
    x1, y1 = xx[~array.mask], yy[~array.mask]
    array_real = array[~array.mask]
    
    if new_coords is not None:
        xx1, yy1 = np.meshgrid(new_coords[0], new_coords[1])
        A_map = griddata((x1, y1), array_real.ravel(), (xx1, yy1), method=method, fill_value=fill_value)
    else:
        A_map = griddata((x1, y1), array_real.ravel(), (xx, yy), method=method, fill_value=fill_value)
    
    return A_map

def guess_roots(x, y, adjust=True, return_coords=False):
    """
    Attempts to guess root positions. Returns empty array if no roots exist.
    Error will be <= diff(x). Excellent for generating starting points for Newton's method.
    
    Arguments:
        x : x coordinate, array-like 1D
        y : y coordinate, array-like 1D
        
    Returns:
        Roots : 1D array-like containing all root guesses
    
    """
    #### finds sign changes
    y_sign = np.sign(y)
    delta = np.roll(y_sign, 1) - y_sign
    loc = delta != 0
    loc[0] = 0 # not periodic
    arr = np.where(loc)[0] # needed for the next part
    if adjust: # adjusts to keep indices *inside* true root
        for i, _arr in enumerate(arr):
            if y[_arr] <= 0:
                if delta[_arr] < 0:
                    arr[i] += 1
                else:
                    arr[i] -= 1
    roots = x[arr]
    if return_coords:
        return roots, arr
    else:
        return roots
    
    
    
def get_flattop_bounds(t, y, thresh_rel=0.85, min_sep=0.1, min_length=1.0, thresh_abs=None):
    
    dt = t[1] - t[0]
    N_sep = int(min_sep / dt)
    N_len = int(min_length / dt)
    
    if thresh_abs is not None:
        idx_ftop = np.where(y >= thresh_abs)[0]
    else:
        idx_ftop = np.where(y >= thresh_rel * np.nanmax(y))[0]
        
    if len(idx_ftop) == 0:
        return None
        
    else:
        diff = np.diff(idx_ftop)
        idx_d = np.where(diff >= N_sep)[0]
        if len(idx_d) > 0:
            t_bounds = np.zeros([len(idx_d)+1,2])
            for i in range(len(idx_d)+1):
                if i == 0:
                    idx_i = idx_ftop[0:idx_d[i]-1]
                elif i > 0 and i < len(idx_d):
                    idx_i = idx_ftop[idx_d[i-1]+1:idx_d[i]-1]
                else:
                    idx_i = idx_ftop[idx_d[-1]+1:]

                if len(idx_i) >= N_len:
                    t_bounds[i,:] = [t[idx_i][0], t[idx_i][-1]]
                else:
                    t_bounds[i,:] = [np.nan, np.nan]
        else:
            t_bounds = np.array([t[idx_ftop][0], t[idx_ftop][-1]]).reshape([1, 2])

        idx = np.where(~np.isnan(t_bounds[:,0]))[0]
        t_bounds = t_bounds[idx,:]

        return t_bounds