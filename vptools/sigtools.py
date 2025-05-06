#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from scipy.fft import fftshift, fftfreq, fft
from scipy.signal import stft, savgol_filter, find_peaks

def basic_fft(x, y):
    """
    Wraps a very basic 1D fft.
    Data should be uniformly spaced with only finite values.
    To retrieve spectral power, compute abs(yf) * (2 / len(x)).
    
    Arguments:
        x : axis
        y : data values
        
    Returns:
        xf : frequency values (upper-half; units of 1/dim(x))
        yf : spectral density (upper-half; complex, un-normalized)
    """
    
    N, dx = len(x), np.mean(np.diff(x))
    xf = fftshift(fftfreq(N, dx))[N//2:]
    yf = fftshift(fft(y))[N//2:]
    
    return xf, yf

def basic_stft(x, y, nperseg=4096):
    """
    Wraps a very basic 2D stft.
    Data should be uniformly spaced with only finite values.
    To retrieve spectral power, compute abs(Zxx) * (2 / len(x)).
    
    Arguments:
        x : axis
        y : data values
        nperseg : width of stft window
        
    Returns:
        F : frequency values, units of 1/dim(x)
        T : time values, units of (x)
        T : 2D spectral density (complex, un-normalized)
    """
    
    fs = 1/np.mean(np.diff(x))
    F, T, Zxx = stft(y, fs, nperseg=nperseg)
    
    return F, T + x[0], Zxx
    

def get_flattop_bounds(t, y, thresh_rel=0.85, min_sep=0.1, min_length=1.0, thresh_abs=None):
    
    """
    Gets the flat-top boundaries of signal y, usually rx.
    Designed to handle mid-shot crashes + recovery.
    
    Arguments:
        t : time axis
        y : signal values
        thresh_rel : relative threshold, flat-top where y / max(y) >= thresh_rel
        min_sep : minimum seperation between flat-tops [ms]
        min_length : minimum length of flat-top [ms]
        thresh_abs : absolute threshold, flat-top where y >= thresh_abs
        
    Returns:
        t_bounds : 2D array. Dim0: flat-top index, Dim1: (t_0, t_1)
    
    """
    
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

def identify_beam_pulses(t, I, t_bounds, t_window=0.5, deriv_threshold=0.5, min_dist=0.5):
    """
    Returns the indices of rising and falling edges of the beam pulses.
    
    Args:
        t : time axis
        I : NB current
        t_bounds : 2-element arraylike, upper and lower bound of time window to consider
    
    Return:
        pks_out : indices of beam edges on time axis
    """
        
    #### get time window
    dt = t[1] - t[0]
    Nt = int(t_window / dt)
    if Nt%2 == 0:
        Nt += 1
    Nd = int(min_dist / dt)
    
    #### find peaks in first derivative
    dI_dt  = savgol_filter(I, window_length=Nt, polyorder=2, deriv=1)
    dI_mag = np.abs(dI_dt)
    dI_max = np.nanmax(dI_mag)
    pks, _ = find_peaks(dI_mag, height=deriv_threshold * dI_max, distance=Nd)
    
    #### throw out peaks out-of-bounds
    idx_t   = [np.argmin(np.abs(t - t_bounds[0])), np.argmin(np.abs(t - t_bounds[1]))]
    pks_out = pks[np.where((pks > idx_t[0]) & (pks < idx_t[1]))[0]]
    
    return pks_out
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    