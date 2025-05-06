#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy import signal, interpolate

def smooth(x, N):
    return signal.savgol_filter(

def identify_events(y, t,
                    mode, t_smooth,
                    **peak_params):
    """
    Given a 1D signal y(t), detect events. Mode determines functionality.
    
    Inputs:
        y        : 1d array-like of signal
        t        : time axis of signal [ms]
        t_smooth : width of smoothing window
        mode     : how to identify events. Options detailed below.
        **peak_params : parameters to pass to the find_peaks function. See [scipy link] for details.
                        Not used for [zero_sig, extr_sig, max_abs_grad]
    
    Returns:
        event_indices : indices of events in y
    
    mode options:
        max_sig      : maxima in the signal
        min_sig      : minima in the signal
        zero_sig     : zeros in the signal
        extr_sig     : all extrema in the signal
        max_abs_grad : maxima in the absolute value of the gradient (all inflection points)
        max_grad     : maxima in the gradient (upward slope)
        min_grad     : minima in the gradient (downward slope)
    """
    
    #### first, smooth the signal & get the gradient
    dt       = t[1] - t[0]
    Nt       = int(t_smooth / dt)
    y_smooth = signal.savgol_filter(y, Nt, 0)
    dy_dt    = signal.savgol_filter(y, Nt, 1, deriv=1)
    
    #### depending on chosen mode, how do we identify extrema?
    if mode == 'max_sig':
        y_id = y_smooth
    elif mode == 'min_sig':
        y_id = -1 * y_smooth
    elif mode == 'zero_sig':
        y_id = y_smooth
    elif mode == 'extr_sig':
        y_id = dy_dt
    elif mode == 'max_abs_grad':
        y_id = np.gradient(smooth(dy_dt, Nt), t)
    elif mode == 'max_grad':
        y_id = dy_dt
    elif mode == 'min_grad':
        y_id = -1 * dy_dt
    else:
        raise Exception('Input a valid mode')
    
    #### zero-finding method
    if mode in ['zero_sig', 'extr_sig', 'max_abs_grad']:
        event_indices = find_zeros(y_id)
        
    #### peak-finding method
    else:
        #### some default options
        if len(peak_params.keys()) == 0:
            peak_params['prominence'] = 0.25 * np.nanmax(y_id)
            peak_params['width']      = 2 * Nt
        
        event_indices, _ = signal.find_peaks(y_id, **peak_params)
    
    return event_indices

def collect_events(y, t, event_indices, t_window, 
                   interp=False, t_new=None):
    
    """
    Given a 1D signal y(t), and a set of event indices, collects the events into a nice format.
    
    Inputs:
        y             : signal, 1d array-like
        t             : time axis of signal [ms]
        event_indices : array-like of indices of events - assumes this is the "middle" of the event
        t_window      : width of event window in milliseconds
        interp        : if True, interpolates events onto new time axis, provided in **relative** time.
                        Useful for inconsistent sampling rates
        t_new         : new time. if interpolate==True, required.
        
    Returns:
        t_rel   : relative time of width t_window [ms]
        t_event : peak times of events, len(event_indices) [ms]
        y_event : array of y-values in event windows, dimension [N_events, len(t_rel)]
        
    """
    if interpolate and t_new is None:
        raise Exception('Provide a new time axis for interpolation.')
        
    #### collect events
    y_event = []
    t_event = []
    t_rel   = []
    for event in event_indices:
        
        idx_tmp = np.where((t >= 0.98 * (t[event] - t_window / 2)) & (t <= 1.02 * (t[event] + t_window / 2)))[0]
        t_event.append(t[event])
        y_event.append(y[idx_tmp])
        t_rel.append(t[idx_tmp] - t[event])
                
    if interp:
        y_event_interp = []
        for i in range(len(t_event)):
            y_fun = interpolate.interp1d(t_rel[i], y_event[i], bounds_error=False, fill_value=np.nan)
            y_event_interp.append(y_fun(t_new))
            
        return t_new, np.array(t_event), np.array(y_event_interp)
    
    else:
        return t_rel, np.array(t_event), np.array(y_event)


def event_param_avg(p, t, event_times, t_window):
    
    """
    Given a 1D parameter signal p(t) and a set of event **times** (not indices), collect the 
    average of the parameter about each event over width of t_window
    
    Inputs:
        p           : parameter signal, 1d-array like
        t           : time axis of signal [ms]
        event_times : array-like of event times [ms] - assumes this is the 'middle' of the event
        t_window    : width of event window [ms]
        
    Returns:
        p_event : array of averaged parameter in event windows, len(event_times)
    
    """
    
    #### collect event averages
    p_event = []
    for time in event_times:
        
        idx_tmp = [np.argmin(np.abs(t - (time - t_window / 2))),
                   np.argmin(np.abs(t - (time + t_window / 2))) + 1]
        p_event.append(np.nanmean(p[idx_tmp[0]:idx_tmp[1]]))
    
    return np.array(p_event)
    
def find_zeros(y):
    """
    Given an array y, find the zero crossings of that array
    
    Inputs:
        y : 1d array like
        
    Returns:
        indices : indices of zero crossings (nearest point)
    
    """
    
    #### identify signs of array values
    ysign = np.sign(y)
    
    #### fixes the convention of 0 for where an array exactly equals 0
    sz = ysign == 0
    while sz.any():
        ysign[sz] = np.roll(ysign, 1)[sz]
        sz = ysign == 0
    
    #### identifies sign change indices
    signchange    = ((np.roll(ysign, 1) - ysign) != 0).astype(int)
    signchange[0] = 0
    indices       = np.where(signchange)[0]
    
    return indices
    
    
    
    
    
