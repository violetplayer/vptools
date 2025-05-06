#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

def make_colorbar(p, ax, ticks=None, colorbar_scaling=[0.04, 0.075]):
    """
    Utility which constructs a colorbar with default scaling values.
    """
    
    if ticks is None:
        c = plt.colorbar(p, ax=ax, fraction=colorbar_scaling[0], pad=colorbar_scaling[1])
    else:
        c = plt.colorbar(p, ax=ax, ticks=ticks, fraction=colorbar_scaling[0], pad=colorbar_scaling[1])
        
    return c

def plot_discrete_mesh(x, y, data, ax=None, colormap='jet', nvals=None):
    """
    Utility which plots a discrete colormap, with centered colorbar labels
    
    Inputs:
        x: array-like, x coordinate
        y: array-like, y coordinate
        data: array like, dimensions [x, y]
        labels: [title, x, y, colorbar] labels
        colorbar_scaling: scaling factors for colorbar, [vertical, horizontal]
        ax: pyplot axis object. If none, opens a new figure.
        colormap: colormap to use in plot
        
    Returns:
        p: plot object
        c: colorbar object
    """
    
    if nvals is None:
        nvals = (np.min(data), np.max(data))
    
    #get discrete colormap
    cmap = plt.get_cmap(colormap, nvals[1]-nvals[0]+1) #np.max(data)-np.min(data)+1)
    
    output_fig = False
    if ax is None:
        output_fig = True
        fig, ax = plt.subplots(1,1,figsize=(12.0, 9.0))
    
    p = ax.pcolormesh(x, y, data, cmap=cmap, vmin=nvals[0]-0.5, vmax=nvals[1]+0.5)
                      #vmin=np.min(data)-.5, vmax=np.max(data)+.5)  
    ticks = np.arange(nvals[0],nvals[1]+1) #np.min(data), np.max(data)+1)
    c =  make_colorbar(p, ax, ticks=ticks)
    
    if output_fig:
        return fig, ax, c, p
    else:
        return c, p

def plot_colorchange_line(x, y, c, axis=None, colorbar=False, cmap='turbo', limits=None):
    
    if axis is None:
        fig, ax = plt.subplots(1,1,figsize=(12.0,9.0))
    else:
        ax = axis
    
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)    
    norm = plt.Normalize(c.min(), c.max())
    lc = LineCollection(segments, cmap=cmap, norm=norm)    
    lc.set_array(c)  
    line = ax.add_collection(lc)
    
    if limits is None:
        dx = np.diff(x)[0]
        dy = 0.05 * np.nanmax(y)
        ax.set_xlim([np.nanmin(x) - 10 * dx, np.nanmax(x) + 10 * dx])
        ax.set_ylim([np.nanmin(y) - 10 * dy,np.nanmax(y) + 10 * dy])
    else:
        ax.set_xlim([limits[0], limits[1]])
        ax.set_ylim([limits[0], limits[1]])
    
    if colorbar:
        c = make_colorbar(line, ax=ax)
    
    if axis is None:
        if colorbar:
            return fig, ax, c, line
        else:
            return fig, ax, line
    else:
        if colorbar:
            return c, line
        else:
            return line
    

def plot_multicolor_lines(x, y, n, axis=None, legend=False):
    """
    
    Utility which plots a series of lines with sequential colors using the turbo colormap.
    
    Inputs:
        x: array-like of x coordinate for lines. Preferred format list of array-likes, or array of dimension [line, data]
        y: array-like of y coordinate for lines. Preferred format list of array-likes, or array of dimension [line, data]
        n: number of lines
        ax: pyplot axis object. If none, opens a new figure
        legend: if true, returns a legend handles list which can be used to populate a legend on the plot
    
    Returns:
        p: plot object
        
    TODO:
        implement legend output
    
    """
    cmap = mpl.cm.turbo(np.linspace(0,1,n))
    
    if axis is None:
        fig, ax = plt.subplots(1,1,figsize=(12.0,9.0))
    
    p_out = []
    for i in range(n):
        p_tmp = ax.plot(x[i], y[i], lw=3, color=cmap[i])
        p_out.append(p_tmp)
    
    if legend:
        
        raise Exception('Not implemented')       # TODO: implement this  
                
    else:
        if axis is None:
            return fig, ax, p_out
        else:
            return p_out

if __name__=='__main__':
    
    mpl.style.use(r'..\basic.mplstyle')
    
    #### show off multicolored line plot
    x = np.linspace(0, 10, 1000)
    y = np.cos(np.pi * x)
    fig, ax, c, line = plot_colorchange_line(x, y, colorbar=True)
    ax.set_title(r'cos($\pi$x)')
    
    #### how about a 2D plot
    x, y = np.linspace(-2, 2, 100), np.linspace(-2, 2, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.ceil(X**2 + Y**2)
    fig, ax, c, p = plot_discrete_mesh(x, y, Z)
    ax.set_title('ceil(x^2 + y^2)')
    
    #### what if we had many plots
    x = np.linspace(0, 10, 1000)
    y = [(i+1) + np.cos(np.pi * x / (i + 1)) for i in range(10)]
    x_out = [x for i in range(10)]
    fig, ax, p = plot_multicolor_lines(x_out, y, len(y))
    ax.set_title(r'n + cos($\pi$x/n)')
    
    
    
    
    

