#!/usr/bin/env python

import numpy as np
from matplotlib import pyplot
from scipy.interpolate import UnivariateSpline
def scatter_plot_from_dict(ax,dict1,dict2,keys='auto'):
    """Function designed to take two dictionaries and convert them into a
    plot, where both dictionaries have the form dict[key] = value. This
    will make a plot of the values at all points where the dictionaries
    have common keys and add them to the axis object ax."""

    if keys == 'auto':
        keys = []
        keys1 = dict1.keys()
        keys2 = dict2.keys()
        testkeys = list(set(keys1).intersection(set(keys2)))
        for testkey in testkeys:
            if not np.isnan(dict1[testkey]): 
                if not np.isnan(dict2[testkey]):
                    keys.append(testkey)
    x,y = np.ones(len(keys)),np.ones(len(keys))
    for index,key in enumerate(keys):
        x[index] = dict1[key]
        y[index] = dict2[key]
    ax.plot(x,y,'.')
    for key in keys:
        ax.text(dict1[key],dict2[key],key)

class PlotStates:
    """Functions to update existing plots with energy levels and connectors.
    Initialize as

    PlotStates(ax,G,halfwidth,textwidth,fontsize,color,textposition,
               text_vspace)

    ax is the axis object to plot to
    G is the free energy dictionary, G[state] = value
    halfwidth is half of the width the plot the horizontal bars.
    textwidth is the amount of space to leave for the text labels.
    fontsize is the fontsize for the text labels
    color is a matplotlib compatiable color specifier
    textposition can be 'inline' for the text to be part of the bar,
        or 'above' to be above the bar, if 'above' is specified, then
        text_vspace may be added to move the text farther above the line
    text_vspace is a float, in units the same as the y axis, of how much
        to move the text above the line
    """

    def __init__(self, ax, G, halfwidth, textwidth, fontsize, color,
                textposition='inline', text_vspace=0.,lw=2.):
        self.ax = ax
        self._G = G
        self._halfwidth = halfwidth
        self._textwidth = textwidth
        self._fontsize = fontsize
        self._color = color
        self._textposition = textposition
        self._text_vspace = text_vspace
        self._lw=lw

    def plotstate(self, state, zone, label=None):
        """Plots a horizontal line at the specified energy state. Pulls 
        the energy from the G dictionary. State should be fed in as a 
        string. zone is the zone number in which to plot it. This is
        generally the total number of electrons transferred to get to 
        the state. zone should be fed in as a number (probably an
        integer). By default, puts a text label of the state (number)
        on the bar, but an alternate label can be fed in via the label
        keyword."""
        if label == None:
            label = state
        color = self._color
        self.ax.plot((zone-self._halfwidth, zone-self._textwidth),
                     (self._G[state],self._G[state]), 
                     color=color, linewidth=self._lw)
        self.ax.plot((zone+self._textwidth, zone+self._halfwidth),
                     (self._G[state],self._G[state]), 
                     color=color, linewidth=self._lw)
        if self._textposition == 'inline':
            self.ax.text(zone, self._G[state], label,
                         fontsize=self._fontsize,
                         horizontalalignment='center',
                         verticalalignment='center',
                         color=color)
                         #backgroundcolor='w')
        elif self._textposition == 'above':
            self.ax.text(zone, self._G[state] + self._text_vspace, label,
                         fontsize=self._fontsize,
                         horizontalalignment='center',
                         verticalalignment='bottom',
                         color=color)
        elif self._textposition == 'below':
            self.ax.text(zone, self._G[state] + -self._text_vspace, label,
                         fontsize=self._fontsize,
                         horizontalalignment='center',
                         verticalalignment='top',
                         color=color)

        elif self._textposition == 'slant':
            self.ax.text(zone-self._halfwidth*0.8, self._G[state]+self._text_vspace, label,
                         fontsize=self._fontsize,
                         horizontalalignment='left',
                         verticalalignment='bottom',
                         color=color,
                         rotation=60)


    def connect(self, state1, zone1, state2, zone2,barrier=0,dash=True):
        """Draws a dashed connector between stated positions and energies.
        Looks up the energies in the G dictionary. state1 and state2
        should be fed in as strings. zone1 and zone2 should be fed
        in as numbers."""
        if dash:
            m = '--'
        else:
            m = '-'
        if barrier==0:
           self.ax.plot((zone1+self._halfwidth, zone2-self._halfwidth),
                     (self._G[state1], self._G[state2]),
                     m, color=self._color,dashes=(2,2))
        if barrier>0:
            A=UnivariateSpline(np.array([zone1+self._halfwidth,0.5*(zone1+zone2),zone2-self._halfwidth]),np.array([self._G[state1],self._G[state1]+barrier,self._G[state2]]),k=2)
            x = np.linspace(zone1+self._halfwidth,zone2-self._halfwidth)
            self.ax.plot(x,A(x),'-', color=self._color)
