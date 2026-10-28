#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 12:56:40 2026

@author: gvill
"""
import LT.box as B
import numpy as np

fname = 'yield_n_xsec/jml_theory_xsec/pm_120_jml_theory.txt'
f = B.get_file(fname)

dtr = np.pi/180

phi_avg = f['phi_avg']
phi_2 = f['phi_wrong']

phi_avg_hist = B.histo(phi_avg,range=(-180,180),bins=100)
phi_2_hist = B.histo(phi_2,range=(0,180),bins=100)

# Set up bins for our distribution (e.g., 16 bins for 16 compass directions)
number_of_bins = 19
bins = np.linspace(-np.pi, np.pi, number_of_bins + 1)

# Compute frequency of data in each bin
counts, integrity_bins = np.histogram(phi_avg*dtr, bins=bins)

# Center the bars appropriately on the bins
bin_centers = (bins[:-1] + bins[1:]) / 2
width = (2 * np.pi) / number_of_bins

# Generate the rose chart
fig, ax = B.pl.subplots(subplot_kw={'projection': 'polar'})
bars = ax.bar(bin_centers, counts, width=width, bottom=0.0, edgecolor='black')

# second histogram
# Set up bins for our distribution (e.g., 16 bins for 16 compass directions)
number_of_bins = 19
bins = np.linspace(0, 2*np.pi, number_of_bins + 1)

# Compute frequency of data in each bin
counts, integrity_bins = np.histogram(phi_2*dtr, bins=bins)

# Center the bars appropriately on the bins
bin_centers = (bins[:-1] + bins[1:]) / 2
width = (2 * np.pi) / number_of_bins

# Generate the rose chart
fig, ax = B.pl.subplots(subplot_kw={'projection': 'polar'})
bars = ax.bar(bin_centers, counts, width=width, bottom=0.0, edgecolor='black')
