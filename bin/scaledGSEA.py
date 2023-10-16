#!/usr/bin/env python
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import glob
import os
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline as spline
import sys
from collections import OrderedDict
import argparse
from pathlib import Path


color_maps = {

    '1': ['#edd717','#ef4348', '#1474a8', .75],
    '2': ["#0EAD69","#29335C","#DD6031", .75],
}


def argparser():
    '''Process command line arguments.'''

    parser = argparse.ArgumentParser(description='Plot stacked GSEA curves.')
    parser.add_argument('-f', required=True, nargs='+', dest='paths', help='Path(s) to "*.xls" file.')
    parser.add_argument('-e', required=True, nargs='+', dest='experiments', help='Experiment names - the order must match "-f" inputs.')
    parser.add_argument('-c', type=str, default='1', dest='colormap', help=f'Colormap (1..{len(color_maps)})')
    parser.add_argument('-t', type=str, default='GSEA', dest='title', help='The desired title at the top of the figure.')
    parser.add_argument('-n', type=int, default=56650, dest='genes', help='Total number of genes quantified.')
    args = parser.parse_args()


    # Make sure the number of experiment names and xlsx files match.
    n_experiments = len(args.experiments)
    n_paths = len(args.paths)
    if n_experiments != n_paths:
        sys.exit(f"Number of experiments ({n_experiments}) does not match number of input files ({n_paths}")
    else:
        print('\nProcessing...\n')
        for experiment, path in zip(args.experiments, args.paths):
            print(f"{experiment}:\t{path}")

    # Ensure all paths are valid.
    for path in args.paths:
        if not Path(path).exists():
            sys.exit(f"{path} not found.")

    # Select a backup color_map if input selection doesn't exist in color_map dictionary.
    if args.colormap not in color_maps:
        print(f"{args.color_map} not a valid color map. Choices are: {[i for i in color_maps]}")
        print(f"Using color_map 1 ({color_maps['1']})")
        args.color_map = '1'
    
    return args


def import_df(path, n_genes):
    '''Return rank and score arrays from input GSEA xlsx file.'''
    
    try:
        df = pd.read_table(path, index_col=0, usecols=["RANK IN GENE LIST", "RUNNING ES"])
    except ValueError:
        print(f"{path} not a valid GSEA output file. Make sure file has the following columns: 'RANK IN GENE LIST' and 'RUNNING ES'")
        sys.exit()

    df.loc[n_genes] = 0  # Ensure the plot ends at 0
    df.loc[0] = 0
    df = df.sort_index()

    rank = np.array(df.index)
    es = np.array(df['RUNNING ES'])
    return rank, es


def plot_curve(ax, color, label, rank, es, n_genes, spline_points=10000):

    x = np.linspace(0, n_genes, spline_points)  
    y = spline(rank, es, k=1)
    ax.plot(x, y(x), linewidth=3, label=label, c=color)

    return x, y(x )

    
def plot_barcode(ax, color, rank, grid_spec, grid_spec_ypos):

    ax2 = plt.subplot(grid_spec[grid_spec_ypos, :], sharex=ax)
    ax2.axis('off')
    plt.plot(rank, [0]*len(rank), '|', ms=30, c=color)
    

def fig_format(ax, xlen, title):
    
    plt.xlim([-1, xlen + 1])
    ax.set_xlabel('Rank', fontsize=23, fontweight='bold')
    ax.set_ylabel('ES', fontsize=23, fontweight='bold')
    ax.set_title(title, fontsize=26, fontweight='bold')
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.axhline(0, c='k', linewidth=4)
    ax.grid()
    ax.legend(loc='upper right', fancybox=True, fontsize=16)
   # plt.subplots_adjust(left=0.2, right=0.8)

plt.rc('axes', linewidth=2)
plt.rc('font', weight='bold')

args = argparser()

grid_ylen = 13 + len(args.paths)
gs = GridSpec(grid_ylen, 1)    
fig = plt.figure(figsize=(12,8))
ax = plt.subplot(gs[:10,:])

spline_points = 1000
ranks = np.zeros([len(args.paths), spline_points]) # Initialize an array
color_map = color_maps[args.colormap]

for n, (path, experiment) in enumerate(zip(args.paths, args.experiments)):
    try:
        color = color_map[n]
    except KeyError:
        color = '.75'

    rank, es = import_df(path, args.genes)
    x, y = plot_curve(ax=ax, color=color, label=experiment, rank=rank, es=es, n_genes=args.genes, spline_points=spline_points)
    ranks[n] = y
    plot_barcode(ax=ax, color=color, rank=rank, grid_spec=gs, grid_spec_ypos=13+n)

sorted_curves = np.argsort(np.min(ranks, 1)) # Sorted indices of minimum values from low to high
for i in range(len(ranks)):
    min_curve = sorted_curves[i]
    y0 = ranks[min_curve]
    color = color_map[min_curve]
    ax.fill_between(x, y0, np.zeros(spline_points), color=color, alpha=.06) # Fill AUC


fig_format(ax=ax, xlen=args.genes, title=args.title)
title = args.title.replace(' ', '_')
plt.savefig(f'{title}.stacked.png', dpi=300)



