from matplotlib.ticker import FuncFormatter, FormatStrFormatter
from matplotlib import ticker, cm, colors, colorbar, gridspec as gs, pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D
from mpl_toolkits.mplot3d import proj3d

import numpy as np
import numpy.polynomial.legendre as npl
import os

import warnings

import util_fns as uf
from pca import pca

data_directory = './data/'

class FormatAxis:
    def __init__(self, ax, has_zaxis=True):
        self._ax = ax
        self._has_zaxis = has_zaxis

    def format(self, axis, data, format_string, offset=0, nticks=5):
        # dictionary of relevant functions/attributes
        d = {'x':{'ax':self._ax.xaxis, 'set-tick-pos':self._ax.set_xticks, 'set-tick-labels':self._ax.set_xticklabels, 'tick-pos':self._ax.xaxis.get_majorticklocs()},
             'y':{'ax':self._ax.yaxis, 'set-tick-pos':self._ax.set_yticks, 'set-tick-labels':self._ax.set_yticklabels, 'tick-pos':self._ax.yaxis.get_majorticklocs()}}

        if self._has_zaxis:
            d['z'] = {'ax':self._ax.zaxis, 'set-tick-pos':self._ax.set_zticks, 'set-tick-labels':self._ax.set_zticklabels, 'tick-pos':self._ax.get_zticks()}

        # tick positions are constant, regardless of offset
        minval, maxval = np.min(data), np.max(data)
        increment = (maxval - minval)/(nticks-1)
        d[axis]['set-tick-pos']([minval + i*increment for i in range(nticks)])
        # subtract offset from data if using
        if offset != 0:
            if offset < 0:
                offset_str =  '- %1.2f' % abs(offset) # unicode dash u'\u2012'
            else:
                offset_str = '+ %1.2f' % abs(offset)
            # go through the terrible process of figuring out where to put the damn thing
            loc = {'x':0, 'y':0, 'z':0}
            for i,key in enumerate(d.keys()):
                if key is axis:
                    if axis is 'x':
                        loc[key] = np.min(d[key]['tick-pos']) - 0.00*(np.max(d[key]['tick-pos']) - np.min(d[key]['tick-pos']))
                    else:
                        loc[key] = np.max(d[key]['tick-pos']) + 0.00*(np.max(d[key]['tick-pos']) - np.min(d[key]['tick-pos']))
                else:
                    if key is 'x':
                        loc[key] = np.max(d[key]['tick-pos']) + 0.2*(np.max(d[key]['tick-pos']) - np.min(d[key]['tick-pos']))
                    else:
                        loc[key] = np.min(d[key]['tick-pos']) - 0.2*(np.max(d[key]['tick-pos']) - np.min(d[key]['tick-pos']))
            if self._has_zaxis:
                self._ax.text(loc['x'], loc['y'], loc['z'], offset_str, fontsize=12) #maxval-0.05*(maxval-minval)
            else:
                self._ax.text(loc['x'], loc['y'], offset_str, fontsize=12) #maxval-0.05*(maxval-minval)
            data = data - offset
        # set axis tick labels
        minval = np.min(data)
        d[axis]['set-tick-labels']([format_string % (minval + i*increment) for i in range(nticks)])

def align_degs(cpi_times, cpi_degs, nocpi_times, nocpi_degs, ci=None):
    mutual_times = []
    cpi_indices = []
    nocpi_indices = []
    cpi_index = 0
    # manual search
    if ci is None:
        n = cpi_times.shape[0]
        m = nocpi_times.shape[0]
        j = 0
        for i in range(n):
            while j < m and cpi_times[i] > nocpi_times[j]:
                j += 1
            if j < m:
                if cpi_times[i] == nocpi_times[j]:
                    cpi_indices.append(i)
                    nocpi_indices.append(j)
                    mutual_times.append(cpi_times[i])
    else:
        for t in cpi_times:
            if t % ci == 0:
                if t <= nocpi_times[-1]:
                    cpi_indices.append(cpi_index)
                    nocpi_indices.append(t/ci - 1)
                    mutual_times.append(t)
            cpi_index = cpi_index + 1
    return [np.array(mutual_times), cpi_degs[cpi_indices, :], nocpi_degs[nocpi_indices, :]]


def degree_evo_figs():
    """Plots degree vs. vertex (projection along time axis) for both lopsided and ER initialized networks, showing that both evolve to the same stationary degree sequence.

    **To generate data:**::
    
        mpirun -n 20 ./pref_attach -n 100 -steps 10000000 -init_type lopsided
        mpirun -n 20 ./pref_attach -n 100 -m 5050 -steps 10000000 -init_type erdos
        cp ./paper-data/erdos*.out ./manuscript-materials/data
        cp ./paper-data/lopsided*.out ./manuscript-materials/data

    """
    # import data
    degsl = np.genfromtxt(data_directory + 'lopsided-degs.out', delimiter=',', skip_header=False)
    degse = np.genfromtxt(data_directory + 'erdos-degs.out', delimiter=',', skip_header=False)
    times = np.genfromtxt(data_directory + 'erdos-times.out', skip_header=False)
    # set number of times to keep, adjusted so that the approach to steady state is clear
    ntimes = 101
    times = times[:ntimes]
    n = degsl.shape[1]

    # create color-related objects
    colornorm = colors.Normalize(vmin=0, vmax=ntimes-1)
    colormap = cm.ScalarMappable(norm=colornorm, cmap='jet')

    # plot lopsided degrees WITHOUT colorbar (only one colorbar is needed among the two plots as they share the same colornorm
    fig = plt.figure()
    ax = fig.add_subplot(111)
    for i in range(ntimes):
        ax.plot(np.arange(n), np.sort(degsl[i]), c=colormap.to_rgba(i))
    ax.set_xlabel('Vertex')
    ax.set_ylabel('Degree')

    # plot erdos degrees
    gspec = gs.GridSpec(6,6)
    fig = plt.figure()
    ax = fig.add_subplot(gspec[:,:5])
    for i in range(ntimes):
        ax.plot(np.arange(n), np.sort(degse[i]), c=colormap.to_rgba(i))
    ax.set_xlabel('Vertex')
    ax.set_ylabel('Degree')
    # add colorbar
    axcb = fig.add_subplot(gspec[:,5])
    colorbarnorm = colors.Normalize(vmin=0, vmax=times[-1])
    cb = colorbar.ColorbarBase(axcb, cmap='jet', norm=colorbarnorm, orientation='vertical')
    cb.ax.get_yaxis().labelpad = 45
    cb.ax.set_ylabel('Step', rotation=270)
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # set tick labels, note the hardcoded scaling
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    warnings.warn("Using hardcoded labels, may need to be manually adjusted")
    ticks = np.empty(6)
    ticks[:5] = times[:100:20]
    ticks[5] = times[-1]
    cb.set_ticks(ticks)
    cb.set_ticklabels(ticks/1000000)
    fig.text(0.85, 0.93, '1e6')

    plt.show()


def newton_figs():
    """Plots both the decrease of the coarse Newton-GMRES residual over iterations, and a comparison of the final Newton solution with the final direct simulatoin state.

    **To generate data (may take repeated attempts):**::

        mpirun -n 10 ./coarse_ng -n 100 -m 50000 -nms 1000000 -omw 100000 -ci 10000 -proj_step 500000
        mpirun -n 20 ./pref_attach -input_filename ./newton_data/xs.csv -n 100
        cp ./newton_data/*.csv ./manuscript-materials/data
        cp ./fromfile_init/*.csv ./manuscript-materials/data

    """
    # import data
    resids = np.genfromtxt(data_directory + 'resids.csv', skip_header=False)
    degn = np.genfromtxt(data_directory + 'xs.csv', delimiter=',', skip_header=False) # newton
    degs = np.genfromtxt(data_directory + 'degs.csv', delimiter=',', skip_header=True) # simulation

    niters = resids.shape[0]
    n = degn.shape[1]
    
    # plot direct simulation and newton-grmes solution
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.arange(n), np.sort(degn[-1]), c='r', label='Newton-GMRES solution', lw=3)
    ax.plot(np.arange(n), degs[-1], c='b', label='Direct simulation stationary state', lw=3)
    ax.legend(loc=2)
    ax.set_xlabel('Vertex')
    ax.set_ylabel('Degree')
    
    # plot error vs. iteration
    fs = 48
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.arange(niters), resids)
    ax.scatter(np.arange(niters), resids, lw=3, s=50)
    ax.set_xlabel('Iteration', fontsize=fs)
    ax.set_ylabel(r'$\parallel F(x^{(k)}) \parallel_2$', fontsize=fs)

    plt.show()

def triangle_fig():
    """Plots the triangle count in three different systems which were initialized with identical degree sequences but different triangle counts, showing that this secondary property is quickly slaved to a slow manifold. Shows evolution over both n^2 and n^3 timescales in one figure with a broken x-axis.

    **To generate data:**::

        ./transients-main -n 99 -m 10000 -s 3000000 -kappa 1
        cp ./paper-data/other-runs/*tris.csv ./manuscript-materials/data
        cp ./paper-data/other-runs/tris_times.csv ./manuscript-materials/data

    """
    # import data
    hh_tris = np.genfromtxt(data_directory + 'hh_tris.csv', delimiter=',', skip_header=False)
    t_tris = np.genfromtxt(data_directory + 'tri_tris.csv', delimiter=',', skip_header=False)
    r_tris = np.genfromtxt(data_directory + 'rando_tris.csv', delimiter=',', skip_header=False)
    times = np.genfromtxt(data_directory + 'tris_times.csv', delimiter=',', skip_header=False)

    # # plot data on single figure with broken x-axis
    # n^2
    fig = plt.figure()
    ax1 = fig.add_subplot(121)
    ifinal = 50 # total number of points to take on for n^2 portion of plot
    ax1.semilogy(times[:ifinal], hh_tris[:ifinal], c='b', lw=7)
    ax1.semilogy(times[:ifinal], t_tris[:ifinal], c='r', lw=7)
    ax1.semilogy(times[:ifinal], r_tris[:ifinal], c='g', lw=7)

    # n^3
    ax2 = fig.add_subplot(122, sharey=ax1)
    istart = ifinal # the index to start the n^3 plot from (makes sense to start from the end of the previous section, i.e. at 'ifinal'
    ax2.semilogy(times[istart:], hh_tris[istart:], lw=7)
    ax2.semilogy(times[istart:], t_tris[istart:], c='r', lw=7)
    ax2.semilogy(times[istart:], r_tris[istart:], c='g', lw=7)

    # formatting
    ax1.set_ylabel('Triangle count')
    ax1.text(0.95, -0.17, 'Step', transform=ax1.transAxes)

    ax2.ticklabel_format(axis='x', style='sci', scilimits=(-1, 1))
    ax1.ticklabel_format(axis='x', style='sci', scilimits=(-1, 1))
    warnings.warn("Using hardcoded axis limits, may need to be manually adjusted")
    ax1.set_ylim((1e4, 1e7))
    ax1.set_xlim(right=times[ifinal-1])
    ax2.set_xlim(left=times[istart])
    ax2.tick_params(axis='y', labelsize=0)
    ax1.spines['right'].set_visible(False)
    ax2.spines['left'].set_visible(False)
    ax1.tick_params(labeltop='off')
    fig.subplots_adjust(left=0.09, bottom=0.17, right=0.97, wspace=0.08, top=0.95)

    # add diagonal lines to indicate broken axis
    d = 0.015
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False, lw=2)
    ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # bottom-right diagonal
    ax1.plot((1 - d, 1 + d), (1-d, 1+d), **kwargs)  # top-right diag
    kwargs.update(transform=ax2.transAxes)  # switch to the right axes
    ax2.plot((-d, +d), (-d, +d), **kwargs)  # bottom-left diagonal
    ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal

    plt.show()


def cpi_figs():
    """Plots 3d evolution of both CPI-accelerated degree sequence and direct simulation, showing visual agreement. Also plots CPI error over time, as measured by discrepancies in degree sequence from direct simulation.

    **To generate data:**::

        mpirun -n 50 ./pref_attach -n 100 -m 50000 -nms 1000000 -omw 100000 -ci 10000 -pStep 500000
        mpirun -n 50 ./pref_attach -n 100 -m 50000 -nms 1000000 -omw 100000 -ci 10000 -pStep 500000 -noinit
        cp ./paper-data/other-runs/withinittimes0.csv ./manuscript-materials/data
        cp ./paper-data/other-runs/withinitpre_proj_degs0.csv ./manuscript-materials/data
        cp ./paper-data/other-runs/noinittimes0.csv ./manuscript-materials/data
        cp ./paper-data/other-runs/noinitpre_proj_degs0.csv ./manuscript-materials/data
    
    """
    init = "noinit"
    times = np.genfromtxt(data_directory + init + 'times0.csv', delimiter=',', skip_header=True)
    degs = np.genfromtxt(data_directory + init + 'pre_proj_degs0.csv', delimiter=',', skip_header=True)

    init = "withinit"
    times_cpi = np.genfromtxt(data_directory + init + 'times0.csv', delimiter=',', skip_header=True)
    degs_cpi = np.genfromtxt(data_directory + init + 'pre_proj_degs0.csv', delimiter=',', skip_header=True)

    gsize = degs.shape[1]

    # # 3d degree comparison
    # lenghthen x-axis
    x_scale=2
    y_scale=1
    z_scale=1
    scale=np.diag([x_scale, y_scale, z_scale, 1.0])
    scale=scale*(1.0/scale.max())
    scale[3,3]=1.0
    def short_proj():
        return np.dot(Axes3D.get_proj(ax), scale)

    # plot 3d degrees
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.get_proj=short_proj
    ntake = 650 # final index for direct simulation degrees
    slice_size = 15 # only plot ever 'slice_size' degree sequences
    ts = []
    tscale = 1e7 # scaling to make labels look nicer
    # plot direct simulation
    for i, t in enumerate(times[:ntake:slice_size]):
        ax.scatter(np.ones(gsize)*t/tscale, np.arange(gsize), degs[i*slice_size], c='b', zorder=2)
        ts.append(t)
    ntake = 400 # final index for cpi simulation degrees (less than previously due to projective step)
    # plot cpi
    for i, t in enumerate(times_cpi[:ntake:slice_size]):
        ax.scatter(np.ones(gsize)*t/tscale, np.arange(gsize), degs_cpi[i*slice_size], c='r', zorder=1)

    # formatting
    warnings.warn("Using hardcoded axis limits and labels, may need to be manually adjusted")
    ax.set_xlim((0, 1))
    ax.set_ylim((0, gsize))
    ax.set_zlim((100, 2500))
    formatter = FormatAxis(ax)
    formatter.format('x', np.array(ts)/tscale, '%1.1f', nticks=5)
    ax.text(0.7, -15, -500, '1e7')
    ax.set_xlabel('\n\n\nStep')
    ax.set_ylabel('\n\nVertex')
    ax.set_zlabel('\n\n\nDegree')
    ax.set_yticks((25, 75))
    ax.set_zticks((100, 1200, 2500))
    ax.set_xlim((0, 0.7))
    plt.show()

    # # cpi error
    # find degree sequences that are collected at identical times
    mutual_times, degs_cpi, degs = align_degs(times_cpi, degs_cpi, times, degs)

    fig = plt.figure()
    ax2 = fig.add_subplot(111)
    scale = np.linalg.norm(degs_cpi[0]) # should be constant throughout
    ax2.plot(mutual_times, np.linalg.norm(np.abs(degs_cpi - degs), axis=1)/scale, color='r')

    # formatting
    fs = 48
    ax2.set_xlabel('Step', fontsize=fs)
    ax2.set_ylabel('Relative error', fontsize=fs)
    ax2.set_xlim((0, mutual_times[-1]))
    ax2.ticklabel_format(axis='y', style='sci', scilimits=(-1, 1))
    fig.subplots_adjust(left=0.1, bottom=0.11)
    ax2.legend(loc=2)

    plt.show()

def coeff_fig():
    """Plots a sorted degree sequence and its polynomial fit (used in lifting/restriction diagram). Also plots evolution of coefficients over time for both a cpi routine and a direct simulation (currently not included in manuscript).

    **To generate data:**::

        mpirun -n 50 ./pref_attach -n 100 -m 50000 -nms 1000000 -omw 100000 -ci 10000 -pStep 500000
        mpirun -n 50 ./pref_attach -n 100 -m 50000 -nms 1000000 -omw 100000 -ci 10000 -pStep 500000 -noinit
        cp ./paper-data/other-runs/withinitfitted_coeffs0.csv ./manuscript-materials/data
        cp ./paper-data/other-runs/withinitpre_proj_degs0.csv ./manuscript-materials/data
        cp ./paper-data/other-runs/withinittimes0.csv ./manuscript-materials/data
        cp ./paper-data/other-runs/noinittimes0.csv ./manuscript-materials/data
        cp ./paper-data/other-runs/noinitfitted_coeffs0.csv ./manuscript-materials/data
        cp ./paper-data/other-runs/noinitpre_proj_degs0.csv ./manuscript-materials/data

    """
    init = "noinit"
    times = np.genfromtxt(data_directory + init + 'times0.csv', delimiter=',', skip_header=True)
    degs = np.genfromtxt(data_directory + init + 'pre_proj_degs0.csv', delimiter=',', skip_header=True)
    fitted_coeffs = np.genfromtxt(data_directory + init + 'fitted_coeffs0.csv', delimiter=',', skip_header=True)

    init = "withinit"
    times_cpi = np.genfromtxt(data_directory + init + 'times0.csv', delimiter=',', skip_header=True)
    degs_cpi = np.genfromtxt(data_directory + init + 'pre_proj_degs0.csv', delimiter=',', skip_header=True)
    fitted_coeffs_cpi = np.genfromtxt(data_directory + init + 'fitted_coeffs0.csv', delimiter=',', skip_header=True)

    gsize = degs.shape[1]

    fig = plt.figure()
    ax = fig.add_subplot(111)
    index = 50 # arbitrary index, should show some curvature
    ax.scatter(np.arange(gsize), degs[index], label='Sorted degree sequence')
    ax.plot(np.arange(gsize), npl.legval(np.arange(gsize), fitted_coeffs[index]), c='r', label='Polynomial fit')
    ax.set_xlim((0,gsize))
    ax.set_xlabel('Vertex')
    ax.set_ylabel('Degree')
    ax.legend(loc=4, fontsize=42)
    fig.subplots_adjust(left=0.1, bottom=0.11)

    ncoeffs = fitted_coeffs.shape[1]
    for i in range(ncoeffs):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.scatter(times, fitted_coeffs[:,i], label='Direct simulation')
        ax.scatter(times_cpi, fitted_coeffs_cpi[:,i], label='CPI')
        ax.set_xlabel('Step')
        ax.set_ylabel(r'$c_' + str(i) + '$')
        ax.legend()

    plt.show()


def pca_rho_kappa_embedding_figs():
    """Performs PCA on a collection of network stationary states arising from a range of $m$ and $\kappa$ values, plots projection of data along PC1 and PC2, also plots variances.

    **To generate data:**::

        ./rho_kappa_embedding 0 2000
        cp ./embedding_data/rho_kappa_graph_embeddings.csv ./manuscript-materials/data
        cp ./embedding_data/rho_kappa_params.csv ./manuscript-materials/data

    """
    # project graph embeddings with PCA
    embeddings = np.genfromtxt(data_directory + 'rho_kappa_graph_embeddings.csv', delimiter=',')
    k = 6
    pcs, variances = pca(embeddings, k)
    projections = np.dot(pcs.T, embeddings.T) # (k, n) array

    # plot variances
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(np.arange(1,k+1), variances)
    ax.semilogy(np.arange(1,k+1), variances)
    ax.set_xlabel(r'$i$')
    ax.set_ylabel(r'$\sigma^2_i$')

    # # plot projection along first and second principal component
    params = np.genfromtxt(data_directory + 'rho_kappa_params.csv', delimiter=',')
    # color by log(kappa)
    fs = 64
    s = 50
    fig = plt.figure()
    ax = fig.add_subplot(111)
    c = ax.scatter(projections[0], projections[1], c=np.log10(params[:,1]), s=s)
    cb = fig.colorbar(c)
    cb.set_label(r'$\log(\kappa)$', fontsize=fs)
    ax.set_xlabel(r'$w_1$', fontsize=fs)
    ax.set_ylabel(r'$w_2$', fontsize=fs)
    ax.set_xlim((1.05*np.min(projections[0]), 1.4*np.max(projections[0])))
    ax.set_ylim(bottom=1.05*np.min(projections[1]))
    formatter = FormatAxis(ax, has_zaxis=False)
    formatter.format('x', projections[0], '%d', nticks=3)
    formatter.format('y', projections[1], '%d', nticks=3)
    fig.subplots_adjust(bottom=0.15)

    # color by rho
    fig = plt.figure()
    ax = fig.add_subplot(111)
    c = ax.scatter(projections[0], projections[1], c=params[:,0], s=s)
    cb = fig.colorbar(c)
    cb.set_label(label=r'$\frac{2m}{n}$', fontsize=1.5*fs)
    ax.set_xlabel(r'$w_1$', fontsize=fs)
    ax.set_ylabel(r'$w_2$', fontsize=fs)
    ax.set_xlim((1.05*np.min(projections[0]), 1.4*np.max(projections[0])))
    ax.set_ylim(bottom=1.05*np.min(projections[1]))
    formatter = FormatAxis(ax, has_zaxis=False)
    formatter.format('x', projections[0], '%d', nticks=3)
    formatter.format('y', projections[1], '%d', nticks=3)
    fig.subplots_adjust(bottom=0.15)
    plt.show()


def dmaps_rho_kappa_embedding_figs():
    """Embeds collection of stationary networks generated over a range of $m$ and $\kappa$ values using DMAPS. Plots resulting embedding in \phi_e1 and \phi_e2.

    **To generate data:**::

        ./rho_kappa_embedding 0 2000
        cp ./embedding_data/dmaps_rho_kappa_embedding_eigvects.csv ./manuscript-materials/data
        cp ./embedding_data/rho_kappa_params.csv ./manuscript-materials/data

    """
    # import data
    eigvects = np.genfromtxt(data_directory + 'dmaps_rho_kappa_embedding_eigvects.csv', delimiter=',')
    params = np.genfromtxt(data_directory +  'rho_kappa_params.csv', delimiter=',')
    e1 = 1 # eigvect-one indeix
    e2 = 4 # eigvect-two indeix
    fs = 64
    s = 50

    # color by log(kappa)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    c = ax.scatter(eigvects[e1], eigvects[e2], c=np.log10(params[:,1]), s=s)
    cb = fig.colorbar(c)
    cb.set_label(r'$\log(\kappa)$', fontsize=fs)
    ax.set_xlabel(r'$\Phi_' + str(e2) + '$', fontsize=fs)
    ax.set_ylabel(r'$\Phi_' + str(e1) + '$', fontsize=fs)
    ax.set_xlim((1.05*np.min(eigvects[e1]), 1.05*np.max(eigvects[e1])))
    ax.set_ylim((1.05*np.min(eigvects[e2]), 1.2*np.max(eigvects[e2])))
    warnings.warn("Using hardcoded axis and labels, may need to be manually adjusted")
    ax.set_xticks((-6e-4, -1e-4, 4e-4))
    ax.set_xticklabels(('-6', '-1', '4'))
    ax.set_yticks((7.5e-4, 0, -7.5e-4))
    ax.set_yticklabels(('-7.5', '0', '7.5'))
    ax.text(0.99, -0.06, r'$10^4$', transform=ax.transAxes) # x
    ax.text(-0.07, 1.0, r'$10^4$', transform=ax.transAxes) # y
    fig.subplots_adjust(bottom=0.15)

    # color by rho
    fig = plt.figure()
    ax = fig.add_subplot(111)
    c = ax.scatter(eigvects[e1], eigvects[e2], c=2*params[:,0]/50.0, s=s)
    cb = fig.colorbar(c)
    cb.set_label(label=r'$\frac{2m}{n}$', fontsize=1.5*fs)
    ax.set_xlim((1.05*np.min(eigvects[e1]), 1.05*np.max(eigvects[e1])))
    ax.set_ylim((1.05*np.min(eigvects[e2]), 1.2*np.max(eigvects[e2])))
    ax.set_xlabel(r'$\Phi_' + str(e2) + '$', fontsize=fs)
    ax.set_ylabel(r'$\Phi_' + str(e1) + '$', fontsize=fs)
    ax.ticklabel_format(axis='x', style='sci', scilimits=(-1, 1))
    warnings.warn("Using hardcoded axis and labels, may need to be manually adjusted")
    ax.set_xticks((-6e-4, -1e-4, 4e-4))
    ax.set_xticklabels(('-6', '-1', '4'))
    ax.set_yticks((7.5e-4, 0, -7.5e-4))
    ax.set_yticklabels(('-7.5', '0', '7.5'))
    ax.text(0.99, -0.06, r'$10^4$', transform=ax.transAxes) # x
    ax.text(-0.07, 1.0, r'$10^4$', transform=ax.transAxes) # y
    fig.subplots_adjust(bottom=0.15)

    plt.show()


def dmaps_trajectory_embedding_figs():
    """Embeds two different system trajectories, one initialized as a lopsided graph, the other as ER, using DMAPS. Plots 3d embedding in $\phi_1$, $\phi_2$ and $\phi_3$.

    **To generate data:**::

        ./graph-embedding-motifs 200 10000 5 1000 1
        cp ./embedding_data/dmaps_many_embedding_eigvects.csv ./manuscript-materials/data
        cp ./embedding_data/dmaps_many_embedding_eigvals.csv ./manuscript-materials/data

    """
    # import dat data
    eigvals = np.genfromtxt('./data/dmaps_many_embedding_eigvals.csv', delimiter=',', skip_header=False)
    eigvects = np.genfromtxt('./data/dmaps_many_embedding_eigvects.csv', delimiter=',', skip_header=False)
    warnings.warn("Using hardcoded 'ntypes' value, may need to be manually adjusted")
    ntypes = 2
    n = eigvects.shape[1]
    n_pertype = n/ntypes
    eigvals = np.abs(eigvals)
    sorted_indices = np.argsort(eigvals)[::-1]
    eigvals = eigvals[sorted_indices]
    eigvects = eigvects[sorted_indices, :]
    nvects = eigvals.shape[0]
    cmaps = ['autumn', 'jet', 'spectral', 'winter', 'summer', 'PrOr', 'RdBu', 'RdYlBu', 'RdYlGn']
    cs = ['b', 'r', 'g', 'y', 'c', 'm', 'k', 'b']
    i = 1 # \phi_i
    j = 2 # \phi_j
    k = 3 # \phi_k
    t = 1 # DMAPS time
    xvals = np.power(eigvals[i], t)*eigvects[i]
    yvals = np.power(eigvals[j], t)*eigvects[j]
    zvals = np.power(eigvals[k], t)*eigvects[k]
    # # plot all trajectories, making the size of the marker increase with time
    # also circle the final state
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, projection='3d')
    for p in range(ntypes):
        ax.scatter(xvals[p*n_pertype:(p+1)*n_pertype], yvals[p*n_pertype:(p+1)*n_pertype], zvals[p*n_pertype:(p+1)*n_pertype], c=cs[p], lw=0, alpha=0.3, s=np.arange(1, n_pertype)*1000/(n_pertype))
        ax.scatter(xvals[(p+1)*n_pertype-1], yvals[(p+1)*n_pertype-1], zvals[(p+1)*n_pertype-1], lw=50, edgecolor=cs[p], c=cs[p], s=1, alpha=1)
    ax.set_xlabel('\n\n\n' + r'$\Phi_ ' + str(i) + '$')
    ax.set_ylabel('\n\n\n' + r'$\Phi_ ' + str(j) + '$')
    ax.set_zlabel('\n\n\n' + r'$\Phi_ ' + str(k) + '$')
    formatter = FormatAxis(ax)
    formatter.format('x', xvals, '%1.2f', nticks=3)
    formatter.format('y', yvals, '%1.2f', nticks=3)
    formatter.format('z', zvals, '%1.2f', nticks=3)

    # # plot an enlarged section of the previous figure, showing the final few states of all trajectories meeting
    # circle average final state
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, projection='3d')

    q = 60 # the number of points to plot in the enlarged figure
    # need to collect data that is actually plotted for axis formatting
    xdata = np.empty(ntypes*q)
    ydata = np.empty(ntypes*q)
    zdata = np.empty(ntypes*q)
    for p in range(ntypes):
        xdata[q*p:q*(p+1)] = xvals[(p+1)*n_pertype-q:(p+1)*n_pertype]
        ydata[q*p:q*(p+1)] = yvals[(p+1)*n_pertype-q:(p+1)*n_pertype]
        zdata[q*p:q*(p+1)] = zvals[(p+1)*n_pertype-q:(p+1)*n_pertype]

        navg = 20
        avg_finalx = np.average(xvals[(p+1)*n_pertype-navg:(p+1)*n_pertype])
        avg_finaly = np.average(yvals[(p+1)*n_pertype-navg:(p+1)*n_pertype])
        avg_finalz = np.average(zvals[(p+1)*n_pertype-navg:(p+1)*n_pertype])

        ax.scatter(xvals[(p+1)*n_pertype-q:(p+1)*n_pertype], yvals[(p+1)*n_pertype-q:(p+1)*n_pertype], zvals[(p+1)*n_pertype-q:(p+1)*n_pertype], c=cs[p], lw=0, alpha=0.3, s=np.arange(1, n_pertype)[-q:]*1000/(n_pertype))
        ax.scatter(avg_finalx, avg_finaly, avg_finalz, lw=75, edgecolor=cs[p], c=cs[p])


    plt.locator_params(nbins=3) # does this do anything?
    ax.set_xlabel('\n\n\n' + r'$\Phi_ ' + str(i) + '$')
    ax.set_ylabel('\n\n\n' + r'$\Phi_ ' + str(j) + '$')
    ax.set_zlabel('\n\n\n' + r'$\Phi_ ' + str(k) + '$')
    warnings.warn("Using hardcoded axis limits, may need to be manually adjusted")
    # ax.set_xlim((-0.0025, 0.0015))
    # ax.set_zlim((-0.0012, 0.0012))

    plt.show()


def plot_all_figs():
    """Plots all the figures"""
    degree_evo_figs()
    newton_figs()
    triangle_fig()
    cpi_figs()
    coeff_fig()
    pca_rho_kappa_embedding_figs()
    dmaps_rho_kappa_embedding_figs()
    dmaps_trajectory_embedding_figs()
    

if __name__=='__main__':
    plot_all_figs()
