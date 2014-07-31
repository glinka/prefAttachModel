import numpy as np
import matplotlib.pyplot as plt
import sys
import time

from mpltools import style
from mpltools import layout

style.use('ggplot')

# integrate the SDE
# dX = h*X*dt + mu*X*dW
# using various stochastic integration schemes

def euler_maruyama(f, g, T, steps, W=None, X0=1, weak=False):
    # set up W and check for dimension agreement
    if W is not None and steps+1 != W.shape[0]:
        print 'input W not of appropriate length'
        exit()
    dt = 1.0*T/steps
    if W is None:
        if weak:
            W = np.power(dt, 0.5)*(2*np.random.randint(low=0, high=2, size=steps+1) - 1)
            print 'weak rw data generated'
        else:
            W = gen_rw(T, dt)[1]
    # run alg
    Xs = np.zeros(steps+1)
    Xs[0] = X0
    for i in range(1, steps+1):
        Xs[i] = Xs[i-1] + f(Xs[i-1])*dt + g(Xs[i-1])*(W[i]-W[i-1])
    # return [times, vals]
    return np.array([dt*np.array(range(steps+1)), Xs])

def milstein(f, g, g1, T, steps, W=None, X0=1):
    # set up W and check for dimension agreement
    if W is not None and steps+1 != W.shape[0]:
        print 'input W not of appropriate length'
        exit()
    dt = 1.0*T/steps
    if W is None:
        W = gen_rw(T, dt)[1]
    # run alg
    Xs = np.zeros(steps+1)
    Xs[0] = X0
    for i in range(1, steps+1):
        Xs[i] = Xs[i-1] + f(Xs[i-1])*dt + g(Xs[i-1])*(W[i]-W[i-1]) + g(Xs[i-1])*g1(Xs[i-1])*(np.power(W[i]-W[i-1], 2) - dt)/2.0
    # return [times, vals]
    return np.array([dt*np.array(range(steps+1)), Xs])
    

def gen_rw(T=100, dt=0.01):
    steps = int(T/dt)
    step_var = np.power(dt, 0.5)
    moves = step_var*np.random.normal(size=steps+1)
    locs = np.zeros(steps+1)
    locs[0] = 0
    locs[1] = moves[1]
    for i in range(2, steps+1):
        locs[i] = locs[i-1] + moves[i]
    # return array of [locations, times]
    return np.array([dt*np.array(range(0, steps+1)), locs])

def test_colors():
    fig = plt.figure()
    ax = fig.add_subplot(111)
    index = 0
    for c in plt.rcParams['axes.color_cycle']:
        ax.scatter(index, index, c=c, s=160, lw=0)
        index = index + 1
    plt.show(fig)

def plot_trajs(times, trajs):
    # expects trajs to have shape (ntrajs, nsteps)
    import matplotlib.cm as cm
    import matplotlib.colors as colors
    import matplotlib.colorbar as colorbar
    import matplotlib.gridspec as gs
    FONTSIZE = 48
    TICKSIZE = 36
    gspec = gs.GridSpec(6,6)
    fig = plt.figure()
    ax = fig.add_subplot(gspec[:6,:5])
    nsteps = times.shape[0]
    ntrajs = trajs.shape[0]
    # plot each step individually, coloring by y-value
    indices = np.array(range(ntrajs))
    for i in range(nsteps):
        ax.scatter(np.ones(ntraj)*times[i], np.sort(n*trajs[:,i]), c=indices, cmap='jet', lw=0, s=12, alpha=0.7)
    ax.set_xlim((times[0], times[-1]))
    ax.set_ylim((0, int(n*1.05*np.max(trajs))))
    ax.set_xlabel('t', fontsize=FONTSIZE)
    ax.set_ylabel('degree (from CIR SDE)', fontsize=FONTSIZE)
    # assume ymajorticklabels is not empty
    txtcolor = ax.get_ymajorticklabels()[0].get_color()
    ax.set_title(r'$\kappa = $' + str(kappa) + ', ' + r'$\rho = $' + str(rho), fontsize=FONTSIZE, color=txtcolor)
    ax.tick_params(axis='both', which='major', labelsize=TICKSIZE)
    ax.tick_params(axis='both', which='minor', labelsize=TICKSIZE)
    greyval = '0.96'
    ax.set_axis_bgcolor(greyval)
    # set up colorbar
    ax_cb = fig.add_subplot(gspec[:,5])
    colorbarnorm = colors.Normalize(vmin=0, vmax=100)
    cb = colorbar.ColorbarBase(ax_cb, cmap='jet', norm=colorbarnorm, orientation='vertical')
    ax_cb.tick_params(axis='both', which='major', labelsize=TICKSIZE)
    ax_cb.tick_params(axis='both', which='minor', labelsize=TICKSIZE)
    fig.text(0.78, 0.93, 'percentile', fontsize=FONTSIZE-4, color=txtcolor)
    plt.show()
    
def progress_bar(current, total, elapsed_time=None):
    perc = int((100.0*current)/total)
    percf = (100.0*current)/total
    bar = '\r['
    for i in range(perc):
        bar = bar + '|'
    for i in range(100-perc):
        bar = bar + ' '
    bar = bar + '] '
    if elapsed_time is not None:
        bar = bar + str(int(elapsed_time/(percf/100.0)) - int(elapsed_time)) + 's remaining'
    print bar,
    sys.stdout.flush()

if __name__=='__main__':
    # test_colors()
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--kappa', '-k', nargs=1, type=float, default=[1.0])
    parser.add_argument('--rho', '-r', nargs=1, type=float, default=[2.0])
    parser.add_argument('--n', '-n', nargs=1, type=int, default=[500])
    args = parser.parse_args()
    kappa = args.kappa[0]
    rho = args.rho[0]
    n = args.n[0]
    f = lambda x: kappa*(1 - x/rho)
    g = lambda x: np.power(2*x, 0.5)
    g1 = lambda x: np.power(2*x, -0.5)
    T = 1
    steps = np.power(2, 9)
    dt = 1.0*T/steps
    ntraj = n
    Xtrajs = np.zeros((ntraj, steps+1))
    start = time.clock()
    for i in range(ntraj):
        W=gen_rw(T, dt)[1]
        times, Xtrajs[i,:] = milstein(f, g, g1, T, steps, W=W, X0=2)
        progress_bar(i+1, ntraj, time.clock() - start)
    plot_trajs(times, Xtrajs)
    # ax.plot(times, XsM, c=c)
    # c1 = plt.rcParams['axes.color_cycle'][5]
    # c2 = plt.rcParams['axes.color_cycle'][4]
    # c3 = plt.rcParams['axes.color_cycle'][1]
    # c4 = plt.rcParams['axes.color_cycle'][0]
    # # ax.scatter(times, XsEM, c=c1, lw=0, s=12, alpha=0.7, label='em')
    # # ax.plot(times, XsEM, c=c1)
    # ax.scatter(times, XsM, c=c3, lw=0, s=12, alpha=0.7, label='m')
    # ax.plot(times, XsM, c=c3)
    # # ax.scatter(times, np.exp((kappa - np.power(mu, 2)/2)*times + mu*W), c=c2, lw=0, s=12, alpha=0.7, label='analytical')
    # # ax.plot(times, np.exp((kappa - np.power(mu, 2)/2)*times + mu*W), c=c2)
    # ax.grid(b=None)
    # ax.legend()
    # plt.show(fig)
