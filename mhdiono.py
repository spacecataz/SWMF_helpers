#!/usr/bin/env python3

from spacepy.pybats import bats, rim
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
import glob
from sys import stdout
import os
import multiprocessing as mp
from argparse import ArgumentParser, RawDescriptionHelpFormatter, BooleanOptionalAction
import matplotlib
import numpy as np
import time
import scipy.interpolate as irp

parser = ArgumentParser(description=__doc__,
                        formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-n", "--nthread", type=int, default=1,
                    help="If >1, plots will be made in parallel using N " +
                    "threads")
parser.add_argument('-mv', '--mhdvars', default=None, nargs='+', type=str, help='Sets vars to plot')
parser.add_argument('-iv', '--ionovars', default=None, nargs='+', type=str, help='Sets vars to plot')
parser.add_argument('-y', '--noy', action='store_true', help='turn off y plotting')
parser.add_argument('-z', '--noz', action='store_true', help='turn off z plotting')
parser.add_argument('-xrng', '--xrng', default=(-64, 16), help='Set the ' +
                    'spatial extent in the GSM X direction via two space-' +
                    'separated numbers with no spaces (e.g., -xrng -10 10)',
                    type=float, nargs=2)
parser.add_argument('-yrng', '--yrng', default=(-32, 32),
                    help='Set the spatial extent in the GSM Y direction via ' +
                    'two space-separated numbers with no spaces ' +
                    '(e.g., -yrng -10 10)', type=float, nargs=2)
parser.add_argument('-zrng', '--zrng', default=(-32, 32),
                    help='Set the spatial extent in the GSM Z direction via ' +
                    'two space-separated numbers with no spaces ' +
                    '(e.g., -zrng -10 10)', type=float, nargs=2)
parser.add_argument("-m", "--mag", action="store_true",
                    help="Add magnetic field lines to y=0 slice")
parser.add_argument("-g", "--grid", action="store_true", help="add grid outline to plots")
parser.add_argument('-hemi', '--hemi', default=None, nargs='+', type=str, help='Set hemisphere to plot')
parser.add_argument('-colat', '--colat_max', default=40, type=int, help='Sets colat max for plot')
parser.add_argument('-tr', '--trace', action='store_true', help='Trace field lines from ionosphere to magnetosphere')
parser.add_argument('-a', '--alpha', default=0.8, type=float, help='Set alpha for traced values')
parser.add_argument('-icm', '--iono_cmap', default='gist_heat', type=str, help='Set ionosphere colors')
parser.add_argument('-lc', '--lat_contour', default=None, type=float, help='Set lat contour of ionosphere for z=0 cut (in degrees)')
# Ex: mhdiono.py -n 56 -mv p ux -iv e-flux sigmah -hemi n_ -tr -m -lc 5 

# Handle arguments:
args = parser.parse_args()

maxz = {'e-flux': 30, 'ave-e': 40, 'sigmah': 80, 'sigmap': 40, 'e-flux-diffe': 30, 'ave-e-diffe': 40,
        'e-flux-diffi': 3, 'ave-e-diffi': 2, 'e-flux-elec': 30, 'ave-e-elec': 40, 'jr': 1.5,
        'e-flux-bbnd': 8E-1, 'ave-e-bbnd': 1, 'rt p': 3E-7, 'rt rho': 1E-17, 'i-long': 100, 'i-short': 100,
        'jphi': 4, 'rt 1/b': 1, 'rt pe': 3E-8, 'rt ppar': 3E-8,
        'rt pepar': 3E-8, 'rt pperp': 3E-8, 'rt peperp': 3E-9, 'electron-loss-cone': 0.25, 'ion-loss-cone': 0.3,
        'efluxdiff': 14E-3, 'aveediff': 100, 'phi': 200, 'im-jr': 1.5}

title = {'e-flux-elec': r'$\Phi_{Energy, Electron}$', 'ave-e-elec': r'$<E_{Electron}>$', 'jr': r'$J_r$', 'sigmah': r'$\sigma_H$',
         'sigmap': r'$\sigma_P$', 'e-flux-diffi': r'$\Phi_{Energy, Ion}$', 'ave-e-diffi': r'$<E_{Ion}>$',
         'e-flux-diffe': r'$\Phi_{Energy, DiffE}$', 'ave-e-diffe': r'$<E_{DiffE}>$', 'rt pe': r'MHD $P_e$',
         'rt pepar': r'MHD $P_{e,||}$', 'rt ppar': r'$MHD P_{i,||}$', 'rt p': r'MHD $P_i$', 'rt rho': r'MHD rho',
         'electron-loss-cone': r'Electron Loss Cones', 'ion-loss-cone': r'Ion Loss Cones', 'jphi': r'$J_{\phi}$',
         'efluxdiff': r'Total - Diffuse $\Phi_{Energy, Electron}$', 'aveediff': r'Total - Diffuse $<E_{Electron}>$',
         'rt 1/b': r'Open vs Closed Field Lines', 'phi': r'$\phi$', 'e-flux-bbnd': r'$\Phi_{Energy, Broadband}$',
         'ave-e-bbnd': r'$<E_{Broadband}>$', 'im-jr': r'IM $J_r$', 'e-flux': r'$\Phi_{Energy}$', 'ave-e': r'$<E>$'}

matplotlib.rcParams['figure.max_open_warning'] = 1E6
matplotlib.use('Agg')

outdir = './mhdiono/'
if not os.path.exists(outdir):
    os.mkdir(outdir)

nmhdvars = len(args.mhdvars)
if args.ionovars is None:
    nionovars = 0
else:
    nionovars = len(args.ionovars)
nvars = max(nmhdvars, nionovars)
if args.hemi is None:
        args.hemi = ['n_', 's_']

def plot_cuts(yfile=None, zfile=None, ionofile=None, iFile=0, nFiles=1):
    if not args.noz:
        mz = bats.Bats2d(zfile)
        mz.calc_temp(units='keV')
    if not args.noy:
        my = bats.Bats2d(yfile)
        my.calc_temp(units='keV')
    if ionofile is not None:
        iono = rim.Iono(ionofile)

    mhd_width = 2
    if args.noy:
        mhd_width = 1
    if args.noz:
        mhd_width = 1
    iono_width = 0
    if ionofile is not None:
        iono_width = len(args.hemi)
    
    ncbars = 1 + 1 * (nionovars > 0)

    fig = plt.figure(figsize=(4 + 6 * (mhd_width + iono_width), 1 + 6 * nvars))
    grid  = gs.GridSpec(nrows=nvars, ncols=4 * (mhd_width + iono_width) + ncbars, figure=fig, 
                        left=0.10 / (mhd_width + iono_width) + 0.02, right=0.99, bottom=0.05, top=0.95, wspace=2, hspace=0.2)
    uselog = True
    ax = []
    for i, var in enumerate(args.mhdvars):
        j = 0
        match var:
            case 't' | 'te':
                lims = (0.00001, 20)
                cmap = 'viridis'
                uselog=True
            case 'Ti/Te':
                lims = (0.001, 100)
                cmap = 'jet'
                if not args.noz:
                    mz['Ti/Te'] = mz['t'] / mz['te']
                if not args.noy:
                    my['Ti/Te'] = my['t'] / my['te']
            case 'p' | 'pe' | 'ppar' | 'pperp' | 'pepar' | 'peperp':
                lims = (0.00001, 100)
                cmap = 'viridis'
                uselog=True
            case 'ppar/p' | 'pepar/pe' | 'ppar/pperp' | 'pepar/peperp':
                lims = (0.01, 100)
                cmap = 'seismic'
                uselog = True
            case 'rho':
                lims = (0.01, 100)
                cmap = 'viridis'
            case 'ux' | 'uy' | 'uz' | 'un':
                lims = (-400, 400)
                cmap = 'seismic'
                uselog = False
            case 'b1z' | 'b1y' | 'b1x' | 'b1n':
                lims = (-20, 20)
                cmap = 'seismic'
                uselog = False
            case 'wn':
                if not args.noz:
                    mz.calc_vort()
                if not args.noy:
                    my.calc_vort()
                lims = (-0.5, 0.5)
                cmap = 'seismic'
                uselog = False
            case _:
                print(f'unrecognized var {var}')
                continue
        if not args.noz:
            if args.noy:
                spacing = 5
                use_cbar = True
            else:
                spacing = 4
                use_cbar = False
            ax.append(fig.add_subplot(grid[i, j:j + spacing]))
            j += spacing
            var_z = var
            if var_z == 'wn':
                var_z = 'wz'
            if var_z == 'un':
                var_z = 'uy'
            if var_z == 'b1n':
                var_z = 'b1z'
            if var_z == 'ppar/p':
                mz['ppar/p'] = mz['ppar'] / mz['p']
            if var_z == 'pepar/pe':
                mz['pepar/pe'] = mz['pepar'] / mz['pe']
            if var_z == 'ppar/pperp':
                mz['ppar/pperp'] = mz['ppar'] / mz['pperp']
            if var_z == 'pepar/peperp':
                mz['pepar/peperp'] = mz['pepar'] / mz['peperp']
            if args.grid:
                mz.add_grid_plot(target=ax[-1], do_fill=False, do_label=False)
            fig_1, ax_1, cont, cbarz = mz.add_contour('y', 'x', var_z, target=ax[-1], cmap=cmap, dolog=uselog, xlim=args.yrng, ylim=args.xrng,
                    zlim=lims, add_cbar=use_cbar, extend='both')
            ax[-1].set_xlabel('Y $(R_E)$', fontsize=18)
            ax[-1].set_ylabel('X $(R_E)$', fontsize=18)
            if use_cbar:
                if lims[0] < 0:
                    cbarz.set_ticks(np.linspace(lims[0], lims[1], 7))
                cbarz.ax.tick_params(labelsize=16)
                if 'units' in mz[var_z].attrs:
                    clabel = f"{var}, ({mz[var_z].attrs['units']})"
                else:
                    clabel = f"{var}"
                cbarz.set_label(clabel, fontsize=18)

        if args.trace and 'status' in mz.keys() and ionofile is not None:
            if iono['n_' + args.ionovars[i]].min() >= 0: 
                interpolator = irp.RegularGridInterpolator((iono['n_theta'][:, 0], iono['n_psi'][0, :]), iono['n_' + args.ionovars[i]])
                mask = mz['status'] == 3
                lats, lons = 90 - mz['lat1'][mask], mz['lon1'][mask]
                interp_values = interpolator((lats, lons))
                if args.ionovars[i] in ['e-flux', 'e-flux-diffe', 'e-flux-elec', 'e-flux-bbnd', 'e-flux-diffi']:
                    interp_values *= 1000
                value_mask = interp_values > maxz[args.ionovars[i]] / 5
                xs, ys = mz['x'][mask][value_mask], mz['y'][mask][value_mask]
                ax[-1].scatter(ys, xs, c=interp_values[value_mask], cmap='gist_heat', s=4, alpha=args.alpha, vmin=0, vmax=maxz[args.ionovars[i]])

        if args.lat_contour is not None and 'status' in mz.keys():
            mask = mz['status'] == 3
            max_lat = np.quantile(mz['lat1'][mask], 0.95)
            levels = np.arange(args.colat_max, max_lat, args.lat_contour)
            plt.tricontour(mz['y'][mask], mz['x'][mask], mz['lat1'][mask], levels=levels, colors='k', linewidths=1, alpha=0.25)

        ax[-1].invert_xaxis()
            
        if not args.noy:
            ax.append(fig.add_subplot(grid[i, j:j + 5]))
            var_y = var
            if var_y == 'wn':
                var_y = 'wy'
            if var_y == 'un':
                var_y = 'uz'
            if var_y == 'b1n':
                var_y = 'b1y'
            if var_y == 'ppar/p':
                my['ppar/p'] = my['ppar'] / my['p']
            if var_y == 'pepar/pe':
                my['pepar/pe'] = my['pepar'] / my['pe']
            if var_y == 'ppar/pperp':
                my['ppar/pperp'] = my['ppar'] / my['pperp']
            if var_y == 'pepar/peperp':
                my['pepar/peperp'] = my['pepar'] / my['peperp']
            if args.grid:
                my.add_grid_plot(target=ax[-1], do_fill=False, do_label=False)
            fig_2, ax_2, cont, cbary = my.add_contour('x', 'z', var_y, target=ax[-1], cmap=cmap, dolog=uselog, xlim=args.xrng, ylim=args.zrng,
                    zlim=lims, add_cbar=True, extend='both')
            ax[-1].set_xlabel('X $(R_E)$', fontsize=18)
            ax[-1].set_ylabel('Z $(R_E)$', fontsize=18)
            cbary.ax.tick_params(labelsize=16)
            if lims[0] < 0:
                cbary.set_ticks(np.linspace(lims[0], lims[1], 7))
            if 'units' in my[var_y].attrs:
                clabel = f"{var}, ({my[var_y].attrs['units']})"
            else:
                clabel = f"{var}"
            cbary.set_label(clabel, fontsize=18)
            if args.mag:
                my.add_b_magsphere(ax[-1], DoOpen=True)

    if ionofile is not None:
        if 'n_jx' in iono.keys():
            iono.calc_j()
            iono['n_jphi'] /= 1000000.
            iono['s_jphi'] /= 1000000.
            iono['n_jphi'].attrs['units'] = r'$\mu A/m^2$'
            iono['s_jphi'].attrs['units'] = r'$\mu A/m^2$'
        if 'n_e-flux-mono' in iono.keys():
            iono['n_e-flux-elec'] = iono['n_e-flux-mono']
            iono['s_e-flux-elec'] = iono['s_e-flux-mono']
            iono['n_ave-e-elec'] = iono['n_ave-e-mono']
            iono['s_ave-e-elec'] = iono['s_ave-e-mono']
        if 'n_e-flux-elec' not in iono.keys():
            iono['n_e-flux-elec'] = iono['n_e-flux']
            iono['s_e-flux-elec'] = iono['s_e-flux']

        for i, var in enumerate(args.ionovars):
            for j, hemi in enumerate(args.hemi):
                use_cbar = hemi == args.hemi[-1]
                kwargs = {'target': fig, 'maxz': maxz[var], 'extend':'both', 'lines': False, 'add_cbar': use_cbar}
                gridstart = 4 * mhd_width + 1 + j * 4
                gridend = gridstart + 4 + (1 * use_cbar)
                if var in ['e-flux-diffe', 'e-flux-elec', 'e-flux', 'e-flux-bbnd', 'e-flux-diffi']:
                    iono[hemi+var] *= 1000
                    iono[hemi+var].attrs['units'] = '$mW/m^2$'
                if var in ['jr', 'jphi', 'efluxdiff', 'aveediff', 'phi', 'im-jr']:
                    ax.append(fig.add_subplot(grid[i, gridstart:gridend], polar=True))
                    _, _, _, cb = iono.add_cont(hemi + var, cmap='seismic', loc=ax[-1], label=title[var], **kwargs, max_colat=args.colat_max)
                elif var in ['i-long', 'i-short']:
                    ax.append(fig.add_subplot(grid[i, gridstart:gridend], polar=True))
                    _, _, _, cb = iono.add_cont(hemi + var, cmap='Greens', loc=ax[-1], label=title[var], **kwargs, max_colat=args.colat_max)
                else:
                    ax.append(fig.add_subplot(grid[i, gridstart:gridend], polar=True))
                    _, _, _, cb = iono.add_cont(hemi + var, cmap=args.iono_cmap, loc=ax[-1], label=title[var], **kwargs, max_colat=args.colat_max)
            cb.ax.tick_params(labelsize=16)
            cb.set_label(iono[hemi+var].attrs['units'], fontsize=18)

    for a in ax:
        a.tick_params(axis='both', which='major', labelsize=16)
        #for tick in a.xaxis.get_major_ticks():
        #    tick.label.set_fontsize(16)
        #for tick in a.yaxis.get_major_ticks():
        #    tick.label.set_fontsize(16)
    fig.suptitle(f'{mz.attrs["time"]}', y=0.99, fontsize=22)
    fig.savefig(outdir + f'{mz.attrs["time"]:%Y%m%d_%H%M%S}_xycuts.png')
    fig.clf()
    plt.close(fig)

    percent_done = float(iFile + 1) / float(nFiles) * 100.0
    return [iFile+1, percent_done]


def print_progress(percent_obj):
    i = percent_obj[0]
    percent = percent_obj[1]
    if (i % args.nthread == 0):
        stdout.write(200 * '\b')
        stdout.write('\t{:5.1f}% complete...'.format(percent))
        stdout.write('time left: {:5.1f} seconds...'.format((100.0 - percent) / percent * (time.time() - start)))
        stdout.write('[' + '\u2588' * int(percent) + ' ' * (100 - int(percent)) + ']')
        stdout.flush()

if __name__ == '__main__':
    if not args.noz:
        zFiles = glob.glob('./GM/z=0_mhd*.out')
        if(len(zFiles) == 0):
            zFiles = glob.glob('./GM/z=0_var*.out')
        if(len(zFiles) == 0):
            zFiles = glob.glob('./GM/z=0_ful*.out')
        zFiles.sort()
        zFiles = zFiles[1:]
    if not args.noy:
        yFiles = glob.glob('./GM/y=0_mhd*.out')
        if(len(yFiles) == 0):
            yFiles = glob.glob('./GM/y=0_var*.out')
        if(len(yFiles) == 0):
            yFiles = glob.glob('./GM/y=0_ful*.out')
        yFiles.sort()
        yFiles = yFiles[1:]
    if args.noz:
        zFiles = yFiles
    if args.noy:
        yFiles = zFiles

    ionoFiles = glob.glob('./IE/it*')
    ie_start = './IE/'
    if(len(ionoFiles) == 0):
        ionoFiles = glob.glob('./IE/plots/it*')
        ie_start = './IE/plots/'
    ionoFiles.sort()

    nFiles = len(yFiles)

    print('\tFiles: working on {} files...'.format(nFiles))
    stdout.write('\t{:5.1f}% complete...'.format(0))

    if args.nthread > 1:
        # Create pool of workers.
        pool = mp.Pool(args.nthread)
        start = time.time()

        for i, files in enumerate(zip(yFiles, zFiles)):
            ionofile = ie_start + 'it' + files[0][-21:-7] + '000.idl'
            ionofile = ionofile.replace('-', '_')
            if ionofile not in ionoFiles:
                ionofile = ionofile + '.gz'
            if ionofile not in ionoFiles and args.ionovars is not None:
                continue
            if args.ionovars is None:
                ionofile = None

            x = pool.apply_async(plot_cuts, args=(files[0], files[1], ionofile),
                             kwds={'iFile': i, 'nFiles': nFiles},
                             callback=print_progress)

        pool.close()
        pool.join()

        stdout.write(200 * '\b')
        stdout.write('\t{:5.1f}% complete...'.format(0))
        stdout.write('time left: 0.0 seconds...')
        stdout.write('[' + '\u2588' * 100 + ']')

    else:
        start = time.time()

        for i, files in enumerate(zip(yFiles, zFiles)):
            ionofile = ie_start + 'it' + files[0][-21:-7] + '000.idl'
            ionofile = ionofile.replace('-', '_')
            if ionofile not in ionoFiles:
                ionofile = ionofile + '.gz'
            if ionofile not in ionoFiles and args.ionovars is not None:
                print(f'No ionosphere file found for {ionofile} and ionovars specified, skipping...')
                continue
            if args.ionovars is None:
                ionofile = None

            per = plot_cuts(files[0], files[1], ionofile, iFile=i, nFiles=nFiles)
            print_progress(per)
