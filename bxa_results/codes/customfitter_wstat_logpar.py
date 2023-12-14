"""

Fits custom models automatically

adapted from xagnfitter.py from https://github.com/JohannesBuchner/BXA/blob/master/examples/sherpa/xagnfitter.py

How to run this script:

1) create filenames.txt:

	contains 3 space-separated columns, giving .pha/.pi filename, lower and upper energy range
	for example:
	
	acis.pi  0.5 8
	nu_A.fits 5 77
	nu_B.fits 5 77

2) load ciao (source path/to/ciao/bin/ciao.sh). You installed ciao and bxa, right?

3) run:

	$ python3 customfitter.py

4) analyse results:

	This will create a bunch of output files starting with multiple_out_*.
	
	intrinsic_photonflux.dist.gz contains the posterior samples as rows 
	(each row is a equally probable solution). 
	
	The columns are:
	
	* redshift 
	* rest-frame intrinsic (unabsorbed) flux in the 2-10keV band
	* absorbed flux in the observed band 2-8keV
	* source normalisation (log)
	* photon index: Prior is 1.95 +- 0.15, so check if it differs from that.
	* log(NH): absorbing column density (from 20 to 26)
	* f_scat: normalisation of soft, unobscured powerlaw
	* apec norm (if did not set WITHAPEC=0)
	* apec temperature (if did not set WITHAPEC=0)
	* redshift (if not fixed to a single value)
	* and one background normalisation parameter for each spectrum
	
	You can use cosmolopy to convert each flux and redshifts to a luminosity.

"""
from sherpa.astro.ui import *
import os
import sys
import numpy

import bxa.sherpa as bxa
# from bxa.sherpa.background.pca import auto_background
from sherpa.models.parameter import Parameter
# from bxa.sherpa.background.models import ChandraBackground #, SwiftXRTBackground, SwiftXRTWTBackground
# from bxa.sherpa.background.fitters import SingleFitter
# from bxa.sherpa.cachedmodel import CachedModel

import logging
logging.basicConfig(filename='bxa.log',level=logging.DEBUG)
logFormatter = logging.Formatter("[%(name)s %(levelname)s]: %(message)s")
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
consoleHandler.setLevel(logging.INFO)
logging.getLogger().addHandler(consoleHandler)
import warnings
logging.getLogger('sherpa.plot').setLevel(logging.ERROR)
warnings.filterwarnings("ignore", message='.*displayed errorbars.*')
warnings.filterwarnings("ignore", message='.*Clearing background convolved model.*', append=True)
warnings.filterwarnings("ignore", message='.*apply_rmf[(]apply_arf.*', append=True)



import argparse
import os

parser = argparse.ArgumentParser(description="customfitter")

parser.add_argument('--filename','-fn',dest='filename',
        help='filename')

parser.add_argument('--model','-md',dest='model',help='model name')

parser.add_argument('--fixnorm','-fxn',dest='fixnorm',default='F',help='fixing normalization of multiple observations?')

args = parser.parse_args()

fielname = args.filename

model = args.model

fixnorm = args.fixnorm

if fixnorm == 'T':
    fixnorm = True
else:
    fixnorm = False

print(fixnorm)
# if not os.path.exists('filenames.txt'):
# 	print("ERROR: No filenames.txt found.")
# 	print()
# 	print(__doc__)

ids = []
galabso = None
clean()
set_xlog()
set_ylog()
set_stat('wstat')
set_xsabund('wilm')
set_xsxsect('vern')

for id, line in enumerate(open(fielname), start=0):
    if id == 0:
        # first line is the directory for the source, not per-obs, save as prefix
        prefix = line.strip()
    else:
        print(id)
        filename, elo, ehi = line.strip().split()
        elo, ehi = float(elo), float(ehi)
        print('loading spectral file "%s" as id=%s ...' % (filename, id))
        load_pha(id, filename)
        ignore_id(id, None, None)
        notice_id(id, elo,  ehi)
        set_analysis(id, 'ener', 'counts')
        # unsubtract(id)
        group_counts(id, num=5)
        ids.append(id)

# model = 'mekal'

prefix = prefix + model + '_out_'

# if len(ids) == 1:
# 	prefix = filename + '_' + model + '_out_'
# else:
# 	prefix = filename + '_' + model + '_multiple_out_'



model_dict = {'mekal':   {'n_comp':1,'xsmodel':'xsmekal.mek', 'model_prefix':'mek', 'model_par':'kT'},
            'powerlaw':  {'n_comp':1,'xsmodel':'xspowerlaw.pl','model_prefix':'pl', 'model_par':'PhoIndex'},
            'bb':        {'n_comp':1,'xsmodel':'xsbbodyrad.bb', 'model_prefix':'bb', 'model_par':'kT'},
            'powerlawbb':{'n_comp':2,'xsmodel1':'xspowerlaw.pl','model_prefix1':'pl', 'model_par1':'PhoIndex', 'xsmodel2':'xsbbodyrad.bb', 'model_prefix2':'bb', 'model_par2':'kT'}, 
            'twomekal':  {'n_comp':2,'xsmodel1':'xsmekal.mek1', 'model_prefix1':'mek1', 'model_par1':'kT','xsmodel2':'xsmekal.mek2', 'model_prefix2':'mek2', 'model_par1':'kT'}}

n_comp = model_dict[model]['n_comp']
priors = []
parameters = []
print('setting source model ...')
print('creating priors')

for id in ids:
    if n_comp ==1:
        model_src = 'xstbabs.tbabs1 * ' + model_dict[model]['xsmodel'] + str(id)
        set_model(id, model_src)
        # convmodel = get_model(id)
        model_comp = get_model_component(model_dict[model]['model_prefix'] + str(id))
        if id==1:
            guess(model_comp)
        # norm_guess = numpy.log10(model_comp.norm.val)
        # print('norm_guess',norm_guess)
        if fixnorm==False:
            set_par(model_comp.norm, val=1., min=1e-8, max=100)
            srcnorm = Parameter(f"src{id}", 'lognorm', 0., -8., 2.,-8., 2.) 
            model_comp.norm = 10**srcnorm
            # print(srcnorm)
            priors += [bxa.create_uniform_prior_for(srcnorm)]
            # priors += [bxa.create_uniform_prior_for(model_comp.norm)]
            parameters += [srcnorm]
        elif fixnorm==True:
            if id==1:
                set_par(model_comp.norm, val=1., min=1e-8, max=100)
                srcnorm = Parameter(f"src{id}", 'lognorm', 0., -8., 2.,-8., 2.) 
                model_comp.norm = 10**srcnorm
                # print(srcnorm)
                priors += [bxa.create_uniform_prior_for(srcnorm)]
                # priors += [bxa.create_uniform_prior_for(model_comp.norm)]
                parameters += [srcnorm]
            elif id>1:
                if model == 'mekal':
                    model_comp.norm = mek1.norm
                if model == 'powerlaw':
                    model_comp.norm = pl1.norm
                if model == 'bb':
                    model_comp.norm = bb1.norm
        # print(mek1.norm)
        if id >1:
            if model == 'mekal':
                model_comp.kT = mek1.kT
            if model == 'powerlaw':
                model_comp.PhoIndex = pl1.PhoIndex
            if model == 'bb':
                model_comp.kT = bb1.kT
    elif n_comp == 2:
        model_src = 'xstbabs.tbabs1 * (' + model_dict[model]['xsmodel1'] + str(id) + ' + ' + model_dict[model]['xsmodel2'] + str(id) + ')'
        set_model(id, model_src)
        # convmodel = get_model(id)
        model_comp1 = get_model_component(model_dict[model]['model_prefix1'] + str(id))
        if id==1:
            guess(model_comp1)
        # norm1_guess = numpy.log10(model_comp1.norm.val)
        set_par(model_comp1.norm, val=1., min=1e-8, max=100)
        srcnorm1 = Parameter(f'src1{id}', 'lognorm', 0., -8., 2.,-8., 2.) 
        model_comp1.norm = 10**srcnorm1
        priors += [bxa.create_uniform_prior_for(srcnorm1)]
        # priors += [bxa.create_uniform_prior_for(model_comp.norm)]
        parameters += [srcnorm1]

        model_comp2 = get_model_component(model_dict[model]['model_prefix2'] + str(id))
        if id==1:
            guess(model_comp2)
        # norm2_guess = numpy.log10(model_comp2.norm.val )
        set_par(model_comp2.norm, val=1., min=1e-8, max=100)
        srcnorm2 = Parameter(f'src2{id}', 'lognorm', 0., -8., 2.,-8., 2.) 
        model_comp2.norm = 10**srcnorm2
        priors += [bxa.create_uniform_prior_for(srcnorm2)]
        # priors += [bxa.create_uniform_prior_for(model_comp.norm)]
        parameters += [srcnorm2]
        if id >1:
            if model == 'twomekal':
                model_comp1.kT = mek11.kT
                model_comp2.kT = mek21.kT
            if model == 'powerlawbb':
                model_comp1.PhoIndex = pl1.PhoIndex
                model_comp2.kT = bb1.kT
         
    # print(get_model(id))

if model == 'mekal':
    set_par(mek1.kT, val=1., min=0.1, max=30)
    mekkT = Parameter('src1', 'logkT', 0, -1, numpy.log10(30), -1,  numpy.log10(30))
    mek1.kT = 10**mekkT
    priors += [bxa.create_uniform_prior_for(mekkT)]
    parameters += [mekkT]
if model == 'powerlaw':
    set_par(pl1.PhoIndex, val=1., min=-3., max=10.)
    priors += [bxa.create_uniform_prior_for(pl1.PhoIndex)]
    parameters += [pl1.PhoIndex]
if model == 'bb':
    set_par(bb1.kT, val=1., min=1e-3, max=10)
    bbkT = Parameter('src1', 'logkT', 0, -3, 1, -3, 1)
    bb1.kT = 10**bbkT
    priors += [bxa.create_uniform_prior_for(bbkT)]
    parameters += [bbkT]
if model == 'twomekal':
    set_par(mek11.kT, val=0.3, min=0.1, max=30)
    mek1kT = Parameter('src11', 'logkT', numpy.log10(0.3), -1, numpy.log10(30), -1,  numpy.log10(30))
    mek11.kT = 10**mek1kT
    priors += [bxa.create_uniform_prior_for(mek1kT)]
    parameters += [mek1kT]

    set_par(mek21.kT, val=2, min=0.1, max=30)
    mek2kT = Parameter('src21', 'logkT', numpy.log10(2), -1, numpy.log10(30), -1,  numpy.log10(30))
    mek21.kT = 10**mek2kT
    priors += [bxa.create_uniform_prior_for(mek2kT)]
    parameters += [mek2kT]
if model == 'powerlawbb':
    set_par(pl1.PhoIndex, val=1., min=-3., max=10.)
    priors += [bxa.create_uniform_prior_for(pl1.PhoIndex)]
    parameters += [pl1.PhoIndex]
    set_par(bb1.kT,  val=1., min=1e-3, max=10)
    bbkT = Parameter('src1', 'logkT', 0, -3, 1, -3, 1)
    bb1.kT = 10**bbkT
    priors += [bxa.create_uniform_prior_for(bbkT)]
    parameters += [bbkT]

set_par(tbabs1.nH, val=0.1, min=1.e-3, max=1.e4)
srcnh = Parameter('src1', 'lognH', 21., 19., 26., 19., 26.)
tbabs1.nh = 10**(srcnh - 22)
parameters += [srcnh]
priors += [bxa.create_uniform_prior_for(srcnh)]


# 
# apec with L(2-10keV) = 1e42 erg/s
# z      norm 10keV    norm 2keV
# 0.1    40e-6         150e-6
# 0.5    2e-6          6e-6
# 1      0.7e-6        2e-6
# 3      0.15e-6       0.5e-6
# ================================
# z -->  0.5e-6/ z**2  2e-6/z**2
# 



assert len(priors) == len(parameters), 'priors: %d parameters: %d' % (len(priors), len(parameters))


#################
# BXA run
priorfunction = bxa.create_prior_function(priors = priors)

print("Priors:", priors)
print("Params:", parameters)
print('running BXA ...')

solver = bxa.BXASolver(prior=bxa.create_prior_function(priors), parameters=parameters, outputfiles_basename=prefix)
results = solver.run(resume=False)

# solver = bxa.BXASolver(id = ids[0], otherids = tuple(ids[1:]),
# 	prior = priorfunction, parameters = parameters, 
# 	outputfiles_basename = prefix)
# results = solver.run(
# 	resume=True, n_live_points = int(os.environ.get('NLIVEPOINTS', 400)),
# 	frac_remain=0.5,
# )

try:
	from mpi4py import MPI
	if MPI.COMM_WORLD.Get_rank() > 0:
		sys.exit(0)
except Exception as e:
	pass

rows = results['samples']

import matplotlib
from matplotlib import pyplot as plt

for id in ids:
    print('plotting spectrum ...')
    set_analysis(id,'ener','counts')
    # subtract(id)
	# b = get_bkg_fit_plot(id)
	# numpy.savetxt(prefix + 'bkg_'+str(id)+'.txt.gz', numpy.transpose([b.dataplot.x, b.dataplot.y, b.modelplot.x, b.modelplot.y]))
    # m = get_fit_plot(id)
    # numpy.savetxt(prefix + 'src_notsubtract_'+str(id)+'.txt.gz', numpy.transpose([m.dataplot.x, m.dataplot.y, m.modelplot.x, m.modelplot.y]))

    # ypreds = []
    # for i, row in enumerate(rows):
    #     sys.stdout.write("%d/%d (%.2f%%)\r" % (i, len(rows), (i + 1)*100./ len(rows)))
    #     sys.stdout.flush()
    #     for p, v in zip(parameters, row):
    #         p.val = v

    #     m = get_fit_plot(id)
    #     ypreds.append(m.modelplot.y)

    # ylo, ymid, yhi = numpy.percentile(ypreds, [15.87, 50, 84.13], axis=0)
    # numpy.savetxt(prefix + '/src_notsubtract_'+str(id)+'.txt.gz', numpy.transpose([m.dataplot.x, m.dataplot.y, m.dataplot.yerr, m.modelplot.x, ylo, ymid, yhi]))


    group_counts(id, 15)
    wstat_model = get_fit_plot(id).modelplot

    # ui.show_model()
    
    # fig = plt.figure(figsize=(6, 4))
    # gs = fig.add_gridspec(2, 1, hspace=0.05,height_ratios=[3, 1])

    # ax1 = fig.add_subplot(gs[0])
    # plt.sca(ax1)
    # # plot_data(id=id, xlog=True, ylog=True, overplot=True, clearwindow=False,)
    # # plot_fit(id=id,  xlog=True, ylog=True, overplot=True, clearwindow=False,)
    # plot_data(id=id, xlog=True, ylog=True, overplot=True, clearwindow=False,)
    # wstat_model.plot(xlog=True, ylog=True, overplot=True, clearwindow=False,)
    # plt.xlim([0.5,8.0])
    # plt.xscale('log')
    # plt.yscale('log')
    # # plt.gcf().axes[0].set_title(plot_title)
    # ax1.set_ylabel('Counts/s/keV')

    # ax2 = fig.add_subplot(gs[1], sharex=ax1)
    # plt.sca(ax2)

    # # exactly replicates plot_ratio
    # dataplot = get_data_plot(id)
    # plt.errorbar(x=dataplot.x, y=dataplot.y/wstat_model.y, xerr=dataplot.xerr/2, yerr=dataplot.yerr/wstat_model.y, fmt='.')
    # plt.yscale('log')
    # ax2.get_xaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    # ax2.set_xticks([0.6, 1, 2, 3, 4, 5, 6, 7], minor=True)
    # ax2.set_xlabel('Energy (keV)')
    # ax2.set_ylabel('Data/Model')
    # plt.savefig(prefix +'/spectrum_nosubtract.png', dpi=300)


    subtract(id)
	# b = get_bkg_fit_plot(id)
	# numpy.savetxt(prefix + 'bkg_'+str(id)+'.txt.gz', numpy.transpose([b.dataplot.x, b.dataplot.y, b.modelplot.x, b.modelplot.y]))
    # m = get_fit_plot(id)
    # numpy.savetxt(prefix + 'src_'+str(id)+'.txt.gz', numpy.transpose([m.dataplot.x, m.dataplot.y, m.modelplot.x, m.modelplot.y]))

    ypreds = []
    for i, row in enumerate(rows):
        sys.stdout.write("%d/%d (%.2f%%)\r" % (i, len(rows), (i + 1)*100./ len(rows)))
        sys.stdout.flush()
        for p, v in zip(parameters, row):
            p.val = v

        m = get_fit_plot(id)
        ypreds.append(m.modelplot.y)

    ylo, ymid, yhi = numpy.percentile(ypreds, [15.87, 50, 84.13], axis=0)
    numpy.savetxt(prefix + '/src_subtract_'+str(id)+'.txt.gz', numpy.transpose([m.dataplot.x, m.dataplot.xerr, m.dataplot.y, m.dataplot.yerr, m.modelplot.x, ylo, ymid, yhi]))

    fig = plt.figure(figsize=(6, 4))
    gs = fig.add_gridspec(2, 1, hspace=0.05,height_ratios=[3, 1])

    ax1 = fig.add_subplot(gs[0])
    plt.sca(ax1)
    plot_data(id=id, xlog=True, ylog=True, overplot=True, clearwindow=False,)
    wstat_model.plot(xlog=True, ylog=True, overplot=True, clearwindow=False,)
    plt.xlim([0.5,8.0])
    plt.xscale('log')
    plt.yscale('log')
    # plt.gcf().axes[0].set_title(plot_title)
    ax1.set_ylabel('Counts/s/keV')

    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    plt.sca(ax2)

    # exactly replicates plot_ratio
    dataplot = get_data_plot(id)
    plt.errorbar(x=dataplot.x, y=dataplot.y/wstat_model.y, xerr=dataplot.xerr/2, yerr=dataplot.yerr/wstat_model.y, fmt='.')
    plt.yscale('log')
    ax2.get_xaxis().set_minor_formatter(matplotlib.ticker.ScalarFormatter())
    ax2.set_xticks([0.6, 1, 2, 3, 4, 5, 6, 7], minor=True)
    ax2.set_xlabel('Energy (keV)')
    ax2.set_ylabel('Data/Model')
    plt.savefig(prefix +f'/spectrum_{id}.png', dpi=300)

# # compute 2-10keV intrinsic luminosities?
# print("calculating intrinsic fluxes and distribution of model spectra")
# # calculate restframe intrinsic flux
# convmodel = get_model(id)
# set_model(id, convmodel)

# r = []
# for i, row in enumerate(rows):
#     sys.stdout.write("%d/%d (%.2f%%)\r" % (i, len(rows), (i + 1)*100./ len(rows)))
#     sys.stdout.flush()
#     # z = 0.
# 	# z = redshift.val if hasattr(redshift, 'val') else redshift
#     for p, v in zip(parameters, row):
#         p.val = v

#     absflux = calc_energy_flux(id=id, lo=2, hi=8, model=convmodel)
#     # srcnh.val = 20
#     # unabsflux = calc_energy_flux(id=id, lo=2/(1+z), hi=10/(1+z), model=torus)
#     r.append([absflux] + list(row))

# print("saving distribution plot data")
# r = numpy.asarray(r)
# assert len(rows) == len(r)
# numpy.savetxt(prefix + "intrinsic_photonflux.dist.gz", r)

# set_full_model(id, get_response(id)(model) + bkg_model * get_bkg_scale(id))
solver.set_best_fit()

#import sys; sys.exit()
#exit()
