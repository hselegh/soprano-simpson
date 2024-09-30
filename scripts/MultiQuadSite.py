#!/usr/bin/env python
# coding: utf-8
# Basic imports
import os

from ase import io as ase_io

from soprano.calculate.nmr.simpson_multi_quad import *

from soprano.selection import AtomSelection
import matplotlib.pyplot as plt

import argparse

# Add arguments to the run the script as command line.
parser = argparse.ArgumentParser()
parser.add_argument('-m', dest='magres_file', type=str, default="input.magres",
                    help='Magres file name/path. Default: "input.magres".')
parser.add_argument('-e', dest='element', type=str,
                    help='Element symbol eg. Na, C, O,... No default value, script will not run without it.')
parser.add_argument('-i', dest='isotope', type=int,
                    help='Isotope number, write as an integer eg. 1, 2, 3,... No default value, script will not run without it.')
parser.add_argument('-b', dest='field', type=float, default=400.0,
                    help='Magnet field strength in 1H Larmor Frequency and MHz.')
parser.add_argument('-c', dest='cs_iso', type=bool, default=False,
                    help='Chemical shift treaded as isotropic, ie. no CSA.')
parser.add_argument('-d', dest='det', type=str, default="I1p",
                    help='Can be either "I1p" (I1+ detection operator) or "I1c" (central transition) detection. Integer spins do not have a central transition, so the detection operator must be "I1p". Default: "I1p".')
parser.add_argument('-np', dest='np', type=str, default="1024",
                    help='Number of points in the simulated FID. Default: "1024".')
parser.add_argument('-sw', dest='sw', type=str, default="40000",
                    help='Spectral window, increase or decrease this value depending on the width of the signals (be careful with folding or really broad signals). Default: "40000".')
parser.add_argument('-mas', dest='mas', type=str, default="10000", help='MAS rate. Default: "10000".')
parser.add_argument('-cr', dest='crystal_file', type=str, default="rep168",
                    help='Crystal file for sampling powder orientation. Default: "rep168".')
parser.add_argument('-g', dest='gamma', type=str, default="40",
                    help='Number of gamma angles, for MAS use between 30-50 and for static use 1. Default: "40".')
parser.add_argument('-gQ', dest='gradQ', type=float, default="1.0",
                    help='Gradient for Cq, used in case of over or underestimation. This value will DIVIDE all the Cq values in the .magres file. Default value: "1.0".')
parser.add_argument('-lb', dest='LB', type=str, default="100", help='Line broadening (in Hz). Default value: "100".')
parser.add_argument('-zf', dest='ZeroFill', type=str, default="8192", help='Zero filling. Default value: "8192".')
parser.add_argument('-xy', dest='outXY', type=str, default="simpson.xy",
                    help='Name and path for the .xy output. Default: "simpson.xy".')
parser.add_argument('-p', dest='plot', type=bool, default=False,
                    help='Choose if plot or not the spectrum using pyplot at the end of the execution, write as True or False. Default: "False".')
parser.add_argument('-core', dest='cores', type=str, default="8",
                    help='Define the number of cores that Simpson will use. Default: "8".')
parser.add_argument('-sodsite', dest='sodsite', type=bool, default=False,
                    help='Choose if write the files XSPEC and SPECTRA for each site for use in the SOD package, write as True or False. Default: "False".')
parser.add_argument('-sodsum', dest='sodsum', type=bool, default=False,
                    help='Choose if write the files XSPEC and SPECTRA for the spectra sum for use in the SOD package, write as True or False. Default: "False".')
parser.add_argument('-sum', dest='sum', type=bool, default=True,
                    help='Choose if sums the .xy files generated for each site creating a .xy file, write as True or False. Default: "True".')
args = parser.parse_args()

# Read .magres file - you can change "args.magres_file" to a default file name using 'filename'.
magres = ase_io.read(args.magres_file)

# Soprano references the chemical shift as ms_shifts = references + gradients * ms_shieldings [y(exp)=a.x(dft)+b].
# References in ppm
references = {'H': 30,
              'C': 175,
              'O': 300,
              'Na': 472,
              'Cs': 5184.1}

# Gradients for the referencing
gradients = {'H': -0.98,
             'C': -1.00,
             'O': -1.00,
             'Na': -0.845,
             'Cs': -0.8869}

# Element selection for the NMR
# You can change "args.element" to a default element using 'element' (string - example 'Na')
# The Multi script cannot use more than one element/isotope at once
labels = magres.get_array('labels')
indices = magres.get_array('indices')
jmol_labels = ["{0}_{1}".format(l, i) for l, i in zip(labels, indices)]
magres_sel = AtomSelection.from_element(magres, args.element).subset(magres)

nsites = len(magres_sel.get_chemical_symbols())

xsite: int
for xsite in range(1,nsites+1):
    xsi = str(args.element + '.' + str(xsite))
    magres_atom = AtomSelection.from_selection_string(magres, xsi).subset(magres)

    # Write the spin system section for the Simpson loop and stores the list of parameters to be written in the
    # write_input_multi_quad command.
    #
    # magres_sel = variable with the selected element information stored on it, wrote in the last section of this script
    #
    # isotope(int) = variable used to set the isotope. Must be an integer.
    #
    # use_ms(bool) = Use magnetic shield - set as True. Script will not run with it as False. Hidden below.
    #
    # field(float) = magnetic field/1e6 for the ppm to frequency (because of SIMPSON/TCL - just accept it) to ppm
    # conversion.
    # MAKE SURE THAT THE FIELD HERE IS THE SAME AS THE PROTON_FREQUENCY IN THE SIMPSON PARAMETERS SECTION.
    #
    # atomsys(str) = element (symbol) used (also used in the ppm-frequency-ppm conversion).
    #
    # ms_iso(bool) = if True, chemical shift is made isotropic.
    #
    # q_order(int) = 2, do not change.
    #
    # dip_sel = None, do not change.
    #
    # path(str) = simpson.spinsys, do not change.
    #
    # ref = references, reading from the list above.
    #
    # grad = gradients, reading from the list above.
    #
    # gradQ(float) = gradient for Cq. Value that will DIVIDE the Cq values in the .magres file.

    write_spinsys_multi_quad(magres_atom, isotope=args.isotope, field=args.field, atomsys=args.element, ms_iso=args.cs_iso,
                             q_order=2, path='simpson_{0}.spinsys'.format(str(xsite)), ref=references, grad=gradients, gradQ=args.gradQ)

    # Add the spin system file to the SIMPSON input.
    simp = SimpsonSequenceMultiQuad('simpson_{0}.spinsys'.format(str(xsite)))
    h_freq = str(args.field * 1000000.0)

    # To edit the parameters section in the SIMPSON input.
    simp.pars = {
        "proton_frequency": h_freq,  # Spectrometer 1H frequency. This value must be equal to the "field" variable above.
        "method": "direct",  # Do not change.
        "start_operator": "I1x",  # Starting operator for an ideal excitation.
        "detect_operator": args.det,
        # Can be either "I1p" (I1+ detection operator) or "I1c" (central transition) detection.
        # Integer spins does not have central transitions,
        # so detection operator must be "I1p".
        "np": args.np,  # Number of points in the simulated FID.
        "sw": args.sw,  # Spectral window, increase or decrease this value depending on
        # the width of the signals (be careful with folding or really broad signals).
        "spin_rate": args.mas,  # MAS rate.
        "num_cores": args.cores,  # Uncomment if need to set a fixed number of cores. It will run in all cores as default.
        "crystal_file": args.crystal_file,  # Crystal file for sampling orientation,
        # use either rep168 or rep320 for good MAS lineshape simulation. Rep678 or rep2000 may
        # be needed depending on the Cq.
        "verbose": "111",  # Do not change.
        "gamma_angles": args.gamma,  # Number of gamma angles, for MAS use between 30-50.
        "variable": "tdwell     1.0e6/sw"  # Do not change.
    }

    # To edit the pulse section in the SIMPSON input. Default gives a perfect excitation (I1x) - no pulse simulation.
    PulseSeq = {
        "pars": {},
        "pulseq": """
       global par
    
        reset
        acq
        for {set i 1} {$i < $par(np)} {incr i} {
             delay $par(tdwell)
             acq
            }
    """,
        "main": """
    """,
    }

    simp.apply_template(PulseSeq)

    # Write the SIMPSON input file as 'simpson.in'.
    #
    # LB = line broadening.
    #
    # ZeroFill = zero filling.
    simp.write_input_multi_quad('simpson_{0}.in'.format(str(xsite)), args.LB, args.ZeroFill)

    # Run SIMPSON, you can change 'simpson' to your SIMPSON path.
    simpsonexec = str('simpson simpson_' + str(xsite) +'.in' + ' > simpson_' + str(xsite) +'.log')
    os.system(simpsonexec)

    # Read the FT transformed SIMPSON output (simpson.spe), convert the x-axis from Hz to ppm
    # (using *1/[field*(gamma_X/gamma_1H)]).
    # Delete the imaginary part of the file (third column) and write as a .xy file in ppm.
    simpson_xy_multi_quad(f'simpson_{str(xsite)}.spe', args.outXY + '_' + str(xsite) + '.xy')

    # Plots the spectrum - good for quickly checking the simulation.
    if args.plot is False:
        print("Job done!")
    else:
        data = np.loadtxt(args.outXY + '_' + str(xsite) + '.xy')
        x = data[:, 0]
        y = data[:, 1]
        plt.plot(x, y, 'r')
        plt.show()

    # Creates the files XSPEC and SPECTRA for each site to be used in SOD (https://github.com/gcmt-group/sod).
    if args.sodsite is True:
        data_sod = np.loadtxt(args.outXY + '_' + str(xsite) + '.xy')
        x_sod = data_sod[:, 0]
        y_sod = data_sod[:, 1]
        np.savetxt('XSPEC', x_sod)
        y_sodT = y_sod[None, :]
        np.savetxt('ysod', y_sodT)
        with open('ysod', 'r+') as ysodt:
            with open('SPECTRA_'+ str(xsite) , 'a+') as f:
                yst = ysodt.read()
                brk = "\n"
                spec = args.ZeroFill + brk + yst
                f.write(spec)
                f.close()
            ysodt.close()
        os.remove('ysod')
    else:
        pass


# Sums all the sites spectra and generates a .xy file
if args.sum is True:
    ymasterlist = []
    for xsite in range(1, nsites + 1):
        globals()[f'fromxy_{xsite}'] = np.loadtxt(args.outXY + '_' + str(xsite) + '.xy')
        globals()[f'yfromxy_{xsite}'] = globals()[f'fromxy_{xsite}'][:, 1]
        ymasterlist.append(globals()[f'yfromxy_{xsite}'])

    ysum = [0] * len(yfromxy_1)
    sumx = fromxy_1[:, 0]
    for lst in ymasterlist:
        for i in range(len(lst)):
            ysum[i] += lst[i]

    ysum_np = np.array(ysum)
    sum_columns = np.column_stack([sumx, ysum_np])
    np.savetxt(args.outXY + '_sum.xy', sum_columns)

else:
    pass

# Creates the files XSPEC and SPECTRA for the spectra sum to be used in SOD (https://github.com/gcmt-group/sod).
if args.sodsite is True:
    data_sod_sum = np.loadtxt(args.outXY + '_sum.xy')
    x_sod_sum = data_sod_sum[:, 0]
    y_sod_sum = data_sod_sum[:, 1]
    np.savetxt('XSPEC', x_sod_sum)
    y_sodT_sum = y_sod_sum[None, :]
    np.savetxt('ysod_sum', y_sodT_sum)
    with open('ysod_sum', 'r+') as ysodt_sum:
        with open('SPECTRA' , 'a+') as f:
            yst_sum = ysodt_sum.read()
            brk = "\n"
            spec_sum = args.ZeroFill + brk + yst_sum
            f.write(spec_sum)
            f.close()
        ysodt_sum.close()
    os.remove('ysod_sum')
else:
    pass