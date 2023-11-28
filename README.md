# soprano* 
> Working version for the implementation of the interface with Simpson.

The original version can be found in https://ccp-nc.github.io/soprano/intro.html and must be cited accordingly to https://ccp-nc.github.io/soprano/citing.html.



## Installation
1. Install the newest version of Python3 from www.python.org.
2. Extract the Soprano folder, in the Terminal and inside the Soprano folder, run the command: `python3.xx setup.py install`.

## Using it

For spin ½, use the script `MultiHalf.py`; for spin > ½, use the script `MultiQuad.py`.

The `MultiQuad.py` (and similarly the `MultiHalf.py` – with the only difference is not having the `-gQ` argument) can be used with the following structure: 

```
usage: MultiQuad.py [-h] [-m MAGRES_FILE] [-e ELEMENT] [-i ISOTOPE] [-b FIELD] [-c CS_ISO] [-d DET] [-np NP] [-sw SW] [-mas MAS] [-cr CRYSTAL_FILE] [-g GAMMA] [-gQ GRADQ] [-lb LB] [-zf ZEROFILL] [-xy OUTXY] [-p PLOT]
```
options:
```
-h	show this help message and exit
-m	Magres file name/path. Default: "input.magres".
-e	Element symbol eg. Na, C, O,... No default value, script will not run without it.
-i	Isotope number, write as an integer eg. 1, 2, 3,... No default value, script will not run without it.
-b	Magnet field strength in 1H Larmor Frequency and MHz.
-c	Chemical shift treaded as isotropic, ie. no CSA.
-d	Can be either "I1p" (I1+ detection operator) or "I1c" (central transition) detection. Integer spins do not have a central transition, so the detection operator must be "I1p". Default: "I1p".
-np	Number of points in the simulated FID. Default: "1024".
-sw	Spectral window, increase or decrease this value depending on the width of the signals (be careful with folding or really broad signals). Default: "40000".
-mas	MAS rate. Default: "10000".
-cr	Crystal file for sampling powder orientation. Default: "rep168".
-g	Number of gamma angles, for MAS use between 30-50 and for static use 1. Default: "40".
-gQ	Gradient for Cq, used in case of over or underestimation. This value will DIVIDE all the Cq values in the .magres file. Default value: "1.0".
-lb	Line broadening (in Hz). Default value: "100".
-zf	Zero filling. Default value: "8192".
-xy	Name and path for the .xy output. Default: "simpson.xy".
-p	Choose if plot or not the spectrum using pyplot at the end of the execution, write as True or False. Default: "False".
-core   Define the number of cores that Simpson will use. Default: "8".
-sod    Choose if write the files XSPEC and SPECTRA for use in the SOD package (https://github.com/gcmt-group/sod), write as True or False. Default: "False".
```

An example can be:
`python3.11 MultiQuad.py -m input.magres -e Na -i 23 -gQ 1.00 -lb 300 -zf 4096 -xy output.xy -p True`

The only thing you will need to edit in the MultiQuad.py file is the references and gradients lists for the isotopes in your system, there is no need for the list to have all the isotopes in your magres file, just the ones you want to simulate with Simpson.

``` 
# Soprano references the chemical shift as ms_shifts = references + gradients * ms_shieldings [y(exp)=a.x(dft)+b].
# List of references in ppm
references = {'H': 30,
              'C': 175,
              'O': 300,
              'Na': 500}

# List of gradients for the referencing
gradients = {'H': -0.98,
             'C': -1.00,
             'O': -1.00,
             'Na': -1.00}
```

If you need to change the path for your Simpson executable, you need to change the following section in the script.
```
# Run SIMPSON, you can change 'simpson' to your SIMPSON path.
os.system("simpson simpson.in > simpson.log")
```
