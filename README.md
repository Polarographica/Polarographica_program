Polarographica_program

(NOW AVAILABLE ALSO FOR PYTHON_3_2020.03.21)

Cyclic voltammetry (CV), Linear-sweep voltammetry (LSV) and electrochemical impedance spectroscopy (EIS) are the prevalently
used techniques in electrochemical investigations. However, quantitative interpretation of experimentally acquired data is usually a non-straight forward task. Polarographica is a graphical user interface program, written in Python 2.7
programming language, for simulating and evaluating electroanalytical experiments. Polarographica supports simulation and data evaluation of the following electroanalytical techniques at systems involving macroporous electrodes as well as non porous macro- or microelectrodes:

1) Voltamerometric techniques:
Cyclic- and linear-sweep voltammetry (CV/LSV)
Fourier transform alternating current cyclic- and linear-sweep voltammetry (FT-ACCV/ FT-ACLSV)
Cyclic- and linear-sweep staircase voltammetry (STCV/STLSV)
Large sine amplitude cyclic voltammetry (LSA-CV)
Single- and multistep Chronoamperometry (CA)

2) Electrochemical Impedance spectroscopy techniques:
Potentiostatic electrochemical impedance spectroscopy (PEIS)
Distribution of relaxation times analysis of EIS data (DRT)

3) Largely automatized classical electrochemical analysis methods like:
Koutecky-Levich analysis
Randles-Sevcik analysis (reversible and irreversible)
Tafel-analysis
Cottrell-analysis

The user is free in choosing electrochemical reaction sequences like CE, EC, CEC mechanisms, diffusion domains (semi-infinite, or finite transmissive/reflective) and electrode geometries (planar, cylindrical, hollow-cylindrical, spherical, hollow-spherical). Furthermore, parasitic reaction can also be considered during data simulation and/or evaluation.


Features of Polarographica

Voltamperimetric techniques:
All voltamperimetric simulations of Polarographica are based on Laplace transformation techniques. The inverse Laplace transformation step for obtaining the time dependent surface concentrations of the electrochemically active 
species is performed numerically, by using the modified Tabbot contour suggested by [1] and its python implementation [2]

EIS
The mathematics of the EIS simulation/evaluation tool are taken from [3].
The distribution of relaxation times (DRT) analysis function (DRT-Tools-NNLS-DRT) of EIS data is basically a translation of the DRT-tools software [4] into Python code. However, it utilizes the NNLS algorithm instead of the quadprog algorithm for the final data fitting step.

Polarographica was created by 

Tim Tichter:  t.tichter@fu-berlin.de 
and
Jonathan Schneider: jonathan.schneider@fu-berlin.de

When using Polarographica please cite our publication : DOI: 10.1002/cphc.201901113

Date:  2020.01.08

References
[1]L.N.Trefethen, J.A.C.Weideman, and T.Schmelzer. Talbot quadraturesand rational approximations. BIT. Numerical Mathematics, 46(3):653 670, 2006.

[2] F.Nieuwveldt, Numerical Inversion of the Laplace Transform using the Talbot method. (Python recipe), 2009, http://code.activestate.com/recipes/576934-numerical-inversion-of-the-laplace-transform-using/?in=user-4172088

[3] A. Lasia, Electrochemical Impedance Spectroscopy and its Applications, Springer-Verlag New York, 2014, 10.1007/978-1-4614-8933-7

[4] T. H. Wan, M. Saccoccio, C. Chen, F. Ciucci, Electrochimica Acta 2015, 184, 483.

