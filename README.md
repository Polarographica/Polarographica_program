$\bf{Polarographica\hspace{2mm}program}$

Polarographica is an open-source program for simulating and evaluating electroanalytical experiments. Initially designed for the simulation of cyclic voltammetry at marcoporous electrodes, Polarographica is steadily expanding. 

-----------------------

Polarographica was initially created by 

Tim Tichter: timtic@gmx.de 

Jonathan Schneider: js-jena@gmx.de

-----------------------

Scientists/colleagues who have contributed to Polarographica are gratefully achnowledged for their indispensible input.

Marius Gernhard: marius.gernhard@uni-bayreuth.de $\hspace{5mm}$ 

--> $\hspace{5mm}$ Significant contribution in the serial-communication interface of PolArStat

-----------------------

$\bf{History\hspace{2mm}of\hspace{2mm}Polarographica}$

Since Polarographica is steadily expanding/improving, there are already a few releases. The following version are available (starting with the most recent one):

$\bf{Polarographica\hspace{2mm}2.2.2.4}$:
The sub-version Polarographica 2.2.2.4 is an improvement of the sub-verion Polarographica 2.2.2.2. It contains a new radial basis function for computing the distribution of relaxation times function from electrochemical impedance spectroscopy data - the Havriliak-Negami-DRT. This function accounts for an asymmetrical shape of the time-constant distribution. Essentially, it includes the Debye-relaxation (classical "spike-DRT"-case), the Cole-Cole relaxation (constant phase element) and the Cole-Davidson relaxation (e.g. relaxation behaviour of EIS at an ideal pore) as special cases and is therefore the most general type of DRT in Polarographica. An additional (yet minor) modification was a bug-fix in all DRT-modules. Please note that until version 2.2.2.2, the annotation of the ordinate in the DRT was not correct. All DRTs were seemingly divided by a $\Delta$ ln($\tau$) term. The "seemingly" refers to the fact that the actual computation of the DRTs was/is correct, however, the annotation was not. The last modification belogs to all impedance modules (referring to simulation, re-loading and fitting). Before version 2.2.2.4, EIS data was not displayed in an equal aspect ratio. This somewhat "bad practice" has been removed. 

An additional feature (not limited to Polarographica 2.2.2.4) is the addition of a "hands-on" guide how to assemble the PolArStat-potentiostat on a breadboard. Owing to the overwhelming amount of requests regarding PolArStat-shields, we have decided to create a detailed guide. This is addressed to all experimentalists who want to check out PolArStat without the need of ordering/manufacturing a PCB-shield and without the need and/or experience of soldering. This guide is included in the version-release of Polarographica 2.2.2.4. Additionally, you can access it as a separate folder.

$\bf{Polarographica\hspace{2mm}2.2.2.2}$:
The Sub-version Polarographica 2.2.2.2 is an improvement of the native Polarographica 2.2.2 version, where the PolArStat modules were thoroughly re-designed. These changes are basically an improvement of the data-acquisition backend, which is now running on a thread and is therefore independent of the main GUI. This removes the flickering of the data acquisition monitor. Note, that the zipped folder of Polarographica 2.2.2.2 does contain the firmware for the PolArStat which can be uploaded to a classical Arduino. The cool feature: when assembled with the shield described in our Publication [5], you have a fully functional potentiostat for about 40 €, which can run CV and CA. Design files are found in the supplementary information of [5].

$\bf{Polarographica\hspace{2mm}2.2.2}$:
The version Polarographica 2.2.2 and its derivatives are running under Python 3.9. In version 2.2.2, the Distribution of relaxation times module was thoroughly re-designed. It contains a new radial basis function, which is based on the analytical DRT of the Cole-Cole-Kernel. Furthermore, it contains additional features such as Kramers-Kronig check.

$\bf{Polarographica\hspace{2mm}1.0.0-2.1.0}$:
Polarographica-Versions 1.0.0 to 2.0.0 have been removed! Polarographica 2.1.0 is an improved version (running under python 3.7). A Bug in Voltammetric modules has been fixed and the GUI has been renewed. 



-----------------------

$\bf{Purpose\hspace{2mm}of\hspace{2mm}Polarographica\}$:

Cyclic voltammetry (CV), Linear-sweep voltammetry (LSV) and electrochemical impedance spectroscopy (EIS) are the prevalently
used techniques in electrochemical investigations. However, quantitative interpretation of experimentally acquired data is usually a non-straight forward task. Polarographica is a graphical user interface program for simulating and evaluating electroanalytical experiments. Polarographica supports simulation and data evaluation of the following electroanalytical techniques at systems involving macroporous electrodes as well as non porous macro- or microelectrodes:

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

-----------------------

$\bf{Features\hspace{2mm}of\hspace{2mm}Polarographica\}$:

Voltamperimetric techniques:
Most voltamperimetric simulations of Polarographica are based on Laplace transformation techniques. The inverse Laplace transformation step for obtaining the time dependent surface concentrations of the electrochemically active 
species is performed numerically, by using the modified Talbot contour suggested by [1] and its python implementation [2]

EIS
The mathematics of the EIS simulation/evaluation tool are taken from [3].
The distribution of relaxation times (DRT) analysis function (DRT-Tools-NNLS-DRT) of EIS data is basically a translation of the DRT-tools software [4] into Python code. However, it utilizes the NNLS algorithm instead of the quadprog algorithm for the final data fitting step. From version >=2.2.2, an additional DRT function was included, which uses the Cole-Cole function as radial basis function and which is indenpendent from DRT-tools. Furthermore, the entire DRT-module was re-designed.


-----------------------

Date:  2024.03.23

References:

[1] L.N.Trefethen, J.A.C.Weideman, and T.Schmelzer. Talbot quadraturesand rational approximations. BIT. Numerical Mathematics, 46(3):653 670, 2006.

[2] F.Nieuwveldt, Numerical Inversion of the Laplace Transform using the Talbot method. (Python recipe), 2009, http://code.activestate.com/recipes/576934-numerical-inversion-of-the-laplace-transform-using/?in=user-4172088

[3] A. Lasia, Electrochemical Impedance Spectroscopy and its Applications, Springer-Verlag New York, 2014, 10.1007/978-1-4614-8933-7

[4] T. H. Wan, M. Saccoccio, C. Chen, F. Ciucci, Electrochimica Acta 2015, 184, 483.

[5] T. Ticher, M. Gernhard, P.C.K. Vesborg, PolArStat: An Arduino based potentiostat for low-power electrochemical applications, Electrochimica Acta, 2023, 143119, 10.1016/j.electacta.2023.143119




