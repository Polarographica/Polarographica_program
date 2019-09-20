# Polarographica_program

#Cyclic voltammetry (CV), Linear-sweep voltammetry (LSV) and electrochemical impedance spectroscopy (EIS) are the predominant
#techniques used in electrochemical investigations.
#Polarographica is a graphical user interface program written in Python 2.7
#for simulating and evaluating CV, LSV and EIS. 


#Polarographica was created by 
#       Tim Tichter 
#email: t.tichter@fu-berlin.de 
#and
#       Jonathan Schneider
#email: jonathan.schneider@fu-berlin.de

#Date:  2019.09.20

#Voltamperimetric techniques:
#Voltamperimetric techniques simulations of Polarographica are based on Laplace transformation techniques.
#The inverse Laplace transformation step for obtaining the time dependent surface concentrations of the electrochemically active 
#species involved in the reaction, required for the CV simulations are obtained via numerical inverse Laplace transform using 
#the modified Tabbot contour suggested by [1] and its python implementation [2]


#EIS
#The mathematics of the EIS simulation/evaluation tool are extracted from [3].
#Polarographica also provides distribution of relaxation times (DRT) analysis for EIS data. 
#The DRT-Tools-NNLS-DRT function is basically a translation of the DRT-tools software [4] into a Python GUI but 
#uses the NNLS algorithm instead of quadprog algorithm for fitting


#References
#[1]L.N.Trefethen, J.A.C.Weideman, and T.Schmelzer. Talbot quadraturesand rational approximations. BIT. Numerical Mathematics, 46(3):653 670, 2006.

#[2] F.Nieuwveldt, Numerical Inversion of the Laplace Transform using the Talbot method. (Python recipe), 2009, http://code.activestate.com/recipes/576934-numerical-inversion-of-the-laplace-transform-using/?in=user-4172088

#[3] A. Lasia, Electrochemical Impedance Spectroscopy and its Applications, Springer-Verlag New York, 2014, 10.1007/978-1-4614-8933-7

#[4] T. H. Wan, M. Saccoccio, C. Chen, F. Ciucci, Electrochimica Acta 2015, 184, 483.

