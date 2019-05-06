
# coding: utf-8


#======================================================================================================================================================================
#Cyclic voltammetry (CV), Linear-sweep voltammetry (LSV)and electrochemical impedance spectroscopy (EIS) are the predominant
#techniques used in electrochemical investigations.
#Polarographica is a graphical user interface program written in Python 2.7
#for simulating and evaluating CV, LSV and EIS. 

#-------------------------------
#Polarographica was created by 
#       Tim Tichter 
#email: t.tichter@fu-berlin.de 
#and
#       Jonathan Schneider
#email: jonathan.schneider@fu-berlin.de

#Date:  2019.05.06

#-------------------------------
#CV/LSV
#CV simulations of Polarographica are based on Laplace transformation techniques.
#The inverse Laplace transformation step for obtaining the time dependent surface concentrations of the electrochemically active 
#species involved in the reaction, required for the CV simulations are obtained via numerical inverse Laplace transform using 
#the modified Tabbot contour suggested by [1] and its python implementation [2]

#-------------------------------
#EIS
#The mathematics of the EIS simulation/evaluation tool are extracted from [3].
#Polarographica also provides distribution of relaxation times (DRT) analysis for EIS data. 
#The DRT-Tools-NNLS-DRT function is basically a translation of the DRT-tools software [4] into a Python GUI but 
#uses the NNLS algorithm instead of quadprog algorithm for fitting


#--------------------------------
#References
#[1]L.N.Trefethen, J.A.C.Weideman, and T.Schmelzer. Talbot quadratures
#   and rational approximations. BIT. Numerical Mathematics,
#   46(3):653 670, 2006.
#[2] F.Nieuwveldt, Numerical Inversion of the Laplace Transform using the Talbot method. (Python recipe), 2009,  
#    http://code.activestate.com/recipes/576934-numerical-inversion-of-the-laplace-transform-using/?in=user-4172088
#[3] A. Lasia, Electrochemical Impedance Spectroscopy and its Applications, Springer-Verlag New York, 2014, 10.1007/978-1-4614-8933-7
#[4] T. H. Wan, M. Saccoccio, C. Chen, F. Ciucci, Electrochimica Acta 2015, 184, 483.

#====================================================================================================================================================================


#Polarographica_program 

# In[1]:


from Tkinter import *

from tkFileDialog   import askopenfilename
from tkFileDialog import asksaveasfilename
import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit
from scipy.optimize import nnls
from scipy.linalg import toeplitz
from scipy.special import kv, iv, gamma
from scipy.integrate import quad

from cmath import *
from math import ceil,floor
import mpmath as mp
mp.dps = 25; mp.pretty = True

global F
global R
global File_Was_Loaded


F = 96485.0
R = 8.314
File_Was_Loaded = 0

def cot(phi):
    return 1.0/tan(phi)

def csc(phi):
    return 1.0/sin(phi)

def coth(x):
    return 1/tanh(x)



# In[2]:


def BaseGetter():
    global BaseData
    Delim = Delimiter
    
    if Delim == 1:
       
        BasePath = askopenfilename()
        BaseData = np.genfromtxt(BasePath, delimiter='\t')

        
    if Delim == 0:
       
        BasePath = askopenfilename()
        BaseData = np.genfromtxt(BasePath, delimiter='')
        
        
    global BasePotenzial
    global BaseStrom
    
    
    if FirstComesxx  == 1:
        BasePotenzial     =  BaseData[FromRowxx:ToRowxx:Readeveryxx,0::2]  * UmrPot              
        BaseStrom         =  BaseData[FromRowxx:ToRowxx:Readeveryxx,1::2]  * UmRStrom

    if FirstComesxx  == 0:
        BasePotenzial     =  BaseData[FromRowxx:ToRowxx:Readeveryxx,1::2]  * UmrPot
        BaseStrom         =  BaseData[FromRowxx:ToRowxx:Readeveryxx,0::2]  * UmRStrom
        
    if BasePotenzial[0,0] > BasePotenzial[1,0]:
        BasePotenzial     = BasePotenzial[::-1]
        BaseStrom         = BaseStrom[::-1]


# In[3]:


def OpenFile():
    root = Toplevel()
    root.title("Your Data")
    if CV == 1 or CA == 1:
        root.geometry("500x500")     
    if IMP == 1:
        root.geometry("1000x500")     
    
    Delim = Delimiter
    
    global data
    
    global File_Was_Loaded
    File_Was_Loaded = 1
    
    if Delim == 1:
       
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='\t')

        
    if Delim == 0:
       
        Path = askopenfilename()
        data = np.genfromtxt(Path, delimiter='')
    
    
    
    
    #============================================================================
    #Plotting of loaded file
    #============================================================================
    
    #============================================================================
    #define entry window
    #============================================================================
    
    if CA == 1 or CV == 1:
        f = Figure(figsize=(5, 5), dpi=100)
        b = f.add_subplot(111)
        
    if IMP == 1:
        f  = Figure(figsize=(5, 10), dpi=100)        
        b1 = f.add_subplot(121)
        b2 = f.add_subplot(122)      
        
        
    #============================================================================
    
    global Potenzial
    global Strom
    global Frequency_Array
    
    if CV == 1:
        if FirstComesxx  == 1:
            Potenzial     =  data[FromRowxx:ToRowxx:Readeveryxx,0::2] * UmrPot            
            Strom         =  data[FromRowxx:ToRowxx:Readeveryxx,1::2] * UmRStrom
        if FirstComesxx  == 0:
            Potenzial     =  data[FromRowxx:ToRowxx:Readeveryxx,1::2] * UmrPot
            Strom         =  data[FromRowxx:ToRowxx:Readeveryxx,0::2] * UmRStrom
            
    if CA == 1:
        if FirstComesxx  == 1:
            Potenzial     =  data[FromRowxx:ToRowxx:Readeveryxx,0::2] * UmrPot            
            Strom         =  data[FromRowxx:ToRowxx:Readeveryxx,1::2] * UmRStrom
        if FirstComesxx  == 0:
            Potenzial     =  data[FromRowxx:ToRowxx:Readeveryxx,1::2] * UmrPot
            Strom         =  data[FromRowxx:ToRowxx:Readeveryxx,0::2] * UmRStrom
            
    #==================================================================================
    if IMP == 1:
        Frequency         =  data[FromRowxx:ToRowxx:Readeveryxx,0:1:1] 
        Potenzial         =  data[FromRowxx:ToRowxx:Readeveryxx,1:2:1] * UmrPot            
        Strom             =  data[FromRowxx:ToRowxx:Readeveryxx,2:3:1] * UmRStrom
        
        #==================================================
        #turn frequencies alsways from high to low
        
        if Frequency[1] > Frequency[0]:            
            Frequency         = Frequency[::-1]
            Potenzial         = Potenzial[::-1]             
            Strom             = Strom[::-1]
        #==================================================
    #==================================================================================                  
    
    if IMP != 1:
        if Potenzial[0,0] > Potenzial[1,0]:
            Potenzial     = Potenzial[::-1]
            Strom         = Strom[::-1]
            
            
    global Weite
    Weite = Potenzial.shape[0]
    global NumMeas
    NumMeas = Potenzial.shape[1]
    
    for i in range(Weite):
        Stromarrays             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))                          
        Potenzialarrays         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:])) 
        LenPotArrays            = len(Potenzialarrays)
        if IMP == 1:
            Frequency_Array     = Frequency  
        
        if CV == 1:
            b.plot(Potenzial, 0.001*Strom, linestyle='-',marker='',color='k')
            b.set_xlabel('E vs. Ref. / V', fontsize=12)
            b.set_ylabel('I / mA', fontsize=12)
            for axis in ['top','bottom','left','right']:
                b.spines[axis].set_linewidth(2)
                b.spines[axis].set_color('k')
            b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
        if CA == 1:
            b.plot(Potenzial, 0.001*Strom, linestyle='-',marker='',color='k')
            b.set_xlabel('t / s', fontsize=12)
            b.set_ylabel('I / mA', fontsize=12)
            for axis in ['top','bottom','left','right']:
                b.spines[axis].set_linewidth(2)
                b.spines[axis].set_color('k')
            b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
    #=========================================================
    #out of for-loop for impedance
    #=========================================================
    
    if IMP == 1:
        b1.set_title("Nyquist-Plot")
        b1.plot(Potenzial, -Strom, linestyle='-',marker='',color='k')
        b1.set_xlabel('Z_real / Ohm', fontsize=12)
        b1.set_ylabel('-Z_imag / Ohm', fontsize=12)
        for axis in ['top','bottom','left','right']:
            b1.spines[axis].set_linewidth(2)
            b1.spines[axis].set_color('k')
        b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            

        b2.set_title("Bode-Plot")    
        b2.plot(np.log10(Frequency), -Strom, linestyle='-',marker='',color='b', label = '-Z.imag')
        b2.plot(np.log10(Frequency), Potenzial, linestyle='-',marker='',color='r',label = 'Z.real')
        b2.set_xlabel('log10(freq /s^-1)', fontsize=12)
        b2.set_ylabel('-Z_imag and Z_real / Ohm  ', fontsize=12)
        b2.plot()
        b2.legend()

        for axis in ['top','bottom','left','right']:
            b2.spines[axis].set_linewidth(2)
            b2.spines[axis].set_color('k')
        b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
    #=========================================================    
    
         
    canvas = FigureCanvasTkAgg(f, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)


# In[4]:


def LSV_DATA():
    global CV
    CV = 1
    global CA
    CA = 0
    global IMP
    IMP = 0
    GetData()
    
def CA_DATA():
    global CV
    CV = 0
    global CA
    CA = 1
    global IMP
    IMP = 0
    GetData()

def IMP_DATA():
    global CV
    CV = 0
    global CA
    CA = 0
    global IMP
    IMP = 1
    GetData()


# In[5]:


def GetData():
    Fenster = Toplevel()                                                         
    Fenster.title("Get-Data")                         
    Fenster.geometry("400x270")                                            

    FromRowxxx_label      = Label(Fenster,text="Read from row")
    FromRowxxx_label.grid(row=0, column=0)
    FromRowxxx_Eingabe    = Entry(Fenster)                                               
    FromRowxxx_Eingabe.grid(row=0, column=1)

    ToRowxxx_label       = Label(Fenster,text="Read to row")
    ToRowxxx_label.grid(row=1, column=0)
    ToRowxxx_Eingabe     = Entry(Fenster)
    ToRowxxx_Eingabe.grid(row=1, column=1)
    
    Readeveryxxx_label   = Label(Fenster,text="Read every")
    Readeveryxxx_label.grid(row=2, column=0)
    Readeveryxxx_Eingabe = Entry(Fenster)
    Readeveryxxx_Eingabe.grid(row=2, column=1)
    
    if CV == 1:
        UmRStrom_Label = Label(Fenster,text="Current Factor to be Microampere")
        UmRStrom_Label.grid(row=3, column=0)
        UmRStrom_Eingabe = Entry(Fenster)
        UmRStrom_Eingabe.grid(row=3, column=1)
        
        UmrPot_Label = Label(Fenster,text="Potential Factor to be Volt")
        UmrPot_Label.grid(row=4, column=0)
        UmrPot_Eingabe = Entry(Fenster)
        UmrPot_Eingabe.grid(row=4, column=1)
    
    
    if CA == 1:
        UmRStrom_Label = Label(Fenster,text="Current Factor to be Microampere")
        UmRStrom_Label.grid(row=3, column=0)
        UmRStrom_Eingabe = Entry(Fenster)
        UmRStrom_Eingabe.grid(row=3, column=1)
    
        UmrPot_Label = Label(Fenster,text="Time Factor to be seconds")
        UmrPot_Label.grid(row=4, column=0)
        UmrPot_Eingabe = Entry(Fenster)
        UmrPot_Eingabe.grid(row=4, column=1)
    
    if IMP == 1:
        UmrPot_Label = Label(Fenster,text="Z_real Factor to be Ohm")
        UmrPot_Label.grid(row=3, column=0)
        UmrPot_Eingabe = Entry(Fenster)
        UmrPot_Eingabe.grid(row=3, column=1)
        
        UmRStrom_Label = Label(Fenster,text="Z_imag Factor to be Ohm")
        UmRStrom_Label.grid(row=4, column=0)
        UmRStrom_Eingabe = Entry(Fenster)
        UmRStrom_Eingabe.grid(row=4, column=1)
    
        
        
    if CV == 1:
        #Desired Reaction
        DesiReactxxx_Label  = Label(Fenster,text="Desired Reaction")
        DesiReactxxx_Label.grid(row=5, column=0)
        var1 = IntVar()
        Checkbutton(Fenster, text="Oxidation", variable=var1).grid(row=5,column=1, sticky=W)
        var2 = IntVar()
        Checkbutton(Fenster, text="Reduction", variable=var2).grid(row=5,column=2, sticky=W)
    if CA == 1:
        #Desired Reaction
        DesiReactxxx_Label  = Label(Fenster,text="Desired Reaction")
        DesiReactxxx_Label.grid(row=5, column=0)
        var1 = IntVar()
        Checkbutton(Fenster, text="Oxidation", variable=var1).grid(row=5,column=1, sticky=W)
        var2 = IntVar()
        Checkbutton(Fenster, text="Reduction", variable=var2).grid(row=5,column=2, sticky=W)

    
    if CV == 1:
        #First Comes
        FirstComesxxx_Label  = Label(Fenster,text="First comes")
        FirstComesxxx_Label.grid(row=6, column=0)
        var3 = IntVar()
        Checkbutton(Fenster, text="Potential", variable=var3).grid(row=6,column=1, sticky=W)
        var4 = IntVar()
        Checkbutton(Fenster, text="Current", variable=var4).grid(row=6,column=2, sticky=W)
    
    if CA == 1:
        #First Comes
        FirstComesxxx_Label  = Label(Fenster,text="First comes")
        FirstComesxxx_Label.grid(row=6, column=0)
        var3 = IntVar()
        Checkbutton(Fenster, text="Time", variable=var3).grid(row=6,column=1, sticky=W)
        var4 = IntVar()
        Checkbutton(Fenster, text="Current", variable=var4).grid(row=6,column=2, sticky=W)
    
    if IMP == 1:
        HasToBe_Label  = Label(Fenster,text="Data Order")
        HasToBe_Label.grid(row=6, column=0)
        HasToBe_Label1  = Label(Fenster,text="Freq___Z.real___Z.Imag")
        HasToBe_Label1.grid(row=6, column=1)
        
        
    

    
    Delimiterxxx_Label  = Label(Fenster,text="Delimiter")
    Delimiterxxx_Label.grid(row=7, column=0)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Tab", variable=var5).grid(row=7,column=1, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Space", variable=var6).grid(row=7,column=2, sticky=W)



    def Accept():
              
        global FromRowxx
        FromRowxx     = (int(FromRowxxx_Eingabe.get()))
        global ToRowxx
        ToRowxx       = (int(ToRowxxx_Eingabe.get()))
        global Readeveryxx
        Readeveryxx   = (int(Readeveryxxx_Eingabe.get()))
        global UmRStrom#
        UmRStrom      = (float(UmRStrom_Eingabe.get()))
        global UmrPot#
        UmrPot        = (float(UmrPot_Eingabe.get()))
        
        global DesiReactxx
        global FirstComesxx

        if CV == 1:
            DesiReactxx   = var1.get()
            FirstComesxx  = var3.get()
        if CA == 1:
            DesiReactxx   = var1.get()
            FirstComesxx  = var3.get()
        
        global Delimiter
        Delimiter         = var5.get()
    
    
    def Next():
        OpenFile()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    Accept = Button(Fenster, text="Accept",command=Accept)
    Accept.grid(row=8, column=0) 
    
    Next = Button(Fenster, text="Next",command=Next)
    Next.grid(row=9, column=0)  


# In[6]:


def KouteckyLevichWindow():
    
    
    Fenster = Toplevel()                                                         
    Fenster.title("Koutecky-Levich Curve-Getter")                         
    Fenster.geometry("400x400")
    
    OCP_Label = Label(Fenster,text="OCPs vs E \nRef in V")
    OCP_Label.grid(row=0, column=0)
    
    RotRates_Label = Label(Fenster,text="Rotation rates \n in rpm")
    RotRates_Label.grid(row=0, column=1)
    
    Concent_Label = Label(Fenster,text="Concentrations as fract.\n of the highest one")
    Concent_Label.grid(row=0, column=2)
    
    #IF ONLY ONE OCP
    #______________________________________
    var1 = IntVar()
    Checkbutton(Fenster, text="One OCP", variable=var1).grid(row=80, column=0, sticky=W)
    OneOCP_Eingabe = Entry(Fenster)                                               
    OneOCP_Eingabe.grid(row=81, column=0)
    
    
    #IF ONLY ONE Rotation rate
    #______________________________________
    var2 = IntVar()
    Checkbutton(Fenster, text="One Rotation", variable=var2).grid(row=80, column=2, sticky=W)
    OneRot_Eingabe = Entry(Fenster)                                               
    OneRot_Eingabe.grid(row=81, column=2)
    

    OCPs     = []
    RotRates = []
    Concents = []

    for i in range(NumMeas):
        
        en = Entry(Fenster)
        en.grid(row=i+1, column=0)
        OCPs.append(en)

        en1 = Entry(Fenster)
        en1.grid(row=i+1, column=1)
        RotRates.append(en1)
        
        en2 = Entry(Fenster)
        en2.grid(row=i+1, column=2)
        Concents.append(en2)
    
    
    
    def GetOCPs():    
        OnlyOneOCP = var1.get()

        OCPots = []
        
        if OnlyOneOCP ==0:
            for entry in OCPs:
                OCPots.append(float(entry.get()))
        
        if OnlyOneOCP ==1:
            for j in range(NumMeas):
                OCPots.append(float(OneOCP_Eingabe.get()))

        global OCPArray
        OCPArray = np.array(OCPots) 

    
    
    def GetRots():    
        Rots= []
        for entry in RotRates:
            Rots.append(float(entry.get()))
        global RotRatesArray
        RotRatesArray = np.array(Rots) 
        global VarParRot
        global VarParCon
        VarParRot = 1
        VarParCon = 0
        
        
    def GetConcents():    
        OnlyOneRot = var2.get()
        
        if OnlyOneRot ==1:
            global RotRate
            RotRate = (float(OneRot_Eingabe.get()))
        
        Concs= []
        for entry in Concents:
            Concs.append(float(entry.get()))
        global ConcentrationsArray
        ConcentrationsArray = np.array(Concs) 
        global VarParRot
        global VarParCon
        VarParRot = 0
        VarParCon = 1
        
   
        
    def Next():
        KL_NextLevel2()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    button=Button(Fenster,text="Accept OCPs",command=GetOCPs).grid(row=82,column=0)
    button=Button(Fenster,text="Accept RotRates",command=GetRots).grid(row=82,column=1)
    button=Button(Fenster,text="Accept Concentrations",command=GetConcents).grid(row=82,column=2)
    button=Button(Fenster,text="Next",command=Next).grid(row=83,column=1)

    
    
    


# In[7]:


def KL_NextLevel2():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Koutecky-Levich Analysis Param-Getter")                         
    Fenster.geometry("400x400")

    
    Getter_Label = Label(Fenster, text="Get From KL")
    Getter_Label.grid(row=0, column=0)
    var1 = IntVar()
    Checkbutton(Fenster, text="get n", variable=var1).grid(row=0, column=1, sticky=W)
    var2 = IntVar()
    Checkbutton(Fenster, text="get c", variable=var2).grid(row=0,column=2, sticky=W)
    var3 = IntVar()
    Checkbutton(Fenster, text="get D", variable=var3).grid(row=0,column=3, sticky=W)
    
        
    n_Label = Label(Fenster, text="n")
    n_Label.grid(row=1, column=0)
    n_Eingabe = Entry(Fenster)                                               
    n_Eingabe.grid(row=1, column=1)

    A_Label = Label(Fenster,text="A in cm^2")
    A_Label.grid(row=2, column=0)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.grid(row=2, column=1)
    
    c_Label = Label(Fenster,text="c in mol/cm^3")
    c_Label.grid(row=3, column=0)
    c_Eingabe = Entry(Fenster)
    c_Eingabe.grid(row=3, column=1)
    
    D_Label = Label(Fenster,text="D in cm^2/s")
    D_Label.grid(row=4, column=0)
    D_Eingabe = Entry(Fenster)
    D_Eingabe.grid(row=4, column=1)
    
    viscosity_Label = Label(Fenster,text="viscosity in cm^2/s")
    viscosity_Label.grid(row=5, column=0)
    viscosity_Eingabe = Entry(Fenster)
    viscosity_Eingabe.grid(row=5, column=1)
    
    T_Label = Label(Fenster,text="T in °C")
    T_Label.grid(row=6, column=0)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.grid(row=6, column=1)
    

    
    def AcceptParams():
        global n
        global A
        global T 
        global Visc
        global D
        global c
        global get_n
        global get_c
        global get_D
        
    
        get_n = var1.get()
        get_c = var2.get()
        get_D = var3.get()
        
        if get_n ==1:
            if get_c ==0:
                if get_D ==0:
              
                    #n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    Visc      = (float(viscosity_Eingabe.get()))
                    D         = (float(D_Eingabe.get()))
                    c         = (float(c_Eingabe.get()))
        
        if get_n ==0:
            if get_c ==1:
                if get_D ==0:
              
                    n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    Visc      = (float(viscosity_Eingabe.get()))
                    D         = (float(D_Eingabe.get()))
                    #c         = (float(c_Eingabe.get()))
        
        if get_n ==0:
            if get_c ==0:
                if get_D ==1:
              
                    n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    Visc      = (float(viscosity_Eingabe.get()))
                    #D         = (float(D_Eingabe.get()))
                    c         = (float(c_Eingabe.get()))
    
    
    
    def Next():
        KL_NextLevel3()
        
        def quit():
            Fenster.destroy()
        quit()   
        
    def Back():
        KouteckyLevichWindow()
        
        def quit():
            Fenster.destroy()
        quit()   
    
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=8,column=0)
    button=Button(Fenster,text="Next",command=Next).grid(row=9,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=9,column=0)
    
    


# In[8]:


def KL_NextLevel3():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Koutecky-Levich Analysis Level 3")                         
    Fenster.geometry("400x500")
    
    
    
    
    KLStart_Label = Label(Fenster,text="KL-StartPot")
    KLStart_Label.grid(row=0, column=0)
    KLStart_Eingabe = Entry(Fenster)
    KLStart_Eingabe.grid(row=0, column=1)
    
    KLEnd_Label = Label(Fenster,text="KL-EndPot")
    KLEnd_Label.grid(row=1, column=0)
    KLEnd_Eingabe = Entry(Fenster)
    KLEnd_Eingabe.grid(row=1, column=1)
    
    REFF_Label = Label(Fenster,text="Refining factor")
    REFF_Label.grid(row=2, column=0)
    REFF_Eingabe = Entry(Fenster)
    REFF_Eingabe.grid(row=2, column=1)
    
    Ru_Label = Label(Fenster,text="Ru in Ohm")
    Ru_Label.grid(row=3, column=0)
    Ru_Eingabe = Entry(Fenster)
    Ru_Eingabe.grid(row=3, column=1)
    
    
    #RefToCotti
    #_________________________
    
    
    if get_n == 1:
    
    
        var1 = IntVar()
        Checkbutton(Fenster, text="Ref to Cottrell", variable=var1).grid(row=4, column=0, sticky=W)

    
        CottiSlope_Label = Label(Fenster,text="Cott Slope microamp/s^0.5")
        CottiSlope_Label.grid(row=5, column=0)
        CottiSlope_Eingabe = Entry(Fenster)
        CottiSlope_Eingabe.grid(row=5, column=1)
    
      
    #Smoothing
    #________________
    var4 = IntVar()
    Checkbutton(Fenster, text="Smooth data", variable=var4).grid(row=6, column=0, sticky=W)
    
    SmoothFact_Label = Label(Fenster,text="Smoothing factor")
    SmoothFact_Label.grid(row=7, column=0)
    SmoothFact_Eingabe = Entry(Fenster)
    SmoothFact_Eingabe.grid(row=7, column=1)
    
    SmoothFrom_Label = Label(Fenster,text="Smooth from")
    SmoothFrom_Label.grid(row=8, column=0)
    SmoothFrom_Eingabe = Entry(Fenster)
    SmoothFrom_Eingabe.grid(row=8, column=1)
    
    SmoothTo_Label = Label(Fenster,text="Smooth to")
    SmoothTo_Label.grid(row=9, column=0)
    SmoothTo_Eingabe = Entry(Fenster)
    SmoothTo_Eingabe.grid(row=9, column=1)
    
    ManuCorrCurr_Label = Label(Fenster,text="Cap. Curr. microamp")
    ManuCorrCurr_Label.grid(row=12, column=0)
    ManuCorrCurr_Eingabe = Entry(Fenster)
    ManuCorrCurr_Eingabe.grid(row=12, column=1)
    
    #Corrections and sep Plots
    #__________________________
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Auto Correction", variable=var5).grid(row=10, column=0, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Separate KL-Plot", variable=var6).grid(row=14, column=0, sticky=W)
    var7 = IntVar()
    Checkbutton(Fenster, text="Manual Correction", variable=var7).grid(row=11, column=0, sticky=W)
    var8 = IntVar()
    Checkbutton(Fenster, text="Baseline Correction", variable=var8).grid(row=13, column=0, sticky=W)
    var9 = IntVar()
    Checkbutton(Fenster, text="Separate LSV-Plot", variable=var9).grid(row=14, column=1, sticky=W)
    var10 = IntVar()
    Checkbutton(Fenster, text="Show uncorr LSV", variable=var10).grid(row=15, column=0, sticky=W)
    var11 = IntVar()
    Checkbutton(Fenster, text="Save as txt", variable=var11).grid(row=15, column=1, sticky=W)
    
    
    
    

    def AcceptParams():
        
        global KLStart 
        global KLEnd        
        global REFF
        global Ru 
        global RefToCotti    
        global CottiSlope    
        global Cotti_D     
        global Cotti_c 
        global Smooth_Data  
        global SmoothFact     
        global SmoothFrom    
        global SmoothTo     
        global AutoCorr 
        global OnlyKLPlot
        global Instead_n_D
        global Instead_n_c
        global ManuCorr
        global ManuCorrCurr
        global OnlyLSVPlot
        global BaseCorr
        global ShowUncorrLSV
        global AsTxtSaver
        
        
        
        KLStart    = (float(KLStart_Eingabe.get()))
        KLEnd      = (float(KLEnd_Eingabe.get()))
        REFF       = (float(REFF_Eingabe.get()))
        Ru         = (float(Ru_Eingabe.get()))

        
        
        Smooth_Data   = var4.get()
        AutoCorr      = var5.get()
        OnlyKLPlot    = var6.get()
        ManuCorr      = var7.get()
        BaseCorr      = var8.get()
        OnlyLSVPlot   = var9.get()
        ShowUncorrLSV = var10.get()
        AsTxtSaver    = var11.get()
        
        
        if ManuCorr == 1:
            ManuCorrCurr  = (float(ManuCorrCurr_Eingabe.get()))
        
        
        
        if  Smooth_Data == 1:
            SmoothFact  = (float(SmoothFact_Eingabe.get()))
            SmoothFrom  = (float(SmoothFrom_Eingabe.get()))
            SmoothTo    = (float(SmoothTo_Eingabe.get()))
        
        
        if get_n == 1:
            
            RefToCotti = var1.get()
            
            if RefToCotti  == 1:
                CottiSlope  = (float(CottiSlope_Eingabe.get()))

                
                
    def Next():
        
        if BaseCorr == 0:
            KL_NextLevel4()
        if BaseCorr == 1:
            BaseGetter()
            KL_NextLevel4()
        
  
        
    def Back():
        KL_NextLevel2()
        
        def quit():
            Fenster.destroy()
        quit()   
    
    
    
    
    
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=17,column=1)
    button=Button(Fenster,text="Next",command=Next).grid(row=18 ,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=19,column=1)
    


# In[9]:


def KL_NextLevel4():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Koutecky-Levich Analysis Level 4")                         
    Fenster.geometry("1200x600")
    
    #DIESE ZWEI FUNKTIONEN MÜSSEN HIER DRIN DEFINIERT WERDEN!!!
        
    def GeradenFit(xWert,ySchnitt,Steig):
        return  ySchnitt + Steig*xWert

    def EchtesPotFinden(WoSieSuchenSoll,WasSieSuchenSoll):                         
        idxStart = (np.abs(WoSieSuchenSoll-WasSieSuchenSoll)).argmin()               
        return WoSieSuchenSoll[idxStart]           
    
    
    f = Figure(figsize=(12, 6), dpi=100)
    #f.subplots_adjust(hspace=0.4)
    
    if OnlyLSVPlot == 1:
        LSVFenster = Toplevel()                                                         
        LSVFenster.title("LSV-Plot")                         
        LSVFenster.geometry("700x700")
        
        LSVPlot = Figure(figsize=(6, 6), dpi=100)
        
    
    global Length
    global ShortLength
    
    global ILimMittel 
    global InvILimMittel
    global InvILimFITTED
    
    #um später den Offset rausholen zu können
    global KL_Offset
    
    Length = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:])))*REFF)
    ShortLength = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:]))))
    
    global KORRPOTENZIALSUPPERARRAY
    global KORRSTROMSUPPERARRAY
    global UNKORRPOTENZIALSUPPERARRAY
    global UNKORRSTROMSUPPERARRAY
    
    KORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,Length])
    KORRSTROMSUPPERARRAY     = np.empty([NumMeas,Length])
    
    UNKORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,ShortLength])
    UNKORRSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    
    
    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN ROT RATES GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    
    if VarParRot == 1:
    
        OCPs = OCPArray
        
        #________________________________________________________________      
        
              
        ILimMittel    = np.empty(int(NumMeas))
        

        for i in range(int(NumMeas)):
            StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
            UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
            if DesiReactxx == 0:
                Stromarrays                = StromarraysROH[::] 
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::] 
                Potenzialarrays            = PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000
            
            
            if DesiReactxx == 1:
                Stromarrays                = StromarraysROH
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
                Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
            if BaseCorr == 1:
                BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
                if DesiReactxx == 0:
                    BaseStromArrays        = BaseStromArraysROH[::]
                if DesiReactxx == 1:
                    BaseStromArrays        = BaseStromArraysROH[::]
                Stromarrays            = Stromarrays - BaseStromArrays
        
            LenPotArrays            = len(Potenzialarrays)
        
            AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=LenPotArrays*REFF)#, endpoint=True)
            Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
            AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
            if Smooth_Data == 1:
            
                from scipy.interpolate import UnivariateSpline
            
                SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
                SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
                IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
                IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
                spl = UnivariateSpline(Potenzialarrays, Stromarrays)
                spl.set_smoothing_factor(SmoothFact)
                AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , REFF*len(Potenzialarrays))
                AuflgelStrarrays = spl(AufgelPotArrays)
            
            
            #___________________________________________________________
            #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
            #___________________________________________________________
            KorrAufgPotArr         = AufgelPotArrays #Widerstandskorrektur erfolgt schon weiter oben!!!
        
            KorrAufgStrArrays      = AuflgelStrarrays  
            if AutoCorr == 1:
                
                OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
                OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
                KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
           
            if ManuCorr == 1:
                KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr 
        
   
            KLStartPotenzial       = EchtesPotFinden(KorrAufgPotArr,KLStart)
            KLEndPotenzial         = EchtesPotFinden(KorrAufgPotArr,KLEnd)
            
            KLStartIndex           = np.asscalar(np.where(KorrAufgPotArr == KLStartPotenzial) [0])
            KLEndIndex             = np.asscalar(np.where(KorrAufgPotArr == KLEndPotenzial) [0])
    

            
            if KLStartIndex < KLEndIndex :
                Strommittelwerte       = np.mean(KorrAufgStrArrays[KLStartIndex:KLEndIndex])
            if KLStartIndex > KLEndIndex :
                Strommittelwerte       = np.mean(KorrAufgStrArrays[KLEndIndex:KLStartIndex])
    
            ILimMittel[i] = Strommittelwerte
            
            
            for j in range(Length):
                
                KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
                KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
            for j in range(ShortLength):
                
                UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
                UNKORRSTROMSUPPERARRAY[i,j]     = StromarraysROH[j]
            
            
        
            
            
        
            bild1 = f.add_subplot(121)
            if ShowUncorrLSV == 1:
                bild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            bild1.axvline(KLStart, color='r', linestyle='--')
            bild1.axvline(KLEnd, color='r', linestyle='--')
            
            
            if OnlyLSVPlot == 1:
                LSVbild1 = LSVPlot.add_subplot(111)
                if ShowUncorrLSV == 1:
                    LSVbild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [$\mu$A]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
                LSVbild1.axvline(KLStart, color='r', linestyle='--')
                LSVbild1.axvline(KLEnd, color='r', linestyle='--')


        global InvWurzelRot
        
        
        InvWurzelRot    = 1/(RotRatesArray**0.5)
        InvILimMittel   = 1/ILimMittel
        
     
        KLFitting, pcov = curve_fit(GeradenFit, InvWurzelRot,InvILimMittel)   
        
        KL_Offset = KLFitting[0]
        
        
        InvILimFITTED = np.empty(int(NumMeas))
        for i in range(int(NumMeas)):
            InvILimFITTED[i] = GeradenFit(InvWurzelRot[i], *KLFitting)
        
        
        bild2 = f.add_subplot(122)
        
        def KLROTPLOTTER():
        
            bild2.plot(InvWurzelRot,GeradenFit(InvWurzelRot, *KLFitting),color='r',linestyle='-',marker='')
            bild2.plot(InvWurzelRot,InvILimMittel,linestyle='',marker='.')
            bild2.set_xlabel(r'$\sqrt{\omega^{-1}}$' '[min$^{-0.5}$]', fontsize=12)
            bild2.set_ylabel('I$^{-1}$ [$\mu$A$^{-1}$]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
            #________________________________________________________
            #________________________________________________________
            #ANNOTOATIONS FOR PARAMETERS AND CALCULATION OF n D or c
            #________________________________________________________
            #________________________________________________________
        
    
        
            #IF SCHLEIFE WENN n BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
            
            global n_KL
            global D_KL
            global c_KL
            
            if get_n ==1:
      
                n_KL = np.absolute(1/(-0.201*1000000*F*A*(D**0.66666666)*Visc**(-0.16666666)*c*KLFitting[1]))

                if RefToCotti == 1:
                    
                    BeschreibefensterCotti = Toplevel()                                                         
                    BeschreibefensterCotti.title("Information")                         
                    BeschreibefensterCotti.geometry("200x200")    
                    T = Text(BeschreibefensterCotti, height=6, width=30)
                    T.pack()
                    T.insert(END, "D got referred to Cottrell")
                    
                    Di = (np.absolute(0.201*CottiSlope*KLFitting[1]*np.pi**0.5 * Visc**(-0.16666666)))**(-6)
                
                    n_KL = np.absolute(1/(-0.201*1000000*F*A*(Di**0.66666666)*Visc**(-0.16666666)*c*KLFitting[1]))
                      
    
                
                bild2.annotate('n     =%8.3f' % n_KL, xy=(0.5, 0.93), xycoords='axes fraction')           
        
        
            #IF SCHLEIFE WENN D BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_D ==1:
            
                D_KL = 1000000*(np.absolute(1/(-0.201*1000000*F*A*n*Visc**(-0.16666666)*c*KLFitting[1])))**1.5
            
                
                bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_KL, xy=(0.5, 0.93), xycoords='axes fraction') 
        
            #IF SCHLEIFE WENN c BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_c ==1:
            
                c_KL = 1000000*np.absolute(1/(-0.201*1000000*F*A*(D**0.66666666)*Visc**(-0.16666666)*n*KLFitting[1]))
            
                
                bild2.annotate('c =%5.3f*10$^{-6} mol/cm^{3}$' % c_KL, xy=(0.5, 0.93), xycoords='axes fraction') 
        
        KLROTPLOTTER()

        
        f.tight_layout()
        canvas = FigureCanvasTkAgg(f, master=Fenster)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, Fenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        if OnlyLSVPlot ==1:
            canvas = FigureCanvasTkAgg(LSVPlot, master=LSVFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, LSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        

 
        
        

        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #ONLY KL PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyKLPlot == 1:
            
            KLFenster = Toplevel()                                                         
            KLFenster.title("Koutecky-Levich Plot")                         
            KLFenster.geometry("700x700")
        
            KLPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = KLPlot.add_subplot(111)
            
            KLROTPLOTTER()
                
            canvas = FigureCanvasTkAgg(KLPlot, master=KLFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, KLFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        if AsTxtSaver == 1:
            KLROTASTXTsaver()

    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN CONCENTRATIONS GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________ 
    if VarParCon == 1:
    
        OCPs = OCPArray
        
        
        
        #________________________________________________________________      
        
               
        ILimMittel    = np.empty(int(NumMeas))
        

        for i in range(int(NumMeas)):
            StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
            UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
            if DesiReactxx == 0:
                Stromarrays                = StromarraysROH[::]
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::]
                Potenzialarrays            = PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000
            
            
            if DesiReactxx == 1:
                Stromarrays                = StromarraysROH
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
                Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
            if BaseCorr == 1:
                BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
                if DesiReactxx == 0:
                    BaseStromArrays        = BaseStromArraysROH[::]
                if DesiReactxx == 1:
                    BaseStromArrays        = BaseStromArraysROH[::]
                Stromarrays            = Stromarrays - BaseStromArrays
            
            LenPotArrays            = len(Potenzialarrays)
        
            AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=LenPotArrays*REFF)#, endpoint=True)
            Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
            AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
            if Smooth_Data == 1:
            
                from scipy.interpolate import UnivariateSpline
            
                SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
                SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
                IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
                IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
                spl = UnivariateSpline(Potenzialarrays, Stromarrays)
                spl.set_smoothing_factor(SmoothFact)
                AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , REFF*len(Potenzialarrays))
                AuflgelStrarrays = spl(AufgelPotArrays)
            
            
            #___________________________________________________________
            #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
            #___________________________________________________________
            KorrAufgPotArr         = AufgelPotArrays # Korrektur erfolgte schon weiter oben
        
            KLStartPotenzial       = EchtesPotFinden(KorrAufgPotArr,KLStart)
            KLEndPotenzial         = EchtesPotFinden(KorrAufgPotArr,KLEnd)
            
            KLStartIndex           = np.asscalar(np.where(KorrAufgPotArr == KLStartPotenzial) [0])
            KLEndIndex             = np.asscalar(np.where(KorrAufgPotArr == KLEndPotenzial) [0])
    
            KorrAufgStrArrays      = AuflgelStrarrays  
            if AutoCorr == 1:
                
                OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
                OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
                KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
                #print OCPotenziale
            
            if ManuCorr == 1:
                KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr    
                
                

            
            if KLStartIndex < KLEndIndex :
                Strommittelwerte       = np.mean(KorrAufgStrArrays[KLStartIndex:KLEndIndex])
            if KLStartIndex > KLEndIndex :
                Strommittelwerte       = np.mean(KorrAufgStrArrays[KLEndIndex:KLStartIndex])
    
            ILimMittel[i] = Strommittelwerte
        
        
        
            for j in range(Length):
                
                KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
                KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
            for j in range(ShortLength):
                
                UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
                UNKORRSTROMSUPPERARRAY[i,j]     = StromarraysROH[j]
            
        
            bild1 = f.add_subplot(121)
            if ShowUncorrLSV == 1:
                bild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            bild1.axvline(KLStart, color='r', linestyle='--')
            bild1.axvline(KLEnd, color='r', linestyle='--')
        
        
            if OnlyLSVPlot == 1:
                LSVbild1 = LSVPlot.add_subplot(111)
                if ShowUncorrLSV == 1:
                    LSVbild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [$\mu$A]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
                LSVbild1.axvline(KLStart, color='r', linestyle='--')
                LSVbild1.axvline(KLEnd, color='r', linestyle='--')
            
        
        
        global InvConcArray
            
        InvConcArray    = 1/ConcentrationsArray
        InvILimMittel   = 1/ILimMittel
        
     
        KLFitting, pcov = curve_fit(GeradenFit, InvConcArray,InvILimMittel)   

        KL_Offset = KLFitting[0]
        

        InvILimFITTED = np.empty(int(NumMeas))
        for i in range(int(NumMeas)):
            InvILimFITTED[i] = GeradenFit(InvConcArray[i], *KLFitting)
        
        bild2 = f.add_subplot(122)
        
        def KLCONCPLOTTER():
            bild2.plot(InvConcArray,GeradenFit(InvConcArray, *KLFitting),color='r',linestyle='-',marker='')
            bild2.plot(InvConcArray,InvILimMittel,linestyle='',marker='.')
            bild2.set_xlabel('$c^{-1}}$' '[cm$^{3}$ mol$^{-1}$]', fontsize=12)
            bild2.set_ylabel('I$^{-1}$ [$\mu$A$^{-1}$]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
            #________________________________________________________
            #________________________________________________________
            #ANNOTATIONS FOR PARAMETERS AND CALCULATION OF n D or c
            #________________________________________________________
            #________________________________________________________
        
        
            #IF SCHLEIFE WENN n BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
            global n_KL
            global D_KL
            global c_KL
            
            if get_n ==1:
         
                n_KL = np.absolute(1/(-0.201*1000000*F*A*(D**0.66666666)*Visc**(-0.16666666)*c *RotRate**0.5 *KLFitting[1]))

                if RefToCotti == 1:
               
                    BeschreibefensterCotti = Toplevel()                                                         
                    BeschreibefensterCotti.title("Information")                         
                    BeschreibefensterCotti.geometry("200x200")    
                    T = Text(BeschreibefensterCotti, height=6, width=30)
                    T.pack()
                    T.insert(END, "D got referred to Cottrell")
                    
                    Di = (np.absolute(0.201*CottiSlope*KLFitting[1]*np.pi**0.5 * Visc**(-0.16666666)))**(-6)
                
                    n_KL = np.absolute(1/(-0.201*1000000*F*A*(Di**0.66666666)*Visc**(-0.16666666) *c *RotRate**0.5 *KLFitting[1]))
      
                bild2.annotate('n     =%8.3f' % n_KL, xy=(0.5, 0.93), xycoords='axes fraction')           
        
        
            #IF SCHLEIFE WENN D BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_D ==1:
            
                D_KL = 1000000*(np.absolute(1/(-0.201*1000000*F*A*n*Visc**(-0.16666666)*c *RotRate**0.5 *KLFitting[1])))**1.5
            
                bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_KL, xy=(0.5, 0.93), xycoords='axes fraction') 
        
            #IF SCHLEIFE WENN c BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_c ==1:
            
                bild2.annotate('GET c NOT POSSIBLE HERE', xy=(0.5, 0.93), xycoords='axes fraction') 
        

        KLCONCPLOTTER()
        
        f.tight_layout()
        canvas = FigureCanvasTkAgg(f, master=Fenster)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, Fenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        if OnlyLSVPlot ==1:
            canvas = FigureCanvasTkAgg(LSVPlot, master=LSVFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, LSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        

        
                
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #ONLY KL PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyKLPlot == 1:
            
            KLFenster = Toplevel()                                                         
            KLFenster.title("Koutecky-Levich Plot")                         
            KLFenster.geometry("700x700")
        
            KLPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = KLPlot.add_subplot(111)
            
            KLCONCPLOTTER()
                
            canvas = FigureCanvasTkAgg(KLPlot, master=KLFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, KLFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        if AsTxtSaver == 1:
            KLCONASTXTsaver()
        
            




# In[10]:


def TafelWindow():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Tafel Curve-Parameter-Getter")                         
    Fenster.geometry("400x400")
    
    OCP_Label = Label(Fenster,text="OCPs vs E \nRef in V")
    OCP_Label.grid(row=0, column=0)
    
    RotRates_Label = Label(Fenster,text="Rotation rates \n in rpm")
    RotRates_Label.grid(row=0, column=1)
    
    Concent_Label = Label(Fenster,text="Concentrations as fract.\n of the highest one")
    Concent_Label.grid(row=0, column=2)
    
    #IF ONLY ONE OCP
    #______________________________________
    var1 = IntVar()
    Checkbutton(Fenster, text="One OCP", variable=var1).grid(row=80, column=0, sticky=W)
    OneOCP_Eingabe = Entry(Fenster)                                               
    OneOCP_Eingabe.grid(row=81, column=0)
    
    
    #IF ONLY ONE Rotation rate
    #______________________________________
    var2 = IntVar()
    Checkbutton(Fenster, text="One Rotation", variable=var2).grid(row=80, column=2, sticky=W)
    OneRot_Eingabe = Entry(Fenster)                                               
    OneRot_Eingabe.grid(row=81, column=2)
    

    OCPs     = []
    RotRates = []
    Concents = []

    for i in range(NumMeas):
        
        en = Entry(Fenster)
        en.grid(row=i+1, column=0)
        OCPs.append(en)

        en1 = Entry(Fenster)
        en1.grid(row=i+1, column=1)
        RotRates.append(en1)
        
        en2 = Entry(Fenster)
        en2.grid(row=i+1, column=2)
        Concents.append(en2)
    
    
    
    def GetOCPs():    
        OnlyOneOCP = var1.get()

        OCPots = []
        
        if OnlyOneOCP ==0:
            for entry in OCPs:
                OCPots.append(float(entry.get()))
        
        if OnlyOneOCP ==1:
            for j in range(NumMeas):
                OCPots.append(float(OneOCP_Eingabe.get()))

        global OCPArray
        OCPArray = np.array(OCPots) 

    
    
    def GetRots():    
        Rots= []
        for entry in RotRates:
            Rots.append(float(entry.get()))
        global RotRatesArray
        RotRatesArray = np.array(Rots) 
        global VarParRot
        global VarParCon
        VarParRot = 1
        VarParCon = 0
        
        
    def GetConcents():    
        OnlyOneRot = var2.get()
        
        if OnlyOneRot ==1:
            global RotRate
            RotRate = (float(OneRot_Eingabe.get()))
        
        Concs= []
        for entry in Concents:
            Concs.append(float(entry.get()))
        global ConcentrationsArray
        ConcentrationsArray = np.array(Concs) 
        global VarParRot
        global VarParCon
        VarParRot = 0
        VarParCon = 1
        
   
        
    def Next():
        TafelWindowLevel2()
        
        def quit():
            Fenster.destroy()
        quit()
        
    
    
  
    
    button=Button(Fenster,text="Accept OCPs",command=GetOCPs).grid(row=82,column=0)
    button=Button(Fenster,text="Accept RotRates",command=GetRots).grid(row=82,column=1)
    button=Button(Fenster,text="Accept Concentrations",command=GetConcents).grid(row=82,column=2)
    button=Button(Fenster,text="Next",command=Next).grid(row=83,column=1)

    
    
    



    


# In[11]:


def TafelWindowLevel2():
    Fenster = Toplevel()                                                         
    Fenster.title("Tafel-Parameter Getter Level 2")                         
    Fenster.geometry("400x680")
    
    
    n_Label = Label(Fenster, text="n")
    n_Label.grid(row=0, column=0)
    n_Eingabe = Entry(Fenster)                                               
    n_Eingabe.grid(row=0, column=1)

    c_Label = Label(Fenster, text="c")
    c_Label.grid(row=1, column=0)
    c_Eingabe = Entry(Fenster)                                               
    c_Eingabe.grid(row=1, column=1)

    A_Label = Label(Fenster,text="A in cm^2")
    A_Label.grid(row=2, column=0)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.grid(row=2, column=1)
    
    T_Label = Label(Fenster,text="T in °C")
    T_Label.grid(row=3, column=0)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.grid(row=3, column=1)
    
    Ru_Label = Label(Fenster,text="Ru in Ohm")
    Ru_Label.grid(row=4, column=0)
    Ru_Eingabe = Entry(Fenster)
    Ru_Eingabe.grid(row=4, column=1)
    

    #Smoothing
    #________________
    var1 = IntVar()
    Checkbutton(Fenster, text="Smooth data", variable=var1).grid(row=5, column=0, sticky=W)
    
    SmoothFact_Label = Label(Fenster,text="Smoothing factor")
    SmoothFact_Label.grid(row=6, column=0)
    SmoothFact_Eingabe = Entry(Fenster)
    SmoothFact_Eingabe.grid(row=6, column=1)
    
    SmoothFrom_Label = Label(Fenster,text="Smooth from")
    SmoothFrom_Label.grid(row=7, column=0)
    SmoothFrom_Eingabe = Entry(Fenster)
    SmoothFrom_Eingabe.grid(row=7, column=1)
    
    SmoothTo_Label = Label(Fenster,text="Smooth to")
    SmoothTo_Label.grid(row=8, column=0)
    SmoothTo_Eingabe = Entry(Fenster)
    SmoothTo_Eingabe.grid(row=8, column=1)
     
    #Corrections
    #__________________________
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Auto Correction", variable=var2).grid(row=9, column=0, sticky=W)
    
    var3 = IntVar()
    Checkbutton(Fenster, text="Manual Correction", variable=var3).grid(row=10, column=0, sticky=W)
    ManuCorr_Label    = Label(Fenster,text="Cap. Curr in microamp.")
    ManuCorr_Label.grid(row=11, column=0)
    ManuCorr_Eingabe  = Entry(Fenster)
    ManuCorr_Eingabe.grid(row=11, column=1)
    
    var4 = IntVar()
    Checkbutton(Fenster, text="Base Correction", variable=var4).grid(row=12, column=0, sticky=W)
        
    
    ManLimPot_Label = Label(Fenster,text="E of I Lim")
    ManLimPot_Label.grid(row=13, column=0)
    ManLimPot_Eingabe = Entry(Fenster)
    ManLimPot_Eingabe.grid(row=13, column=1)
    
    Tafelreg_Label = Label(Fenster,text="Tafel Region")
    Tafelreg_Label.grid(row=14, column=0)
    Tafelreg_Eingabe = Entry(Fenster)
    Tafelreg_Eingabe.grid(row=14, column=1)
    
    TKSTP_Label = Label(Fenster,text="Shifter")
    TKSTP_Label.grid(row=15, column=0)
    TKSTP_Eingabe = Entry(Fenster)
    TKSTP_Eingabe.grid(row=15, column=1)
    
    TafyLimLow_Label = Label(Fenster,text="Plot y-Lim low")
    TafyLimLow_Label.grid(row=16, column=0)
    TafyLimLow_Eingabe = Entry(Fenster)
    TafyLimLow_Eingabe.grid(row=16, column=1)
    
    TafyLimUp_Label = Label(Fenster,text="Plot y-Lim up")
    TafyLimUp_Label.grid(row=17, column=0)
    TafyLimUp_Eingabe = Entry(Fenster)
    TafyLimUp_Eingabe.grid(row=17, column=1)
    
    ShowTafTo_Label = Label(Fenster,text="Plot x-Lim up")
    ShowTafTo_Label.grid(row=18, column=0)
    ShowTafTo_Eingabe = Entry(Fenster)
    ShowTafTo_Eingabe.grid(row=18, column=1)
    
    REFF_Label = Label(Fenster,text="Refining")
    REFF_Label.grid(row=19, column=0)
    REFF_Eingabe = Entry(Fenster)
    REFF_Eingabe.grid(row=19, column=1)
    
    EZero_Label = Label(Fenster,text="E_Zero in V")
    EZero_Label.grid(row=24, column=0)
    EZero_Eingabe = Entry(Fenster)
    EZero_Eingabe.grid(row=24, column=1)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Show Uncorr LSV", variable=var5).grid(row=20, column=0, sticky=W)
    
    var6 = IntVar()
    Checkbutton(Fenster, text="Separate LSV Plot", variable=var6).grid(row=22, column=0, sticky=W)
    var7 = IntVar()
    Checkbutton(Fenster, text="Separate Tafel Plot", variable=var7).grid(row=22, column=1, sticky=W)
    var8 = IntVar()
    Checkbutton(Fenster, text="Show Shape Plot", variable=var8).grid(row=20, column=1, sticky=W)
    var9 = IntVar()
    Checkbutton(Fenster, text="Draw Tafel Taker", variable=var9).grid(row=21, column=0, sticky=W)
    var10 = IntVar()
    Checkbutton(Fenster, text="Draw Shape Taker", variable=var10).grid(row=21, column=1, sticky=W)
    var11 = IntVar()
    Checkbutton(Fenster, text="Save as txt", variable=var11).grid(row=25, column=0, sticky=W)
    var12 = IntVar()
    Checkbutton(Fenster, text="Calculate k_zero", variable=var12).grid(row=23, column=0, sticky=W)
    

    def AcceptParams():
    
    
        global n
        global c
        global A
        global T
        global Ru
        global Smooth_Data
        global AutoCorr
        global ManLimPot
        global Tafelreg
        global TKSTP
        global TafyLimLow
        global TafyLimUp
        global ShowTafTo
        global REFF
        global ManuCorr
        global BaseCorr
        global ShowUncorrLSV
        global SeparateLSV
        global SeparateTafel
        global WaveShapeShower
        global WaveShapeTaker
        global TafelTaker
        global AsTxtSaver
        global KZeroCalcer
        global EZero
        
        
        
    
        n           = float(n_Eingabe.get())
        c           = float(c_Eingabe.get())
        A           = float(A_Eingabe.get())
        T           = float(T_Eingabe.get()) + 273.15
        Ru          = float(Ru_Eingabe.get())
        
        Smooth_Data        = var1.get()
        if Smooth_Data     == 1:
            global SmoothFact
            global SmoothFrom
            global SmoothTo
            SmoothFact     = float(SmoothFact_Eingabe.get())
            SmoothFrom     = float(SmoothFrom_Eingabe.get())
            SmoothTo       = float(SmoothTo_Eingabe.get())
        
        
        ManuCorr           = var3.get()
        if ManuCorr        == 1:
            global ManuCorrCurr
            ManuCorrCurr   = float(ManuCorr_Eingabe.get())
        
        AutoCorr         = var2.get()
        BaseCorr         = var4.get()
        ShowUncorrLSV    = var5.get()
        SeparateLSV      = var6.get()
        SeparateTafel    = var7.get()
        WaveShapeShower  = var8.get()
        TafelTaker       = var9.get()
        WaveShapeTaker   = var10.get()
        AsTxtSaver       = var11.get()
        KZeroCalcer      = var12.get()
        ManLimPot        = float(ManLimPot_Eingabe.get())
        Tafelreg         = float(Tafelreg_Eingabe.get())
        TKSTP            = float(TKSTP_Eingabe.get())
        TafyLimLow       = float(TafyLimLow_Eingabe.get())
        TafyLimUp        = float(TafyLimUp_Eingabe.get())
        ShowTafTo        = float(ShowTafTo_Eingabe.get())
        REFF             = float(REFF_Eingabe.get())
    
    
        if KZeroCalcer == 1:
            EZero = float(EZero_Eingabe.get())
    
    
    
    
    def Next():
        
        if BaseCorr == 0:
            TafelWindowLevel3()
        if BaseCorr == 1:
            BaseGetter()
            TafelWindowLevel3()
            
    def Back():
        TafelWindow()
        
        def quit():
            Fenster.destroy()
        quit()   
        

    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=26,column=1)
    button=Button(Fenster,text="Next",command=Next).grid(row=27,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=28,column=1)
    
    
    

    
    


# In[12]:


def TafelWindowLevel3():

    
    Fenster = Toplevel()                                                         
    Fenster.title("Tafel-Plot")                         
    Fenster.geometry("1200x600")
    
    #DIESE ZWEI FUNKTIONEN MÜSSEN HIER DRIN DEFINIERT WERDEN!!!
        
    def GeradenFit(xWert,ySchnitt,Steig):
        return  ySchnitt + Steig*xWert

    def EchtesPotFinden(WoSieSuchenSoll,WasSieSuchenSoll):                         
        idxStart = (np.abs(WoSieSuchenSoll-WasSieSuchenSoll)).argmin()               
        return WoSieSuchenSoll[idxStart] 
    
    
    
    global Length
    global ShortLength
    
    Length = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:])))*REFF)
    ShortLength = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:]))))
    
    global KORRPOTENZIALSUPPERARRAY
    global KORRSTROMSUPPERARRAY
    global UNKORRPOTENZIALSUPPERARRAY
    global UNKORRSTROMSUPPERARRAY
    
    KORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,Length])
    KORRSTROMSUPPERARRAY     = np.empty([NumMeas,Length])
    
    UNKORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,ShortLength])
    UNKORRSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    
    
    f = Figure(figsize=(12, 6), dpi=100)
    bild1 = f.add_subplot(121)
    
    if SeparateLSV == 1:
        LSVFenster = Toplevel()                                                         
        LSVFenster.title("LSV-Plot")                         
        LSVFenster.geometry("700x700")
        
        LSVPlot = Figure(figsize=(6, 6), dpi=100)
        
    if SeparateTafel == 1:
        SepTafelFenster = Toplevel()                                                         
        SepTafelFenster.title("Tafel-Plot")                         
        SepTafelFenster.geometry("700x700")
        
        SepTafelPlot = Figure(figsize=(6, 6), dpi=100)
        
    if WaveShapeShower == 1:
        WaveShapeFenster = Toplevel()                                                         
        WaveShapeFenster.title("Wave-Shape-Plot")                         
        WaveShapeFenster.geometry("700x700")
        
        WaveShapePlot = Figure(figsize=(6, 6), dpi=100)
        
        
    
    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #ERSTE FOR SCHLEIFE FÜR TAFEL PLOT DATENEXTRAKTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    
    
    OCPs = OCPArray
    

    for i in range(int(NumMeas)):
        StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
        UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
        PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
        if DesiReactxx == 0:
            Stromarrays                = StromarraysROH[::]
            UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::]
            Potenzialarrays            = PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000
            
            
        if DesiReactxx == 1:
            Stromarrays                = StromarraysROH
            UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
            Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
        if BaseCorr == 1:
            BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
            if DesiReactxx == 0:
                BaseStromArrays        = BaseStromArraysROH[::]
            if DesiReactxx == 1:
                BaseStromArrays        = BaseStromArraysROH[::]
            Stromarrays            = Stromarrays - BaseStromArrays
        

        
        LenPotArrays            = len(Potenzialarrays)
        
        AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=LenPotArrays*REFF)#, endpoint=True)
        Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
        AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
        if Smooth_Data == 1:
            
            from scipy.interpolate import UnivariateSpline
            
            SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
            SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
            IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
            IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
            spl = UnivariateSpline(Potenzialarrays, Stromarrays)
            spl.set_smoothing_factor(SmoothFact)
            AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , REFF*len(Potenzialarrays))
            AuflgelStrarrays = spl(AufgelPotArrays)
            
            
        #___________________________________________________________
        #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
        #___________________________________________________________
        KorrAufgPotArr         = AufgelPotArrays #Ohmsche Korrektur erfolgte weiter oben

        KorrAufgStrArrays      = AuflgelStrarrays  
        if AutoCorr == 1:
                
            OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
            OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
            KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
        
        if ManuCorr == 1:
            KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr
        
        
        for j in range(Length):
                
            KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
            KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
        for j in range(ShortLength):
                
            UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
            UNKORRSTROMSUPPERARRAY[i,j]     = StromarraysROH[j]
        
        
        
        
        bild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
        if ShowUncorrLSV == 1:
            bild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='', color='k')
        bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
        bild1.set_ylabel('I [$\mu$A]', fontsize=12)
        for axis in ['top','bottom','left','right']:
            bild1.spines[axis].set_linewidth(2)
        bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        bild1.axvline(OCPs[i],color='k')
        bild1.axhline(0, color='k')
        
        if SeparateLSV == 1:
            LSVbild1 = LSVPlot.add_subplot(111)
            LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
            if ShowUncorrLSV == 1:
                LSVbild1.plot(UnkorrPotenzialarrays,StromarraysROH, linestyle='-',marker='', color='k')
            LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            LSVbild1.set_ylabel('I [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                LSVbild1.spines[axis].set_linewidth(2)
            LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            LSVbild1.axvline(OCPs[i],color='k')
            LSVbild1.axhline(0, color='k')

    
    
        #________________________________________
        #________________________________________
                   #JETZT TAFEL!!!!!
        #________________________________________
        #________________________________________
        
        
    global LogFitStromArray
    global LogKinStromArray
    global EtaArray
    
    
    LogFitStromArray = np.empty([len(OCPs),len(KorrAufgPotArr)])
    LogKinStromArray = np.empty([len(OCPs),len(KorrAufgPotArr)])
    EtaArray         = np.empty([len(OCPs),len(KorrAufgPotArr)])
    
    global TafelbetaArray
    global TafelIEQArray
    TafelbetaArray = np.empty(len(OCPs))
    TafelIEQArray  = np.empty(len(OCPs))

    global FormBetaArray
    global FormIEQArray
    FormBetaArray  = np.empty(len(OCPs))
    FormIEQArray   = np.empty(len(OCPs))

    global PotVonHalbLim
    global PotVonDreivLim
    PotVonHalbLim  = np.empty(len(OCPs))
    PotVonDreivLim = np.empty(len(OCPs))
    
    global EtaStartEchtARRAY
    EtaStartEchtARRAY  = np.empty(len(OCPs))
    
    global EtaEndEchtARRAY
    EtaEndEchtARRAY  = np.empty(len(OCPs))
    
    global TafelStartImLSVARRAY
    TafelStartImLSVARRAY  = np.empty(len(OCPs))
    
    global TafelEndImLSVARRAY
    TafelEndImLSVARRAY  = np.empty(len(OCPs))

    if DesiReactxx   == 0:
        EtaStart       = np.log(0.01)*(R*T)/(n*F)  +  TKSTP
        
    if DesiReactxx   == 1:
        EtaStart       = -np.log(0.01)*(R*T)/(n*F)  +  TKSTP
    
      
    EtaEnd         = EtaStart   +    Tafelreg
    
    if TafelTaker == 1:
        for i in range(int(NumMeas)):
            bild1.axvline(OCPs[i]+EtaEnd, color='k', linestyle='--')
            bild1.axvline(OCPs[i]+EtaStart, color='k', linestyle='--')
            bild1.axvline(ManLimPot, color='k', linestyle='-')
    
    
    if SeparateLSV == 1:
        if TafelTaker == 1:
            for i in range(int(NumMeas)):
                LSVbild1.axvline(OCPs[i]+EtaEnd, color='k', linestyle='--')
                LSVbild1.axvline(OCPs[i]+EtaStart, color='k', linestyle='--')
                LSVbild1.axvline(ManLimPot, color='k', linestyle='-')
    

    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #Zweite FOR SCHLEIFE FÜR TAFEL PLOT Berechnungen von allem Tafeligen!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    
    
    
      
    for i in range(int(NumMeas)):
        StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
        UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
        PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
        if DesiReactxx == 0:
            Stromarrays                = StromarraysROH[::]
            UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::]
            Potenzialarrays            = PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000
            
            
        if DesiReactxx == 1:
            Stromarrays                = StromarraysROH
            UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
            Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
        if BaseCorr == 1:
            BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
            if DesiReactxx == 0:
                BaseStromArrays        = BaseStromArraysROH[::]
            if DesiReactxx == 1:
                BaseStromArrays        = BaseStromArraysROH[::]
            Stromarrays            = Stromarrays - BaseStromArrays
        
        LenPotArrays            = len(Potenzialarrays)
        
        AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=LenPotArrays*REFF)#, endpoint=True)
        Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
        AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
        if Smooth_Data == 1:
            
            from scipy.interpolate import UnivariateSpline
            
            SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
            SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
            IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
            IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
            spl = UnivariateSpline(Potenzialarrays, Stromarrays)
            spl.set_smoothing_factor(SmoothFact)
            AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , REFF*len(Potenzialarrays))
            AuflgelStrarrays = spl(AufgelPotArrays)
            
            
        #___________________________________________________________
        #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
        #___________________________________________________________
        KorrAufgPotArr         = AufgelPotArrays # Ohmsche Korrektur erfolgte schon weiter oben

        KorrAufgStrArrays      = AuflgelStrarrays  
        if AutoCorr == 1:
                
            OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
            OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
            KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
        
        if ManuCorr == 1:
            KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr
        
        
        Eta = KorrAufgPotArr - OCPs[i]
        
        HydroStrom     = KorrAufgStrArrays
        
        EtaTafCurLim    = ManLimPot-OCPs[i]                                          #Überspannung von Limiting Pot
        EtaTafCurLimEcht= EchtesPotFinden(Eta,EtaTafCurLim)  
        EtaTafCurLiEIdx = np.asscalar(np.where(Eta == EtaTafCurLimEcht) [0])

        LimStrom        = np.absolute(HydroStrom[EtaTafCurLiEIdx])
        
        #print LimStrom #BIS HIER HIN GEHT ALLES
        
        
        if DesiReactxx        == 0:
            LimStrom        = -LimStrom
        if DesiReactxx        == 1:
            LimStrom        =  LimStrom
        
        
        #BETA AUS KURVENFORM
        
        
        LimStromHalb        = LimStrom*0.5
        LimStromDrVi        = LimStrom*0.75
        EchtLimStromHalb    = EchtesPotFinden(HydroStrom,LimStromHalb)
        EchtLimStromDrVi    = EchtesPotFinden(HydroStrom,LimStromDrVi)
        IdxEchtLimStromHalb = np.asscalar(np.where(HydroStrom == EchtLimStromHalb) [0])
        IdxEchtLimStromDrVi = np.asscalar(np.where(HydroStrom == EchtLimStromDrVi) [0])
    
        EtaEHalb            = Eta[IdxEchtLimStromHalb]
        EtaDreiviertel      = Eta[IdxEchtLimStromDrVi]
        KurvenFormbeta      = np.log(3)*R*T/(n*F*(EtaEHalb-EtaDreiviertel))
        
        FormBetaArray[i]    = np.absolute(KurvenFormbeta)
        
        #print KurvenFormbeta
        
    
        
        PotVonHalbLim[i]        = EtaEHalb+OCPs[i]
        PotVonDreivLim[i]       = EtaDreiviertel+OCPs[i]
    
        EtaStartEcht            = EchtesPotFinden(Eta,EtaStart)
        EtaStartEchtARRAY[i]    = EtaStartEcht
        EtaStartIndex           = np.asscalar(np.where(Eta == EtaStartEcht) [0])
        TafelStartImLSVARRAY[i] = EtaStartEcht + OCPs[i]
    
    
        EtaEndEcht              = EchtesPotFinden(Eta,EtaEnd)
        EtaEndEchtARRAY[i]      = EtaEndEcht
        EtaEndIndex             = np.asscalar(np.where(Eta == EtaEndEcht) [0])
        TafelEndImLSVARRAY[i]   = EtaEndEcht + OCPs[i]
     
        LogKinStrom         = np.log(np.absolute(LimStrom*HydroStrom/(LimStrom-HydroStrom+10**(-20))))
           
        if DesiReactxx          == 0:
            TafelFitting, pcov  = curve_fit(GeradenFit, Eta[EtaEndIndex:EtaStartIndex], LogKinStrom[EtaEndIndex:EtaStartIndex])
            Tafelsteigung       = TafelFitting[1]
            Tafelbeta           = -Tafelsteigung*R*T/(n*F)
        
        if DesiReactxx          == 1:
            TafelFitting, pcov  = curve_fit(GeradenFit, Eta[EtaStartIndex:EtaEndIndex], LogKinStrom[EtaStartIndex:EtaEndIndex])
            Tafelsteigung       = TafelFitting[1]
            Tafelbeta           = -Tafelsteigung*R*T/(n*F)
        
        
        TafelbetaArray[i]   = np.absolute(Tafelbeta)
    
        IEQ_Tafel           = np.exp(TafelFitting[0])
        TafelIEQArray[i]    = IEQ_Tafel*0.000001/A
    
        IEQFORM             = np.absolute(LimStrom*np.exp(KurvenFormbeta*n*F*(EtaEHalb)/(R*T)))
        FormIEQArray[i]     = IEQFORM*0.000001/A
        

        LogKinStromArray[i] = LogKinStrom
        LogFitStromArray[i] = GeradenFit(Eta, *TafelFitting)
        EtaArray[i]         = Eta
        


    #_______________________________________________
    #MITTELWERTE DER BETAS BERECHNEN
    #_______________________________________________
    
    global BETA_T
    global BETA_F
    global IEQT
    global IEQF
    global KNULLTARRAY
    global KNULLFARRAY
    global KNULLT
    global KNULLF
    
    if VarParCon ==1:
        global CDecayArray
        CDecayArray = ConcentrationsArray*c
    
    BETA_T = np.absolute(np.mean(TafelbetaArray))
    BETA_F = np.absolute(np.mean(FormBetaArray))

    IEQT   = np.mean(TafelIEQArray)
    IEQTPLOTANNO = IEQT*1000000
    IEQF   = np.mean(FormIEQArray)
    IEQFPLOTANNO = IEQF*1000000
    

    if KZeroCalcer ==1:
        if DesiReactxx == 0:
            if VarParRot ==1:
                KNULLTARRAY = (TafelIEQArray/(n*F*c))*np.exp(BETA_T*n*F*(OCPArray-EZero)/(R*T))
                KNULLT      = np.mean(KNULLTARRAY)
                KNULLFARRAY = (FormIEQArray/(n*F*c))*np.exp(BETA_F*n*F*(OCPArray-EZero)/(R*T))
                KNULLF      = np.mean(KNULLFARRAY)
            if VarParCon ==1:
                KNULLTARRAY = (TafelIEQArray/(n*F*CDecayArray))*np.exp(BETA_T*n*F*(OCPArray-EZero)/(R*T))
                KNULLT      = np.mean(KNULLTARRAY)
                KNULLFARRAY = (FormIEQArray/(n*F*CDecayArray))*np.exp(BETA_F*n*F*(OCPArray-EZero)/(R*T))
                KNULLF      = np.mean(KNULLFARRAY)
            
        
        if DesiReactxx == 1:
            if VarParRot ==1:
                KNULLTARRAY = (TafelIEQArray/(n*F*c))*np.exp(-BETA_T*n*F*(OCPArray-EZero)/(R*T))
                KNULLT      = np.mean(KNULLTARRAY)
                KNULLFARRAY = (FormIEQArray/(n*F*c))*np.exp(-BETA_F*n*F*(OCPArray-EZero)/(R*T))
                KNULLF      = np.mean(KNULLFARRAY)
            if VarParCon ==1:
                KNULLTARRAY = (TafelIEQArray/(n*F*CDecayArray))*np.exp(-BETA_T*n*F*(OCPArray-EZero)/(R*T))
                KNULLT      = np.mean(KNULLTARRAY)
                KNULLFARRAY = (FormIEQArray/(n*F*CDecayArray))*np.exp(-BETA_F*n*F*(OCPArray-EZero)/(R*T))
                KNULLF      = np.mean(KNULLFARRAY)

   

    
    #________________________________________________
    #WAVESHAPESHOWER IM LSV PLOT
    #________________________________________________
    if WaveShapeShower == 1:
        if WaveShapeTaker == 1:
            for i in range(int(NumMeas)):
                bild1.axvline(PotVonHalbLim[i], color='g', linestyle='-')
                bild1.axvline(PotVonDreivLim[i], color='g', linestyle='-')
                
                
    if SeparateLSV == 1:
        if WaveShapeShower == 1:
            if WaveShapeTaker == 1:
                for i in range(int(NumMeas)):
                    LSVbild1.axvline(PotVonHalbLim[i], color='g', linestyle='-')
                    LSVbild1.axvline(PotVonDreivLim[i], color='g', linestyle='-')
        
    
    #__________________________________________________
    #TAFELPLOTS
    #__________________________________________________
    
    Tafelbild = f.add_subplot(122)
    for i in range(NumMeas):
            Tafelbild.plot(EtaArray[i],LogKinStromArray[i], color='k',linestyle='-',marker='')
            Tafelbild.plot(EtaArray[i],LogFitStromArray[i],color='r',linestyle='-',marker='')
    
    Tafelbild.set_xlabel('$\eta$'  '[V]', fontsize=12)
    Tafelbild.set_ylabel('ln(I$_k$)', fontsize=12)
    for axis in ['top','bottom','left','right']:
        Tafelbild.spines[axis].set_linewidth(2)
    Tafelbild.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    if DesiReactxx          == 0:
        Tafelbild.set_xlim([ShowTafTo,0]) 
    if DesiReactxx          == 1:
        Tafelbild.set_xlim([0,ShowTafTo]) 
    Tafelbild.set_ylim([TafyLimLow,TafyLimUp])
    
    if DesiReactxx == 0:
        Tafelbild.annotate(u'\u03B2' '$_{Av}$' '=%5.3f'    % BETA_T, xy=(0.4, 0.93), xycoords='axes fraction',fontsize=12)
        Tafelbild.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % IEQTPLOTANNO, xy=(0.4, 0.85), xycoords='axes fraction',fontsize=12)
        
    if DesiReactxx == 1:
        Tafelbild.annotate(u'\u03B2' '$_{Av}$' '=%5.3f'  % BETA_T, xy=(0.4, 0.15), xycoords='axes fraction',fontsize=12)  
        Tafelbild.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % IEQTPLOTANNO, xy=(0.4, 0.10), xycoords='axes fraction',fontsize=12)
    
    

    
    if SeparateTafel == 1:
        SepTafelbild = SepTafelPlot.add_subplot(111)
        for i in range(NumMeas):
                SepTafelbild.plot(EtaArray[i],LogKinStromArray[i], color='k',linestyle='-',marker='')
                SepTafelbild.plot(EtaArray[i],LogFitStromArray[i],color='r',linestyle='-',marker='')
    
        SepTafelbild.set_xlabel('$\eta$'  '[V]', fontsize=12)
        SepTafelbild.set_ylabel('ln(I$_k$)', fontsize=12)
        for axis in ['top','bottom','left','right']:
            SepTafelbild.spines[axis].set_linewidth(2)
        SepTafelbild.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        if DesiReactxx          == 0:
            SepTafelbild.set_xlim([ShowTafTo,0]) 
        if DesiReactxx          == 1:
            SepTafelbild.set_xlim([0,ShowTafTo]) 
        SepTafelbild.set_ylim([TafyLimLow,TafyLimUp])
        
        if DesiReactxx == 0:
            SepTafelbild.annotate(u'\u03B2' '$_{Av}$' '=%5.3f'    % BETA_T, xy=(0.4, 0.93), xycoords='axes fraction',fontsize=12) 
            SepTafelbild.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % IEQTPLOTANNO, xy=(0.4, 0.85), xycoords='axes fraction',fontsize=12)
        if DesiReactxx == 1:
            SepTafelbild.annotate(u'\u03B2' '$_{Av}$' '=%5.3f'  % BETA_T, xy=(0.4, 0.13), xycoords='axes fraction',fontsize=12)  
            SepTafelbild.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % IEQTPLOTANNO, xy=(0.4, 0.10), xycoords='axes fraction',fontsize=12)
    
    
    
    
    if WaveShapeShower == 1: 
        WaveShapebild = WaveShapePlot.add_subplot(111)
        
        if VarParRot == 1:
            WaveShapebild.plot((1/(RotRatesArray**0.5)),np.absolute(FormBetaArray),linestyle='-',marker='.',color='g')
            WaveShapebild.set_xlabel('1/' r'$\sqrt{\omega} $' ' ' r'[1/$\sqrt{min} $]' , fontsize=12)
            
        if VarParCon == 1:
            WaveShapebild.plot(ConcentrationsArray,np.absolute(FormBetaArray),linestyle='-',marker='.',color='g')
            WaveShapebild.set_xlabel('c' '/C' '$_{max}$'  , fontsize=12)
        
        
        
        for axis in ['top','bottom','left','right']:
            WaveShapebild.spines[axis].set_linewidth(2)
        WaveShapebild.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        WaveShapebild.set_ylim([0,1])
        WaveShapebild.set_ylabel(u'\u03B2' '$_{shape}$' , fontsize=12)
        WaveShapebild.annotate(u'\u03B2' '$_{Av}$' '=%5.3f'    % BETA_F, xy=(0.4, 0.93), xycoords='axes fraction',fontsize=12) 
        WaveShapebild.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % IEQFPLOTANNO, xy=(0.4, 0.85), xycoords='axes fraction',fontsize=12)
    

        
    #ALLE CANVAS EINSTELLUNGEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    f.tight_layout()
    canvas = FigureCanvasTkAgg(f, master=Fenster)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, Fenster)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
       
    
    if SeparateLSV == 1:
        canvas = FigureCanvasTkAgg(LSVPlot, master=LSVFenster)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, LSVFenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    if SeparateTafel == 1:
        canvas = FigureCanvasTkAgg(SepTafelPlot, master=SepTafelFenster)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, SepTafelFenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
    if WaveShapeShower == 1:
        canvas = FigureCanvasTkAgg(WaveShapePlot, master=WaveShapeFenster)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, WaveShapeFenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
    
    if AsTxtSaver ==1:
    
        TAFELASTXTsaver()
    
    
    


# In[13]:


def RandlesRevWindow():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Curve-Getter")                         
    Fenster.geometry("400x400")
    
    OCP_Label = Label(Fenster,text="OCPs vs E \nRef in V")
    OCP_Label.grid(row=0, column=0)
    
    Scanrates_Label = Label(Fenster,text="Scanrates \n in mV/s")
    Scanrates_Label.grid(row=0, column=1)
    
    Concent_Label = Label(Fenster,text="Concentrations as fract.\n of the highest one")
    Concent_Label.grid(row=0, column=2)
    
    #IF ONLY ONE OCP
    #______________________________________
    var1 = IntVar()
    Checkbutton(Fenster, text="One OCP", variable=var1).grid(row=80, column=0, sticky=W)
    OneOCP_Eingabe = Entry(Fenster)                                               
    OneOCP_Eingabe.grid(row=81, column=0)
    
    
    #IF ONLY ONE Rotation rate
    #______________________________________
    var2 = IntVar()
    Checkbutton(Fenster, text="One Scanrate", variable=var2).grid(row=80, column=2, sticky=W)
    OneSca_Eingabe = Entry(Fenster)                                               
    OneSca_Eingabe.grid(row=81, column=2)
    

    OCPs     = []
    Scanrates = []
    Concents = []

    for i in range(NumMeas):
        
        en = Entry(Fenster)
        en.grid(row=i+1, column=0)
        OCPs.append(en)

        en1 = Entry(Fenster)
        en1.grid(row=i+1, column=1)
        Scanrates.append(en1)
        
        en2 = Entry(Fenster)
        en2.grid(row=i+1, column=2)
        Concents.append(en2)
    
    
    
    def GetOCPs():    
        OnlyOneOCP = var1.get()

        OCPots = []
        
        if OnlyOneOCP ==0:
            for entry in OCPs:
                OCPots.append(float(entry.get()))
        
        if OnlyOneOCP ==1:
            for j in range(NumMeas):
                OCPots.append(float(OneOCP_Eingabe.get()))

        global OCPArray
        OCPArray = np.array(OCPots) 

    
    
    def GetScanrates():    
        Scanrs= []
        for entry in Scanrates:
            Scanrs.append(float(entry.get()))
        global ScanratesArray
        ScanratesArray = np.array(Scanrs) 
        global VarParSca
        global VarParCon
        VarParSca = 1
        VarParCon = 0
        
        
    def GetConcents():    
        OnlyOneSca = var2.get()
        
        if OnlyOneSca ==1:
            global Scanrate
            Scanrate = (float(OneSca_Eingabe.get()))
        
        Concs= []
        for entry in Concents:
            Concs.append(float(entry.get()))
        global ConcentrationsArray
        ConcentrationsArray = np.array(Concs) 
        global VarParSca
        global VarParCon
        VarParSca = 0
        VarParCon = 1
        
   
        
    def Next():
        RS_REV_NextLevel2()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    button=Button(Fenster,text="Accept OCPs",command=GetOCPs).grid(row=82,column=0)
    button=Button(Fenster,text="Accept Scanrates",command=GetScanrates).grid(row=82,column=1)
    button=Button(Fenster,text="Accept Concentrations",command=GetConcents).grid(row=82,column=2)
    button=Button(Fenster,text="Next",command=Next).grid(row=83,column=1)

    
    
    

    


# In[14]:


def RS_REV_NextLevel2():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Reversible Parameter-Getter")                         
    Fenster.geometry("400x400")
    
    
    Getter_Label = Label(Fenster, text="Get From RS Rev")
    Getter_Label.grid(row=0, column=0)
    var1 = IntVar()
    Checkbutton(Fenster, text="get n", variable=var1).grid(row=0, column=1, sticky=W)
    var2 = IntVar()
    Checkbutton(Fenster, text="get c", variable=var2).grid(row=0,column=2, sticky=W)
    var3 = IntVar()
    Checkbutton(Fenster, text="get D", variable=var3).grid(row=0,column=3, sticky=W)
    
        
    n_Label = Label(Fenster, text="n")
    n_Label.grid(row=1, column=0)
    n_Eingabe = Entry(Fenster)                                               
    n_Eingabe.grid(row=1, column=1)

    A_Label = Label(Fenster,text="A in cm^2")
    A_Label.grid(row=2, column=0)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.grid(row=2, column=1)
    
    c_Label = Label(Fenster,text="c in mol/cm^3")
    c_Label.grid(row=3, column=0)
    c_Eingabe = Entry(Fenster)
    c_Eingabe.grid(row=3, column=1)
    
    D_Label = Label(Fenster,text="D in cm^2/s")
    D_Label.grid(row=4, column=0)
    D_Eingabe = Entry(Fenster)
    D_Eingabe.grid(row=4, column=1)
    
    
    T_Label = Label(Fenster,text="T in °C")
    T_Label.grid(row=6, column=0)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.grid(row=6, column=1)
    

    
    def AcceptParams():
        global n
        global A
        global T 
        global D
        global c
        global get_n
        global get_c
        global get_D
        
    
        get_n = var1.get()
        get_c = var2.get()
        get_D = var3.get()
        
        if get_n ==1:
            if get_c ==0:
                if get_D ==0:
              
                    #n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    D         = (float(D_Eingabe.get()))
                    c         = (float(c_Eingabe.get()))
        
        if get_n ==0:
            if get_c ==1:
                if get_D ==0:
              
                    n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    D         = (float(D_Eingabe.get()))
                    #c         = (float(c_Eingabe.get()))
        
        if get_n ==0:
            if get_c ==0:
                if get_D ==1:
              
                    n         = (float(n_Eingabe.get()))
                    A         = (float(A_Eingabe.get()))
                    T         = (float(T_Eingabe.get())) + 273.15
                    #D         = (float(D_Eingabe.get()))
                    c         = (float(c_Eingabe.get()))
    
    
    
    def Next():
        RS_REV_NextLevel3()
        
        def quit():
            Fenster.destroy()
        quit()   
        
    def Back():
        RandlesRevWindow()
        
        def quit():
            Fenster.destroy()
        quit()   
    
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=8,column=0)
    button=Button(Fenster,text="Next",command=Next).grid(row=9,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=9,column=0)
    
    


# In[15]:


def RS_REV_NextLevel3():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Reversible Parameter-Getter-2")                         
    Fenster.geometry("400x400")
    
    
    REFF_Label = Label(Fenster,text="Refining factor")
    REFF_Label.grid(row=2, column=0)
    REFF_Eingabe = Entry(Fenster)
    REFF_Eingabe.grid(row=2, column=1)
    
    Ru_Label = Label(Fenster,text="Ru in Ohm")
    Ru_Label.grid(row=3, column=0)
    Ru_Eingabe = Entry(Fenster)
    Ru_Eingabe.grid(row=3, column=1)
    
    
    #RefToCotti
    #_________________________
    
    
    if get_n == 1:
    
    
        var1 = IntVar()
        Checkbutton(Fenster, text="Ref to Cottrell", variable=var1).grid(row=4, column=0, sticky=W)

    
        CottiSlope_Label = Label(Fenster,text="Cott Slope microamp/s^0.5")
        CottiSlope_Label.grid(row=5, column=0)
        CottiSlope_Eingabe = Entry(Fenster)
        CottiSlope_Eingabe.grid(row=5, column=1)
    
      
    #Smoothing
    #________________
    var4 = IntVar()
    Checkbutton(Fenster, text="Smooth data", variable=var4).grid(row=6, column=0, sticky=W)
    
    SmoothFact_Label = Label(Fenster,text="Smoothing factor")
    SmoothFact_Label.grid(row=7, column=0)
    SmoothFact_Eingabe = Entry(Fenster)
    SmoothFact_Eingabe.grid(row=7, column=1)
    
    SmoothFrom_Label = Label(Fenster,text="Smooth from")
    SmoothFrom_Label.grid(row=8, column=0)
    SmoothFrom_Eingabe = Entry(Fenster)
    SmoothFrom_Eingabe.grid(row=8, column=1)
    
    SmoothTo_Label = Label(Fenster,text="Smooth to")
    SmoothTo_Label.grid(row=9, column=0)
    SmoothTo_Eingabe = Entry(Fenster)
    SmoothTo_Eingabe.grid(row=9, column=1)
     
    
    #Corrections
    #__________________________
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Auto Correction", variable=var5).grid(row=10, column=0, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Separate RS-Plot", variable=var6).grid(row=13, column=1, sticky=W)
    
    var7 = IntVar()
    Checkbutton(Fenster, text="Manual Correction", variable=var7).grid(row=11, column=0, sticky=W)
    ManuCorr_Label    = Label(Fenster,text="Cap. Curr in microamp.")
    ManuCorr_Label.grid(row=12, column=0)
    ManuCorr_Eingabe  = Entry(Fenster)
    ManuCorr_Eingabe.grid(row=12, column=1)
    
    var8 = IntVar()
    Checkbutton(Fenster, text="Base Correction", variable=var8).grid(row=10, column=1, sticky=W)
    var9 = IntVar()
    Checkbutton(Fenster, text="Separate LSV-Plot", variable=var9).grid(row=13, column=0, sticky=W)
    var10 = IntVar()
    Checkbutton(Fenster, text="Show Uncorr LSV", variable=var10).grid(row=14, column=0, sticky=W)
    var11 = IntVar()
    Checkbutton(Fenster, text="Save Results as Txt", variable=var11).grid(row=14, column=1, sticky=W)
    

    def AcceptParams():
               
        global REFF
        global Ru 
        global RefToCotti    
        global CottiSlope    
        global Cotti_D     
        global Cotti_c 
        global Smooth_Data  
        global SmoothFact     
        global SmoothFrom    
        global SmoothTo     
        global AutoCorr 
        global OnlyRSPlot
        global Instead_n_D
        global Instead_n_c
        global SeparateLSV
        global BaseCorr
        global ManuCorr
        global ManuCorrCurr
        global ShowUncorrLSV
        global AsTxtSaver
        
        
        REFF            = (float(REFF_Eingabe.get()))
        Ru              = (float(Ru_Eingabe.get()))


        Smooth_Data     = var4.get()
        AutoCorr        = var5.get()
        OnlyRSPlot      = var6.get()
        ManuCorr        = var7.get()
        BaseCorr        = var8.get()
        SeparateLSV     = var9.get()
        ShowUncorrLSV   = var10.get()
        AsTxtSaver      = var11.get()
        
        
        
        if ManuCorr == 1:
            ManuCorrCurr   =(float(ManuCorr_Eingabe.get()))
        
        
        if  Smooth_Data == 1:
            SmoothFact  = (float(SmoothFact_Eingabe.get()))
            SmoothFrom  = (float(SmoothFrom_Eingabe.get()))
            SmoothTo    = (float(SmoothTo_Eingabe.get()))
        
        
        if get_n == 1:
            
            RefToCotti = var1.get()
            
            if RefToCotti  == 1:
                CottiSlope  = (float(CottiSlope_Eingabe.get()))

                
    def Next():
        
        if BaseCorr == 0:
            RS_REV_NextLevel4()
        if BaseCorr == 1:
            BaseGetter()
            RS_REV_NextLevel4()
             
        
    def Back():
        RS_REV_NextLevel2()
        
        def quit():
            Fenster.destroy()
        quit()   
    

    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=15,column=1)
    button=Button(Fenster,text="Next",command=Next).grid(row=16 ,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=17,column=1)
    


# In[16]:


def RS_REV_NextLevel4():
    Fenster = Toplevel()                                                         
    Fenster.title("Reversible Randles-Sevcik Analysis")                         
    Fenster.geometry("1200x600")
    
    #DIESE ZWEI FUNKTIONEN MÜSSEN HIER DRIN DEFINIERT WERDEN!!!
        
    def GeradenFit(xWert,ySchnitt,Steig):
        return  ySchnitt + Steig*xWert

    def EchtesPotFinden(WoSieSuchenSoll,WasSieSuchenSoll):                         
        idxStart = (np.abs(WoSieSuchenSoll-WasSieSuchenSoll)).argmin()               
        return WoSieSuchenSoll[idxStart]           
    

    
    global Length
    global ShortLength
    
    Length = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:])))*REFF)
    ShortLength = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:]))))
    
    
    global KORRPOTENZIALSUPPERARRAY
    global KORRSTROMSUPPERARRAY
    global UNKORRPOTENZIALSUPPERARRAY
    global UNKORRSTROMSUPPERARRAY
    global I_Peak
    global FITTED_I_Peak
    global WurzelScanratesarray
    
    
    KORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,Length])
    KORRSTROMSUPPERARRAY     = np.empty([NumMeas,Length])
    
    UNKORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,ShortLength])
    UNKORRSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
 
    
    f = Figure(figsize=(12, 6), dpi=100)
    
    if SeparateLSV == 1:
        LSVFenster = Toplevel()                                                         
        LSVFenster.title("LSV-Plot")                         
        LSVFenster.geometry("700x700")
        
        LSVPlot = Figure(figsize=(6, 6), dpi=100)
    
    #f.subplots_adjust(hspace=0.4)
    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN ROT RATES GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    
    if VarParSca == 1:
    
        OCPs = OCPArray
        
        #________________________________________________________________      
        
               
        I_Peak    = np.empty(int(NumMeas))
        
        
        for i in range(int(NumMeas)):
            
            StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
            UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
            if DesiReactxx == 0:
                Stromarrays                = StromarraysROH[::-1]
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::-1]
                Potenzialarrays            = PotenzialarraysROH[::-1]   -  Ru*Stromarrays/1000000
            
            
            if DesiReactxx == 1:
                Stromarrays                = StromarraysROH
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
                Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
            if BaseCorr == 1:
                BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
                if DesiReactxx == 0:
                    BaseStromArrays        = BaseStromArraysROH[::-1]
                if DesiReactxx == 1:
                    BaseStromArrays        = BaseStromArraysROH[::]
                Stromarrays            = Stromarrays - BaseStromArrays
            
  
            
            LenPotArrays            = len(Potenzialarrays)
        
            AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=LenPotArrays*REFF)#, endpoint=True)
            Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
            AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
            if Smooth_Data == 1:
            
                from scipy.interpolate import UnivariateSpline
            
                SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
                SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
                IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
                IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
                spl = UnivariateSpline(Potenzialarrays, Stromarrays)
                spl.set_smoothing_factor(SmoothFact)
                AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , REFF*len(Potenzialarrays))
                AuflgelStrarrays = spl(AufgelPotArrays)
            
            
            #___________________________________________________________
            #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
            #___________________________________________________________
            KorrAufgPotArr         = AufgelPotArrays #Ohmsche Korrektur erfolgte schon weiter oben
            
 
            KorrAufgStrArrays      = AuflgelStrarrays  
            if AutoCorr == 1:
                
                OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
                OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
                KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
             
            if ManuCorr == 1:
                KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr
              
            
            for j in range(Length):
                
                KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
                KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
            for j in range(ShortLength):
                
                UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
                UNKORRSTROMSUPPERARRAY[i,j]     = Stromarrays[j]
            
            
            if DesiReactxx == 0:
                PeakMax    = np.min(KorrAufgStrArrays) 
            
            if DesiReactxx == 1:
                PeakMax    = np.max(KorrAufgStrArrays) 
            
            

    
            I_Peak[i] = PeakMax
        
            
        
            bild1 = f.add_subplot(121)
            if ShowUncorrLSV ==1:
                bild1.plot(UnkorrPotenzialarrays,Stromarrays, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            
            if SeparateLSV == 1:
                LSVbild1 = LSVPlot.add_subplot(111)
                if ShowUncorrLSV ==1:
                    LSVbild1.plot(UnkorrPotenzialarrays,Stromarrays, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [$\mu$A]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
        
        
        

        WurzelScanratesarray = (ScanratesArray*0.001)**0.5

        RSREVFitting, pcov = curve_fit(GeradenFit, WurzelScanratesarray,I_Peak)   
        
        FITTED_I_Peak = GeradenFit(WurzelScanratesarray, *RSREVFitting)
        
        bild2 = f.add_subplot(122)
        
        def RSREVSCPLOTTER():
        
            bild2.plot(WurzelScanratesarray,GeradenFit(WurzelScanratesarray, *RSREVFitting),color='r',linestyle='-',marker='')
            bild2.plot(WurzelScanratesarray,I_Peak,linestyle='',marker='.')
            bild2.set_xlabel(r'$\sqrt{\nu} $' ' ' r'[$\sqrt{V} $]' , fontsize=12)
            bild2.set_ylabel('I$_p$ [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
            #________________________________________________________
            #________________________________________________________
            #ANNOTOATIONS FOR PARAMETERS AND CALCULATION OF n D or c
            #________________________________________________________
            #________________________________________________________
        
        
            #IF SCHLEIFE WENN n BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
            if get_n ==1:
         
                global n_RS
                n_RS = (np.absolute(RSREVFitting[1]*0.000001*(R*T)**0.5 / (0.4463*F*A*c*(F*D)**0.5)))**0.666666667

                if RefToCotti == 1:
               
                    BeschreibefensterCotti = Toplevel()                                                         
                    BeschreibefensterCotti.title("Information")                         
                    BeschreibefensterCotti.geometry("200x200")    
                    Te = Text(BeschreibefensterCotti, height=6, width=30)
                    Te.pack()
                    Te.insert(END, "Randles Slope got\ndivided by Cottrell Slope")
                    
                    
                
                    n_RS = np.absolute((RSREVFitting[1]*(R*T)**0.5)/(np.pi**0.5 * 0.4463 * CottiSlope * F**0.5))**2
      
                
                if DesiReactxx == 0:
                    bild2.annotate('n     =%8.3f' % n_RS, xy=(0.5, 0.93), xycoords='axes fraction')     
                if DesiReactxx == 1:
                    bild2.annotate('n     =%8.3f' % n_RS, xy=(0.5, 0.13), xycoords='axes fraction')  
        
        
            #IF SCHLEIFE WENN D BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_D ==1:
            
                global D_RS
                D_RS = 1000000*(0.000001*RSREVFitting[1]*(R*T)**0.5 / (0.4463*n*F*A*c*(n*F)**0.5))**2
            
                
                if DesiReactxx == 0:
                    bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_RS, xy=(0.5, 0.93), xycoords='axes fraction') 
                if DesiReactxx == 1:
                    bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_RS, xy=(0.5, 0.13), xycoords='axes fraction') 
                

            #IF SCHLEIFE WENN c BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_c ==1:
            
                global c_RS
                c_RS = np.absolute(1000000*(RSREVFitting[1]*0.000001*(R*T)**0.5) / (0.4463*n*F*A*(n*F*D)**0.5))
            
                if DesiReactxx == 0:
                    bild2.annotate('c =%5.3f*10$^{-6} mol/cm^{3}$' % c_RS, xy=(0.5, 0.93), xycoords='axes fraction') 
                if DesiReactxx == 1:
                    bild2.annotate('c =%5.3f*10$^{-6} mol/cm^{3}$' % c_RS, xy=(0.5, 0.13), xycoords='axes fraction') 
                

        
        RSREVSCPLOTTER()

        f.tight_layout()
        canvas = FigureCanvasTkAgg(f, master=Fenster)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, Fenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        


        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #ONLY RS PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyRSPlot == 1:
            
            RSFenster = Toplevel()                                                         
            RSFenster.title("Reversible Randles-Sevcik Plot")                         
            RSFenster.geometry("700x700")
        
            RSREVPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = RSREVPlot.add_subplot(111)
            
            RSREVSCPLOTTER()
                
            canvas = FigureCanvasTkAgg(RSREVPlot, master=RSFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, RSFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        if SeparateLSV == 1:
            canvas = FigureCanvasTkAgg(LSVPlot, master=LSVFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, LSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        if AsTxtSaver ==1:
            RSREVSCASTXTsaver()
    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN CONCENTRATIONS GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________ 
    if VarParCon == 1:
    
        OCPs = OCPArray
        
        #________________________________________________________________      
        
               
        I_Peak    = np.empty(int(NumMeas))
        

        for i in range(int(NumMeas)):
            
            StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
            UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
            if DesiReactxx == 0:
                Stromarrays                = StromarraysROH[::-1]
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::-1]
                Potenzialarrays            = PotenzialarraysROH[::-1]   -  Ru*Stromarrays/1000000
            
            
            if DesiReactxx == 1:
                Stromarrays                = StromarraysROH
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
                Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
            if BaseCorr == 1:
                BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
                if DesiReactxx == 0:
                    BaseStromArrays        = BaseStromArraysROH[::-1]
                if DesiReactxx == 1:
                    BaseStromArrays        = BaseStromArraysROH[::]
                Stromarrays            = Stromarrays - BaseStromArrays
            
            LenPotArrays            = len(Potenzialarrays)
        
            AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=LenPotArrays*REFF)#, endpoint=True)
            Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
            AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
            if Smooth_Data == 1:
            
                from scipy.interpolate import UnivariateSpline
            
                SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
                SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
                IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
                IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
                spl = UnivariateSpline(Potenzialarrays, Stromarrays)
                spl.set_smoothing_factor(SmoothFact)
                AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , REFF*len(Potenzialarrays))
                AuflgelStrarrays = spl(AufgelPotArrays)
            
            
            #___________________________________________________________
            #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
            #___________________________________________________________
            KorrAufgPotArr         = AufgelPotArrays #Ohmsche Korrektur erfolgte schon weiter oben
            
 
            KorrAufgStrArrays      = AuflgelStrarrays  
            if AutoCorr == 1:
                
                OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
                OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
                KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
            
            if ManuCorr == 1:
                KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr
            
            for j in range(Length):
                
                KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
                KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
            for j in range(ShortLength):
                
                UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
                UNKORRSTROMSUPPERARRAY[i,j]     = Stromarrays[j]
            
            
            if DesiReactxx == 0:
                PeakMax    = np.min(KorrAufgStrArrays) 
            
            if DesiReactxx == 1:
                PeakMax    = np.max(KorrAufgStrArrays) 
            

            I_Peak[i] = PeakMax
        

            bild1 = f.add_subplot(121)
            if ShowUncorrLSV ==1:
                bild1.plot(UnkorrPotenzialarrays,Stromarrays, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            
        
            if SeparateLSV == 1:
                LSVbild1 = LSVPlot.add_subplot(111)
                if ShowUncorrLSV ==1:
                    LSVbild1.plot(UnkorrPotenzialarrays,Stromarrays, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [$\mu$A]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
                
                
                
                
                

        Concentrationsarray      = ConcentrationsArray * c
        

        RSREVFitting, pcov = curve_fit(GeradenFit, Concentrationsarray,I_Peak)   
        
        cMax = c*1000000
        
        FITTED_I_Peak = GeradenFit(Concentrationsarray, *RSREVFitting)
        
        bild2 = f.add_subplot(122)
        
        def RSREVCONCPLOTTER():
        
            bild2.plot(ConcentrationsArray,GeradenFit(Concentrationsarray, *RSREVFitting),color='r',linestyle='-',marker='')
            bild2.plot(ConcentrationsArray,I_Peak,linestyle='',marker='.')
            bild2.set_xlabel('fract.(c$_{max}$)', fontsize=12)
            bild2.set_ylabel('I$_p$ [$\mu$A]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
            if DesiReactxx == 0:
                bild2.annotate( 'c$_{max}$ =%5.3f*10$^{-6} mol/cm^{3}$' % cMax , xy=(0.5, 0.83), xycoords='axes fraction')
            if DesiReactxx == 1:
                bild2.annotate( 'c$_{max}$ =%5.3f*10$^{-6} mol/cm^{3}$' % cMax , xy=(0.5, 0.23), xycoords='axes fraction')          
            
        
        
            #________________________________________________________
            #________________________________________________________
            #ANNOTOATIONS FOR PARAMETERS AND CALCULATION OF n D or c
            #________________________________________________________
            #________________________________________________________
        
            
            #IF SCHLEIFE WENN n BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
            if get_n ==1:
         
                global n_RS
                n_RS = (np.absolute(RSREVFitting[1]*0.000001*(R*T)**0.5 / (0.4463*F*A*((Scanrate*0.001)**0.5)*(F*D)**0.5)))**0.666666667

                if RefToCotti == 1:
               
                    BeschreibefensterCotti = Toplevel()                                                         
                    BeschreibefensterCotti.title("Information")                         
                    BeschreibefensterCotti.geometry("200x200")    
                    Te = Text(BeschreibefensterCotti, height=6, width=30)
                    Te.pack()
                    Te.insert(END, "Randles Slope got\ndivided by Cottrell Slope")
                    
                    
                
                    n_RS = (np.absolute(RSREVFitting[1]*c*(R*T)**0.5/(0.4463*np.pi**0.5 *(Scanrate*0.001)**0.5 *CottiSlope *F**0.5)))**2
      
                if DesiReactxx == 0:
                    bild2.annotate('n     =%8.3f' % n_RS, xy=(0.5, 0.93), xycoords='axes fraction')     
                if DesiReactxx == 1:
                    bild2.annotate('n     =%8.3f' % n_RS, xy=(0.5, 0.13), xycoords='axes fraction')           
        
        
            #IF SCHLEIFE WENN D BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_D ==1:
            
                global D_RS
                D_RS = 1000000*(0.000001*RSREVFitting[1]*(R*T)**0.5 / (0.4463*n*F*A*(n*F*Scanrate*0.001)**0.5))**2
            
                if DesiReactxx == 0:
                    bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_RS, xy=(0.5, 0.93), xycoords='axes fraction') 
                if DesiReactxx == 1:
                    bild2.annotate('D =%5.3f*10$^{-6} cm^{2}/s$' % D_RS, xy=(0.5, 0.13), xycoords='axes fraction') 
        
            #IF SCHLEIFE WENN c BESTIMMT WERDEN SOLL
            #__________________________________________________________________________________________________________________
        
            if get_c ==1:

                if DesiReactxx == 0:
                    bild2.annotate('GETTING c NOT POSSIBLE', xy=(0.5, 0.93), xycoords='axes fraction') 
                if DesiReactxx == 1:
                    bild2.annotate('GETTING c NOT POSSIBLE', xy=(0.5, 0.13), xycoords='axes fraction')  
        
        
        
        
        
        RSREVCONCPLOTTER()

        f.tight_layout()
        canvas = FigureCanvasTkAgg(f, master=Fenster)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, Fenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        
        
        
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #ONLY RSREV PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyRSPlot == 1:
            
            RSFenster = Toplevel()                                                         
            RSFenster.title("Reversible Randles-Sevcik Plot")                         
            RSFenster.geometry("700x700")
        
            RSREVPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = RSREVPlot.add_subplot(111)
            
            RSREVCONCPLOTTER()
                
            canvas = FigureCanvasTkAgg(RSREVPlot, master=RSFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, RSFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        if SeparateLSV == 1:
            canvas = FigureCanvasTkAgg(LSVPlot, master=LSVFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, LSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        
        if AsTxtSaver ==1:
            RSREVCONCASTXTsaver()
      


# In[17]:


def RS_IRR_Window():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Curve-Getter")                         
    Fenster.geometry("400x400")
    
    OCP_Label = Label(Fenster,text="OCPs vs E \nRef in V")
    OCP_Label.grid(row=0, column=0)
    
    Scanrates_Label = Label(Fenster,text="Scanrates \n in mV/s")
    Scanrates_Label.grid(row=0, column=1)
    
    Concent_Label = Label(Fenster,text="Concentrations as fract.\n of the highest one")
    Concent_Label.grid(row=0, column=2)
    
    #IF ONLY ONE OCP
    #______________________________________
    var1 = IntVar()
    Checkbutton(Fenster, text="One OCP", variable=var1).grid(row=80, column=0, sticky=W)
    OneOCP_Eingabe = Entry(Fenster)                                               
    OneOCP_Eingabe.grid(row=81, column=0)
    
    
    #IF ONLY ONE Rotation rate
    #______________________________________
    var2 = IntVar()
    Checkbutton(Fenster, text="One Scanrate", variable=var2).grid(row=80, column=2, sticky=W)
    OneSca_Eingabe = Entry(Fenster)                                               
    OneSca_Eingabe.grid(row=81, column=2)
    

    OCPs     = []
    Scanrates = []
    Concents = []

    for i in range(NumMeas):
        
        en = Entry(Fenster)
        en.grid(row=i+1, column=0)
        OCPs.append(en)

        en1 = Entry(Fenster)
        en1.grid(row=i+1, column=1)
        Scanrates.append(en1)
        
        en2 = Entry(Fenster)
        en2.grid(row=i+1, column=2)
        Concents.append(en2)
    
    
    
    def GetOCPs():    
        OnlyOneOCP = var1.get()

        OCPots = []
        
        if OnlyOneOCP ==0:
            for entry in OCPs:
                OCPots.append(float(entry.get()))
        
        if OnlyOneOCP ==1:
            for j in range(NumMeas):
                OCPots.append(float(OneOCP_Eingabe.get()))

        global OCPArray
        OCPArray = np.array(OCPots) 

    
    
    def GetScanrates():    
        Scanrs= []
        for entry in Scanrates:
            Scanrs.append(float(entry.get()))
        global ScanratesArray
        ScanratesArray = np.array(Scanrs) 
        global VarParSca
        global VarParCon
        VarParSca = 1
        VarParCon = 0
        
        
    def GetConcents():    
        OnlyOneSca = var2.get()
        
        if OnlyOneSca ==1:
            global Scanrate
            Scanrate = (float(OneSca_Eingabe.get()))
        
        Concs= []
        for entry in Concents:
            Concs.append(float(entry.get()))
        global ConcentrationsArray
        ConcentrationsArray = np.array(Concs) 
        global VarParSca
        global VarParCon
        VarParSca = 0
        VarParCon = 1
        
   
        
    def Next():
        RS_IRR_WindowLevel2()
        
        def quit():
            Fenster.destroy()
        quit()
    
    
    
    
    button=Button(Fenster,text="Accept OCPs",command=GetOCPs).grid(row=82,column=0)
    button=Button(Fenster,text="Accept Scanrates",command=GetScanrates).grid(row=82,column=1)
    button=Button(Fenster,text="Accept Concentrations",command=GetConcents).grid(row=82,column=2)
    button=Button(Fenster,text="Next",command=Next).grid(row=83,column=1)

    
    
    

    


# In[18]:


def RS_IRR_WindowLevel2():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Irreversible Parameter-Getter")                         
    Fenster.geometry("400x400")

    
    n_Label = Label(Fenster, text="n")
    n_Label.grid(row=1, column=0)
    n_Eingabe = Entry(Fenster)                                               
    n_Eingabe.grid(row=1, column=1)

    A_Label = Label(Fenster,text="A in cm^2")
    A_Label.grid(row=2, column=0)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.grid(row=2, column=1)
    
    c_Label = Label(Fenster,text="c in mol/cm^3")
    c_Label.grid(row=3, column=0)
    c_Eingabe = Entry(Fenster)
    c_Eingabe.grid(row=3, column=1)
    
    D_Label = Label(Fenster,text="D in cm^2/s")
    D_Label.grid(row=4, column=0)
    D_Eingabe = Entry(Fenster)
    D_Eingabe.grid(row=4, column=1)
    
    
    T_Label = Label(Fenster,text="T in °C")
    T_Label.grid(row=6, column=0)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.grid(row=6, column=1)
    

    
    def AcceptParams():
        global n
        global A
        global T 
        global D
        global c


        n         = (float(n_Eingabe.get()))
        A         = (float(A_Eingabe.get()))
        T         = (float(T_Eingabe.get())) + 273.15
        D         = (float(D_Eingabe.get()))
        c         = (float(c_Eingabe.get()))
        

    def Next():
        RS_IRR_WindowLevel3()
        
        def quit():
            Fenster.destroy()
        quit()   
        
    def Back():
        RS_IRR_Window()
        
        def quit():
            Fenster.destroy()
        quit()   
    
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=8,column=0)
    button=Button(Fenster,text="Next",command=Next).grid(row=9,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=9,column=0)
    
    


# In[19]:


def RS_IRR_WindowLevel3():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Irreversible Parameter-Getter-2")                         
    Fenster.geometry("400x520")
    
    
    REFF_Label = Label(Fenster,text="Refining factor")
    REFF_Label.grid(row=2, column=0)
    REFF_Eingabe = Entry(Fenster)
    REFF_Eingabe.grid(row=2, column=1)
    
    Ru_Label = Label(Fenster,text="Ru in Ohm")
    Ru_Label.grid(row=3, column=0)
    Ru_Eingabe = Entry(Fenster)
    Ru_Eingabe.grid(row=3, column=1)
    
    
    #RefToCotti
    #_________________________

    if VarParSca  == 1:
        var1 = IntVar()
        Checkbutton(Fenster, text="Ref to Cottrell", variable=var1).grid(row=4, column=0, sticky=W)

    
        CottiSlope_Label = Label(Fenster,text="Cott Slope microamp/s^0.5")
        CottiSlope_Label.grid(row=5, column=0)
        CottiSlope_Eingabe = Entry(Fenster)
        CottiSlope_Eingabe.grid(row=5, column=1)
    
      
    #Smoothing
    #________________
    var4 = IntVar()
    Checkbutton(Fenster, text="Smooth data", variable=var4).grid(row=6, column=0, sticky=W)
    
    SmoothFact_Label = Label(Fenster,text="Smoothing factor")
    SmoothFact_Label.grid(row=7, column=0)
    SmoothFact_Eingabe = Entry(Fenster)
    SmoothFact_Eingabe.grid(row=7, column=1)
    
    SmoothFrom_Label = Label(Fenster,text="Smooth from")
    SmoothFrom_Label.grid(row=8, column=0)
    SmoothFrom_Eingabe = Entry(Fenster)
    SmoothFrom_Eingabe.grid(row=8, column=1)
    
    SmoothTo_Label = Label(Fenster,text="Smooth to")
    SmoothTo_Label.grid(row=9, column=0)
    SmoothTo_Eingabe = Entry(Fenster)
    SmoothTo_Eingabe.grid(row=9, column=1)
    

     
    
    #Corrections
    #__________________________
    
    var5 = IntVar()
    Checkbutton(Fenster, text="Auto Correction", variable=var5).grid(row=10, column=0, sticky=W)
    var6 = IntVar()
    Checkbutton(Fenster, text="Separate RS-Plot", variable=var6).grid(row=13, column=1, sticky=W)
    
    var7 = IntVar()
    Checkbutton(Fenster, text="Manual Correction", variable=var7).grid(row=11, column=0, sticky=W)
    ManuCorr_Label    = Label(Fenster,text="DS-cap.(+ox or -red)\nin (microamp.*s)/mV")
    ManuCorr_Label.grid(row=12, column=0)
    ManuCorr_Eingabe  = Entry(Fenster)
    ManuCorr_Eingabe.grid(row=12, column=1)
    
    var8 = IntVar()
    Checkbutton(Fenster, text="Base Correction", variable=var8).grid(row=10, column=1, sticky=W)
    var9 = IntVar()
    Checkbutton(Fenster, text="Separate LSV-Plot", variable=var9).grid(row=13, column=0, sticky=W)
    var10 = IntVar()
    Checkbutton(Fenster, text="Show Uncorr LSV", variable=var10).grid(row=15, column=0, sticky=W)
    var11 = IntVar()
    Checkbutton(Fenster, text="Matsuda-Ayabe-1", variable=var11).grid(row=14, column=0, sticky=W)
    var12 = IntVar()
    Checkbutton(Fenster, text="Matsuda-Ayabe-2", variable=var12).grid(row=14, column=1, sticky=W)
    var13 = IntVar()
    Checkbutton(Fenster, text="Save Results as Txt", variable=var13).grid(row=18, column=0, sticky=W)
    var14 = IntVar()
    Checkbutton(Fenster, text="Calculate k_zero", variable=var14).grid(row=16, column=0, sticky=W)
    
    EZero_Label = Label(Fenster,text="E_0 vs ref. in V (for calc k_0)")
    EZero_Label.grid(row=17, column=0)
    EZero_Eingabe = Entry(Fenster)
    EZero_Eingabe.grid(row=17, column=1)
    
    
    
    

    def AcceptParams():
               
        global REFF
        global Ru 
        if VarParSca  == 1:
            global RefToCotti    
            global CottiSlope    
            global Cotti_D     
            global Cotti_c 
        global Smooth_Data  
        global SmoothFact     
        global SmoothFrom    
        global SmoothTo     
        global AutoCorr 
        global OnlyRSPlot
        global SeparateLSV
        global SeparateMA1
        global SeparateMA2
        global BaseCorr
        global ManuCorr
        global ManuCorrCurr
        global ShowUncorrLSV
        global AsTxtSaver
        global KZeroCalcer
        global EZero
        
        
        
        REFF        = (float(REFF_Eingabe.get()))
        Ru          = (float(Ru_Eingabe.get()))


        Smooth_Data   = var4.get()
        AutoCorr      = var5.get()
        OnlyRSPlot    = var6.get()
        ManuCorr      = var7.get()
        BaseCorr      = var8.get()
        SeparateLSV   = var9.get()
        ShowUncorrLSV = var10.get()
        SeparateMA1   = var11.get()
        SeparateMA2   = var12.get()
        AsTxtSaver    = var13.get()
        KZeroCalcer   = var14.get()
        
        if KZeroCalcer ==1:
            EZero = (float(EZero_Eingabe.get()))

        
        if ManuCorr == 1:
            ManuCorrCurr = (float(ManuCorr_Eingabe.get()))
        
        
        if  Smooth_Data == 1:
            SmoothFact  = (float(SmoothFact_Eingabe.get()))
            SmoothFrom  = (float(SmoothFrom_Eingabe.get()))
            SmoothTo    = (float(SmoothTo_Eingabe.get()))
        
        
        if VarParSca  == 1:    
            RefToCotti = var1.get()
            
            if RefToCotti  == 1:
                CottiSlope  = (float(CottiSlope_Eingabe.get()))

                
    def Next():
        
        if BaseCorr == 0:
            RS_IRR_WindowLevel4()
        if BaseCorr == 1:
            BaseGetter()
            RS_IRR_WindowLevel4()
        

        
        
    def Back():
        RS_IRR_WindowLevel2()
        
        def quit():
            Fenster.destroy()
        quit()   
    

    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=19,column=1)
    button=Button(Fenster,text="Next",command=Next).grid(row=20 ,column=1)
    button=Button(Fenster,text="Back",command=Back).grid(row=21,column=1)


# In[20]:


def RS_IRR_WindowLevel4():
    Fenster = Toplevel()                                                         
    Fenster.title("Randles-Sevcik Irreversible Level4")                         
    
    Fenster.geometry("700x700")
    
    
    
    #DIESE FUNKTIONEN MÜSSEN HIER DRIN DEFINIERT WERDEN!!!
        
    def GeradenFit(xWert,ySchnitt,Steig):
        return  ySchnitt + Steig*xWert

    def EchtesPotFinden(WoSieSuchenSoll,WasSieSuchenSoll):                         
        idxStart = (np.abs(WoSieSuchenSoll-WasSieSuchenSoll)).argmin()               
        return WoSieSuchenSoll[idxStart]     
    
    def find_nearest(AuflStromBisPeak,HalbMaxAuflStrom):                         
        idx = (np.abs(AuflStromBisPeak-HalbMaxAuflStrom)).argmin()               
        return AuflStromBisPeak[idx]

    
    global Length
    global ShortLength
    
    Length = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:])))*REFF)
    ShortLength = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:]))))
    
    
    global KORRPOTENZIALSUPPERARRAY
    global KORRSTROMSUPPERARRAY
    global UNKORRPOTENZIALSUPPERARRAY
    global UNKORRSTROMSUPPERARRAY
    global I_Peak
    global FITTED_I_Peak
    global WurzelScanratesarray
    global LnScanratesarray 
    global PeakPotentiale
    global FIT_PeakPotentiale
    global HalbPeakPotenziale
    
    
    
    
    KORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,Length])
    KORRSTROMSUPPERARRAY     = np.empty([NumMeas,Length])
    
    UNKORRPOTENZIALSUPPERARRAY = np.empty([NumMeas,ShortLength])
    UNKORRSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    
    
        
    if VarParSca == 1:
        f = Figure(figsize=(7, 7), dpi=80)
    
    if VarParCon == 1:
        f = Figure(figsize=(7, 7), dpi=80)
        
    
    if SeparateLSV == 1:
        SepLSVFenster = Toplevel()                                                         
        SepLSVFenster.title("LSV-Plot")                         
        SepLSVFenster.geometry("700x700")
        
        SepLSVPlot = Figure(figsize=(6, 6), dpi=100)
    
    
    if SeparateMA1 == 1:
        MA1Fenster = Toplevel()                                                         
        MA1Fenster.title("Matsuda-Ayabe-Plot 1")                         
        MA1Fenster.geometry("700x700")
        
        MA1Plot = Figure(figsize=(6, 6), dpi=100)
    
    if SeparateMA2 == 1:
        MA2Fenster = Toplevel()                                                         
        MA2Fenster.title("Matsuda-Ayabe-Plot 2")                         
        MA2Fenster.geometry("700x700")
        
        MA2Plot = Figure(figsize=(6, 6), dpi=100)
        
    

    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN SCANRATES GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    
    if VarParSca == 1:
    
        global SHAPE_BETAARRAY
        global RS_IEQARRAY
        global RS_KNULLARRAY
        global LNPLOT_IEQARRAY
        global LNPLOT_KNULLARRAY
        global SHAPE_IEQARRAY
        global SHAPE_KNULLARRAY

   
    
        PeakPotentiale     = np.empty(int(NumMeas))
        HalbPeakPotenziale = np.empty(int(NumMeas))
        SHAPE_BETAARRAY    = np.empty(int(NumMeas))
        LNPLOT_IEQARRAY    = np.empty(int(NumMeas))
        RS_IEQARRAY        = np.empty(int(NumMeas))
        
        
    
        
        OCPs = OCPArray
        
        #________________________________________________________________      
        
               
        I_Peak                     = np.empty(int(NumMeas))
        Potenzialmaximumsarray     = np.empty(int(NumMeas))
        Potenzialhalbmaximumsarray = np.empty(int(NumMeas))
        

        for i in range(int(NumMeas)):
            
            StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
            UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
            if DesiReactxx == 0:
                Stromarrays                = StromarraysROH[::-1]
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::-1]
                Potenzialarrays            = PotenzialarraysROH[::-1]   -  Ru*Stromarrays/1000000
            
            
            if DesiReactxx == 1:
                Stromarrays                = StromarraysROH
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
                Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
            if BaseCorr == 1:
                BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
                if DesiReactxx == 0:
                    BaseStromArrays        = BaseStromArraysROH[::-1]
                if DesiReactxx == 1:
                    BaseStromArrays        = BaseStromArraysROH[::]
                Stromarrays            = Stromarrays - BaseStromArrays
            
            LenPotArrays            = len(Potenzialarrays)
        
            AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=LenPotArrays*REFF)#, endpoint=True)
            Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
            AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
            if Smooth_Data == 1:
            
                from scipy.interpolate import UnivariateSpline
            
                SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
                SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
                IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
                IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
                spl = UnivariateSpline(Potenzialarrays, Stromarrays)
                spl.set_smoothing_factor(SmoothFact)
                AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , REFF*len(Potenzialarrays))
                AuflgelStrarrays = spl(AufgelPotArrays)
            
            
            #___________________________________________________________
            #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
            #___________________________________________________________
            KorrAufgPotArr         = AufgelPotArrays #Korrektur war schon weiter oben
            
 
            KorrAufgStrArrays      = AuflgelStrarrays  
            if AutoCorr == 1:
                
                OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
                OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
                KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
             
            if ManuCorr == 1:
                KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr*ScanratesArray[i] #MauCorrCurr ist hier eine Kapazität von oben
                #Es ist nur noch als Curr benannt, weil ich zu faul zum Umdefinieren bin                                                                    
            
            for j in range(Length):
                
                KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
                KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
            for j in range(ShortLength):
                
                UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
                UNKORRSTROMSUPPERARRAY[i,j]     = Stromarrays[j]
            
            
            
            
            
            
            if DesiReactxx == 0:
                PeakMax    = np.min(KorrAufgStrArrays) 
            
            if DesiReactxx == 1:
                PeakMax    = np.max(KorrAufgStrArrays) 
            
            
    
            I_Peak[i]       = PeakMax
            
            #AB HIER NEUES
        
            IndexPeak                    = np.asscalar(np.where(KorrAufgStrArrays == PeakMax) [0]) 
            Peakpotenzial                = KorrAufgPotArr[IndexPeak]
            Potenzialmaximumsarray[i]    = Peakpotenzial
            PeakPotentiale[i]            = Peakpotenzial  #globales zum Rausnehmen
            HalbMaxAuflStrom             = PeakMax*0.5
            AuflStromBisPeak             = KorrAufgStrArrays[0:IndexPeak]
            HalbMaxAuflStromBisPeak      = (find_nearest(AuflStromBisPeak, HalbMaxAuflStrom))
            IndexHalbMaxAuflStromBisPeak = np.asscalar(np.where(AuflStromBisPeak ==    HalbMaxAuflStromBisPeak) [0])
            Halbpeakpotenzial            = (np.asscalar(KorrAufgPotArr[IndexHalbMaxAuflStromBisPeak]))
            Potenzialhalbmaximumsarray[i]= Halbpeakpotenzial
            HalbPeakPotenziale[i]        = Halbpeakpotenzial #globales zum Rausnehmen    
            
            #JETZT MÜSSTE ALLES WICHTIGE DEFINIERT SEIN
        
        
        
            bild1 = f.add_subplot(221)
            if ShowUncorrLSV ==1:
                bild1.plot(UnkorrPotenzialarrays,Stromarrays*0.001, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays*0.001, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [mA]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            bild1.axvline(Potenzialhalbmaximumsarray[i],color='r')
            bild1.axvline(Potenzialmaximumsarray[i],color='b')    
            
            
            if SeparateLSV == 1:
                LSVbild1 = SepLSVPlot.add_subplot(111)
                if ShowUncorrLSV ==1:
                    LSVbild1.plot(UnkorrPotenzialarrays,Stromarrays*0.001, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays*0.001, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [mA]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
                LSVbild1.axvline(Potenzialhalbmaximumsarray[i],color='r')
                LSVbild1.axvline(Potenzialmaximumsarray[i],color='b')    
        

        

        WurzelScanratesarray = (ScanratesArray*0.001)**0.5
        LnScanratesarray     = np.log(ScanratesArray*0.001)
        
        #FITFUNKTIONEN NOCH OHNE SCANRATELIMITIERUNG
        
        RSIRRFitting, pcov = curve_fit(GeradenFit, WurzelScanratesarray,I_Peak)  
        
        FITTED_I_Peak = GeradenFit(WurzelScanratesarray, *RSIRRFitting)
        
        #RSIRR-FIT Beta und IEQ
        #____________________________________________________________________________________________________________
        
        beta_RS = (RSIRRFitting[1]*0.000001*(R*T)**0.5 / (0.4967*n*F*A*c*(n*F*D)**0.5))**2
        global RS_BETA
        RS_BETA = beta_RS
        
        #Wenn Cottiref
        if RefToCotti == 1:
            beta_n_RS = np.absolute((RSIRRFitting[1]*(R*T)**0.5)/(np.pi**0.5 * 0.4967 * CottiSlope * F**0.5))**2      
            global nBETA_RS
            nBETA_RS = beta_n_RS  #globales Rausholen
            
        #IEQs Berechnen
        global RS_IEQ
        global RS_KNULL
        if DesiReactxx == 0:
            for i in range(int(NumMeas)):
                RS_IEQARRAY[i] = np.absolute(2.1815*n*F*c*(beta_RS*n*F*D*0.001*ScanratesArray[i]/(R*T))**0.5 * np.exp(beta_RS*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            RS_IEQ      = np.mean(RS_IEQARRAY)
            RS_IEQma    = RS_IEQ*1000000
            
            if KZeroCalcer ==1:
                RS_KNULLARRAY = (RS_IEQARRAY/(n*F*c))*np.exp(RS_BETA*n*F*(OCPArray-EZero)/(R*T))
                RS_KNULL      = np.mean(RS_KNULLARRAY)
            
        if DesiReactxx == 1:
            for i in range(int(NumMeas)):
                RS_IEQARRAY[i] = np.absolute(2.1815*n*F*c*(beta_RS*n*F*D*0.001*ScanratesArray[i]/(R*T))**0.5 * np.exp(-beta_RS*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            RS_IEQ      = np.mean(RS_IEQARRAY)
            RS_IEQma    = RS_IEQ*1000000
            
            if KZeroCalcer ==1:
                RS_KNULLARRAY = (RS_IEQARRAY/(n*F*c))*np.exp(-RS_BETA*n*F*(OCPArray-EZero)/(R*T))
                RS_KNULL      = np.mean(RS_KNULLARRAY)
   
        #___________________________________________________________________________________________________________
        
        
        #AYABE1-FIT Beta und IEQ
        #____________________________________________________________________________________________________________
        AyabeFit,pcov     = curve_fit(GeradenFit, LnScanratesarray,Potenzialmaximumsarray)
        betaVonAY = -(R*T/(2*n*F*AyabeFit[1])) 
        betaAy1 = np.absolute(betaVonAY)
        global LNPLOT_BETA
        LNPLOT_BETA = betaAy1 #globales Rausholen
        FIT_PeakPotentiale = GeradenFit(LnScanratesarray, *AyabeFit*1000)
        global LNPLOT_IEQ
        global LNPLOT_KNULL
        if DesiReactxx == 0:
            for i in range(int(NumMeas)):
                LNPLOT_IEQARRAY[i]= np.absolute(2.1815*n*F*c*(betaAy1*n*F*D*0.001*ScanratesArray[i]/(R*T))**0.5 * np.exp(betaAy1*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            LNPLOT_IEQ    = np.mean(LNPLOT_IEQARRAY)
            LNPLOT_IEQma  = LNPLOT_IEQ*1000000
            
            if KZeroCalcer ==1:
                LNPLOT_KNULLARRAY = (LNPLOT_IEQARRAY/(n*F*c))*np.exp(LNPLOT_BETA*n*F*(OCPArray-EZero)/(R*T))
                LNPLOT_KNULL    = np.mean(LNPLOT_KNULLARRAY)
        
        if DesiReactxx == 1:
            for i in range(int(NumMeas)):
                LNPLOT_IEQARRAY[i] = np.absolute(2.1815*n*F*c*(betaAy1*n*F*D*0.001*ScanratesArray[i]/(R*T))**0.5 * np.exp(-betaAy1*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            LNPLOT_IEQ    = np.mean(LNPLOT_IEQARRAY)
            LNPLOT_IEQma  = LNPLOT_IEQ*1000000
            
            if KZeroCalcer ==1:
                LNPLOT_KNULLARRAY = (LNPLOT_IEQARRAY/(n*F*c))*np.exp(-LNPLOT_BETA*n*F*(OCPArray-EZero)/(R*T))
                LNPLOT_KNULL    = np.mean(LNPLOT_KNULLARRAY)
        #_____________________________________________________________________________________________________________
        
        
        
        
        
        #AYABE2-Beta und IEQ
        #_____________________________________________________________________________________________________________
        Delta_EP_EPHalb = np.absolute(1.85*R*T/(F*n*(Potenzialhalbmaximumsarray - Potenzialmaximumsarray)))
        SHAPE_BETAARRAY = Delta_EP_EPHalb
        beta = np.absolute(np.mean(Delta_EP_EPHalb))
        betaAy2 = beta
        global SHAPE_BETA
        SHAPE_BETA = betaAy2 #globales Rausholen
        global SHAPE_IEQ
        global SHAPE_KNULL
        if DesiReactxx == 0:
            SHAPE_IEQARRAY = 2.1815*n*F*c*(Delta_EP_EPHalb*n*F*D*0.001*ScanratesArray/(R*T))**0.5 * np.exp(Delta_EP_EPHalb*n*F*(Potenzialmaximumsarray-OCPArray)/(R*T))
            SHAPE_IEQ      = np.mean(SHAPE_IEQARRAY)
            SHAPE_IEQma    = SHAPE_IEQ*1000000
 
            if KZeroCalcer ==1:
                SHAPE_KNULLARRAY = (SHAPE_IEQARRAY/(n*F*c))*np.exp(SHAPE_BETA*n*F*(OCPArray-EZero)/(R*T))
                SHAPE_KNULL    = np.mean(SHAPE_KNULLARRAY)
        
        if DesiReactxx == 1:
            SHAPE_IEQARRAY = 2.1815*n*F*c*(Delta_EP_EPHalb*n*F*D*0.001*ScanratesArray/(R*T))**0.5 * np.exp(-Delta_EP_EPHalb*n*F*(Potenzialmaximumsarray-OCPArray)/(R*T))
            SHAPE_IEQ      = np.mean(SHAPE_IEQARRAY)
            SHAPE_IEQma    = SHAPE_IEQ*1000000
            

            if KZeroCalcer ==1:
                SHAPE_KNULLARRAY = (SHAPE_IEQARRAY/(n*F*c))*np.exp(-SHAPE_BETA*n*F*(OCPArray-EZero)/(R*T))
                SHAPE_KNULL    = np.mean(SHAPE_KNULLARRAY)
            
        #_____________________________________________________________________________________________________________
        
        
        
        
        
        
        
        bild3 = f.add_subplot(223)
        bild3.plot(LnScanratesarray,1000*Potenzialmaximumsarray,linestyle='',marker='.')
        bild3.plot(LnScanratesarray,FIT_PeakPotentiale,linestyle='-',marker='')
        bild3.set_xlabel(r'$\ln{\nu} $'  , fontsize=12)
        bild3.set_ylabel('E$_p$' '-E$_{p/2}$' '[mV]' , fontsize=12)
        for axis in ['top','bottom','left','right']:
            bild3.spines[axis].set_linewidth(2)
        bild3.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        if DesiReactxx == 0:
            bild3.annotate(u'\u03B2=%5.3f'    % betaAy1, xy=(0.4, 0.93), xycoords='axes fraction')
            bild3.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % LNPLOT_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
            
        if DesiReactxx == 1:
            bild3.annotate(u'\u03B2=%5.3f'  % betaAy1, xy=(0.4, 0.13), xycoords='axes fraction')
            bild3.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % LNPLOT_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
        
        
        if SeparateMA1 ==1:
            MA1bild3 = MA1Plot.add_subplot(111)
            MA1bild3.plot(LnScanratesarray,1000*Potenzialmaximumsarray,linestyle='',marker='.')
            MA1bild3.plot(LnScanratesarray,FIT_PeakPotentiale,linestyle='-',marker='')
            MA1bild3.set_xlabel(r'$\ln{\nu} $'  , fontsize=12)
            MA1bild3.set_ylabel('E$_p$' '-E$_{p/2}$' '[mV]' , fontsize=12)
            for axis in ['top','bottom','left','right']:
                MA1bild3.spines[axis].set_linewidth(2)
            MA1bild3.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
            if DesiReactxx == 0:
                MA1bild3.annotate(u'\u03B2=%5.3f'    % betaAy1, xy=(0.4, 0.93), xycoords='axes fraction')
                MA1bild3.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % LNPLOT_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
            if DesiReactxx == 1:
                MA1bild3.annotate(u'\u03B2=%5.3f'  % betaAy1, xy=(0.4, 0.13), xycoords='axes fraction')
                MA1bild3.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'   % LNPLOT_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
            
            canvas = FigureCanvasTkAgg(MA1Plot, master=MA1Fenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, MA1Fenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        
        
        bild4  = f.add_subplot(224)
        bild4.plot(WurzelScanratesarray,Delta_EP_EPHalb,linestyle='-',marker='.')
        bild4.set_xlabel(r'$\sqrt{\nu} $' ' ' r'[$\sqrt{V} $]' , fontsize=12)
        bild4.set_ylabel(u'\u03B2' , fontsize=12)
        for axis in ['top','bottom','left','right']:
            bild4.spines[axis].set_linewidth(2)
        bild4.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        bild4.set_ylim([0,1])
        
        if DesiReactxx == 0:
            bild4.annotate(u'\u03B2=%5.3f'    % betaAy2, xy=(0.4, 0.93), xycoords='axes fraction') 
            bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
        if DesiReactxx == 1:
            bild4.annotate(u'\u03B2=%5.3f'  % betaAy2, xy=(0.4, 0.13), xycoords='axes fraction')  
            bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
        
        if SeparateMA2 ==1:
            MA2bild4  = MA2Plot.add_subplot(111)
            MA2bild4.plot(WurzelScanratesarray,Delta_EP_EPHalb,linestyle='-',marker='.')
            MA2bild4.set_xlabel(r'$\sqrt{\nu} $' ' ' r'[$\sqrt{V} $]' , fontsize=12)
            MA2bild4.set_ylabel(u'\u03B2' , fontsize=12)
            for axis in ['top','bottom','left','right']:
                MA2bild4.spines[axis].set_linewidth(2)
            MA2bild4.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            MA2bild4.set_ylim([0,1])
            
            if DesiReactxx == 0:
                MA2bild4.annotate(u'\u03B2=%5.3f'    % betaAy2, xy=(0.4, 0.93), xycoords='axes fraction')
                MA2bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
            if DesiReactxx == 1:
                MA2bild4.annotate(u'\u03B2=%5.3f'  % betaAy2, xy=(0.4, 0.13), xycoords='axes fraction')  
                MA2bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
            
            canvas = FigureCanvasTkAgg(MA2Plot, master=MA2Fenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, MA2Fenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        
        
        bild2 = f.add_subplot(222)
        
        
        
        
        def RSIRRSCPLOTTER():
        
            bild2.plot(WurzelScanratesarray,0.001*GeradenFit(WurzelScanratesarray, *RSIRRFitting),color='r',linestyle='-',marker='')
            bild2.plot(WurzelScanratesarray,0.001*I_Peak,linestyle='',marker='.')
            bild2.set_xlabel(r'$\sqrt{\nu} $' ' ' r'[$\sqrt{V} $]' , fontsize=12)
            bild2.set_ylabel('I$_p$ [mA]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
               
            
            #___________________________________________________________________________________________________________________
            
            if RefToCotti == 0:
            
                if DesiReactxx == 0:
                    bild2.annotate(u'\u03B2=%8.3f' % beta_RS, xy=(0.4, 0.93), xycoords='axes fraction')
                    bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
                
                if DesiReactxx == 1:
                    bild2.annotate(u'\u03B2=%8.3f' % beta_RS, xy=(0.4, 0.13), xycoords='axes fraction')
                    bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
            
            if RefToCotti == 1:
               
                BeschreibefensterCotti = Toplevel()                                                         
                BeschreibefensterCotti.title("Information")                         
                BeschreibefensterCotti.geometry("200x200")    
                Te = Text(BeschreibefensterCotti, height=6, width=30)
                Te.pack()
                Te.insert(END, "Randles Slope got\ndivided by Cottrell Slope")
                    

                if DesiReactxx == 0:
                    bild2.annotate('n'  u'\u03B2=%5.3f'    % beta_n_RS, xy=(0.4, 0.93), xycoords='axes fraction') 
                    bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
                if DesiReactxx == 1:
                    bild2.annotate('n'  u'\u03B2=%5.3f'  % beta_n_RS, xy=(0.4, 0.13), xycoords='axes fraction') 
                    bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
        

           
        
        RSIRRSCPLOTTER()

        f.tight_layout()
        
        canvas = FigureCanvasTkAgg(f, master=Fenster)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, Fenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        


        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #ONLY RSIRR PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyRSPlot == 1:
            
            RSIRRFenster = Toplevel()                                                         
            RSIRRFenster.title("Irreversible Randles-Sevcik Plot")                         
            RSIRRFenster.geometry("700x700")
        
            RSIRRPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = RSIRRPlot.add_subplot(111)
            
            RSIRRSCPLOTTER()
                
            canvas = FigureCanvasTkAgg(RSIRRPlot, master=RSIRRFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, RSIRRFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        #Das ist das Canvas für den Separaten LSV Plot, was vllt in der for schleife nicht geht???
        
        if SeparateLSV == 1:
            canvas = FigureCanvasTkAgg(SepLSVPlot, master=SepLSVFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, SepLSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
        if AsTxtSaver ==1:
            RSIRRSCASTXTsaver()
    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #FOR SCHLEIFE WENN CONCENTRATIONS GEÄNDERT WERDEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________    
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________   
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________
    #______________________________________________________________________________________________________________ 
    

    
    if VarParCon == 1:
    
        
        global SHAPE_BETAARRAYCONC
        global RS_IEQARRAYCONC
        global RS_KNULLARRAYCONC
        global SHAPE_IEQARRAYCONC
        global SHAPE_KNULLARRAYCONC

   
    
        PeakPotentiale     = np.empty(int(NumMeas))
        HalbPeakPotenziale = np.empty(int(NumMeas))
        SHAPE_BETAARRAYCONC= np.empty(int(NumMeas))
        RS_IEQARRAYCONC    = np.empty(int(NumMeas))
        
        
        OCPs = OCPArray
        
        #________________________________________________________________      
        
               
        I_Peak                     = np.empty(int(NumMeas))
        Potenzialmaximumsarray     = np.empty(int(NumMeas))
        Potenzialhalbmaximumsarray = np.empty(int(NumMeas))
        
        

        for i in range(int(NumMeas)):
            
            StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
            UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
            
            
            if DesiReactxx == 0:
                Stromarrays                = StromarraysROH[::-1]
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::-1]
                Potenzialarrays            = PotenzialarraysROH[::-1]   -  Ru*Stromarrays/1000000
            
            
            if DesiReactxx == 1:
                Stromarrays                = StromarraysROH
                UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH
                Potenzialarrays            = PotenzialarraysROH   -  Ru*Stromarrays/1000000
            
            
            
            if BaseCorr == 1:
                BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
                if DesiReactxx == 0:
                    BaseStromArrays        = BaseStromArraysROH[::-1]
                if DesiReactxx == 1:
                    BaseStromArrays        = BaseStromArraysROH[::]
                Stromarrays            = Stromarrays - BaseStromArrays
            
            
            LenPotArrays            = len(Potenzialarrays)
        
            AufgelPotArrays        = np.linspace(Potenzialarrays[0], Potenzialarrays[-1], num=LenPotArrays*REFF)#, endpoint=True)
            Interpolation          = interp1d(Potenzialarrays, Stromarrays, kind='cubic')                  
            AuflgelStrarrays       = Interpolation(AufgelPotArrays)                                        
        
            if Smooth_Data == 1:
            
                from scipy.interpolate import UnivariateSpline
            
                SmoothFromPotential = EchtesPotFinden(Potenzialarrays,SmoothFrom)
                SmoothToPotential   = EchtesPotFinden(Potenzialarrays,SmoothTo)
            
                IdxSmoothFromPotential = np.asscalar(np.where(Potenzialarrays == SmoothFromPotential) [0])
                IdxSmoothToPotential   = np.asscalar(np.where(Potenzialarrays == SmoothToPotential) [0])
            
            
                spl = UnivariateSpline(Potenzialarrays, Stromarrays)
                spl.set_smoothing_factor(SmoothFact)
                AufgelPotArrays  = np.linspace(Potenzialarrays[IdxSmoothFromPotential],Potenzialarrays[IdxSmoothToPotential] , REFF*len(Potenzialarrays))
                AuflgelStrarrays = spl(AufgelPotArrays)
            
            
            #___________________________________________________________
            #AUFGELÖSTES POTENZIALARRAY UM OHMSCHEN ANTEIL KORRIGIEREN
            #___________________________________________________________
            KorrAufgPotArr         = AufgelPotArrays #Korrektur war schon weiter oben
            
 
            KorrAufgStrArrays      = AuflgelStrarrays  
            if AutoCorr == 1:
                
                OCPotenziale       = EchtesPotFinden(KorrAufgPotArr,OCPs[i])
                OCPIndices         = np.asscalar(np.where(KorrAufgPotArr == OCPotenziale) [0])
                KorrAufgStrArrays  = AuflgelStrarrays - AuflgelStrarrays[OCPIndices]
             
            if ManuCorr == 1:
                KorrAufgStrArrays  = AuflgelStrarrays - ManuCorrCurr*Scanrate#MauCorrCurr ist hier eine Kapazität von oben
                #Es ist nur noch als Curr benannt, weil ich zu faul zum Umdefinieren bin 
                
            
            
            for j in range(Length):
                
                KORRPOTENZIALSUPPERARRAY[i,j] = KorrAufgPotArr[j] 
                KORRSTROMSUPPERARRAY[i,j]     = KorrAufgStrArrays[j]
                      
            for j in range(ShortLength):
                
                UNKORRPOTENZIALSUPPERARRAY[i,j] = UnkorrPotenzialarrays[j] 
                UNKORRSTROMSUPPERARRAY[i,j]     = Stromarrays[j]
            
 
            
            
            if DesiReactxx == 0:
                PeakMax    = np.min(KorrAufgStrArrays) 
            
            if DesiReactxx == 1:
                PeakMax    = np.max(KorrAufgStrArrays) 
            
    
            I_Peak[i]       = PeakMax
            
            
            #AB HIER NEUES
        
            IndexPeak                    = np.asscalar(np.where(KorrAufgStrArrays == PeakMax) [0]) 
            Peakpotenzial                = KorrAufgPotArr[IndexPeak]
            Potenzialmaximumsarray[i]    = Peakpotenzial
            PeakPotentiale[i]            = Peakpotenzial  #globales zum Rausnehmen
            HalbMaxAuflStrom             = PeakMax*0.5
            AuflStromBisPeak             = KorrAufgStrArrays[0:IndexPeak]
            HalbMaxAuflStromBisPeak      = (find_nearest(AuflStromBisPeak, HalbMaxAuflStrom))
            IndexHalbMaxAuflStromBisPeak = np.asscalar(np.where(AuflStromBisPeak ==    HalbMaxAuflStromBisPeak) [0])
            Halbpeakpotenzial            = (np.asscalar(KorrAufgPotArr[IndexHalbMaxAuflStromBisPeak]))
            Potenzialhalbmaximumsarray[i]= Halbpeakpotenzial
            HalbPeakPotenziale[i]        = Halbpeakpotenzial #globales zum Rausnehmen    
        
            #JETZT MÜSSTE ALLES WICHTIGE DEFINIERT SEIN
        
        
        
            bild1 = f.add_subplot(221)
            if ShowUncorrLSV ==1:
                bild1.plot(UnkorrPotenzialarrays,Stromarrays*0.001, linestyle='-',marker='',color='k')
            bild1.plot(KorrAufgPotArr,KorrAufgStrArrays*0.001, linestyle='-',marker='', color='r')
            bild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
            bild1.set_ylabel('I [mA]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild1.spines[axis].set_linewidth(2)
            bild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            bild1.axvline(OCPs[i],color='k')
            bild1.axhline(0, color='k')
            
            if SeparateLSV == 1:
                LSVbild1 = SepLSVPlot.add_subplot(111)
                if ShowUncorrLSV ==1:
                    LSVbild1.plot(UnkorrPotenzialarrays,Stromarrays*0.001, linestyle='-',marker='',color='k')
                LSVbild1.plot(KorrAufgPotArr,KorrAufgStrArrays*0.001, linestyle='-',marker='', color='r')
                LSVbild1.set_xlabel('E vs. Ref. [V]', fontsize=12)
                LSVbild1.set_ylabel('I [mA]', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    LSVbild1.spines[axis].set_linewidth(2)
                LSVbild1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
                LSVbild1.axvline(OCPs[i],color='k')
                LSVbild1.axhline(0, color='k')
                

        global CDecayArray

        CDecayArray = ConcentrationsArray*c
       

        #FITFUNKTIONEN NOCH OHNE SCANRATELIMITIERUNG
        
        RSIRRFitting, pcov = curve_fit(GeradenFit, CDecayArray,I_Peak)  
         
        
        FITTED_I_Peak = GeradenFit(CDecayArray, *RSIRRFitting)
        
        #RSIRR-FIT Beta und IEQ
        #____________________________________________________________________________________________________________
        
        beta_RS = (RSIRRFitting[1]*0.000001*(R*T)**0.5 / (0.4967*n*F*A*(n*F*D*Scanrate*0.001)**0.5))**2
        global RS_BETACONC
        RS_BETACONC = beta_RS
        
        #HIER GEHT KEIN REF TO COTTI
            
        #IEQs Berechnen
        global RS_IEQCONC
        global RS_KNULLCONC
        if DesiReactxx == 0:
            for i in range(int(NumMeas)):
                RS_IEQARRAYCONC[i] = np.absolute(2.1815*n*F*CDecayArray[i]*(beta_RS*n*F*D*0.001*Scanrate/(R*T))**0.5 * np.exp(beta_RS*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            RS_IEQCONC      = np.mean(RS_IEQARRAYCONC)
            RS_IEQma    = RS_IEQCONC*1000000
            
            if KZeroCalcer ==1:
                RS_KNULLARRAYCONC = (RS_IEQARRAYCONC/(n*F*CDecayArray[i]))*np.exp(RS_BETACONC*n*F*(OCPArray-EZero)/(R*T))
                RS_KNULLCONC      = np.mean(RS_KNULLARRAYCONC)
            
        if DesiReactxx == 1:
            for i in range(int(NumMeas)):
                RS_IEQARRAYCONC[i] = np.absolute(2.1815*n*F*CDecayArray[i]*(beta_RS*n*F*D*0.001*Scanrate/(R*T))**0.5 * np.exp(-beta_RS*n*F*(Potenzialmaximumsarray[i]-OCPArray[i])/(R*T)))
            RS_IEQCONC      = np.mean(RS_IEQARRAYCONC)
            RS_IEQma    = RS_IEQCONC*1000000
            
            if KZeroCalcer ==1:
                RS_KNULLARRAYCONC = (RS_IEQARRAYCONC/(n*F*CDecayArray[i]))*np.exp(-RS_BETACONC*n*F*(OCPArray-EZero)/(R*T))
                RS_KNULLCONC      = np.mean(RS_KNULLARRAYCONC)
   
        #___________________________________________________________________________________________________________
        
        
        #AYABE2-Beta und IEQ
        #_____________________________________________________________________________________________________________
        Delta_EP_EPHalb = np.absolute(1.85*R*T/(F*n*(Potenzialhalbmaximumsarray - Potenzialmaximumsarray)))
        SHAPE_BETAARRAYCONC = Delta_EP_EPHalb
        beta = np.absolute(np.mean(Delta_EP_EPHalb))
        betaAy2 = beta
        global SHAPE_BETACONC
        SHAPE_BETACONC = betaAy2 #globales Rausholen
        global SHAPE_IEQCONC
        global SHAPE_KNULLCONC
        if DesiReactxx == 0:
            SHAPE_IEQARRAYCONC = 2.1815*n*F*CDecayArray*(Delta_EP_EPHalb*n*F*D*0.001*Scanrate/(R*T))**0.5 * np.exp(Delta_EP_EPHalb*n*F*(Potenzialmaximumsarray-OCPArray)/(R*T))
            SHAPE_IEQCONC      = np.mean(SHAPE_IEQARRAYCONC)
            SHAPE_IEQma    = SHAPE_IEQCONC*1000000
 
            if KZeroCalcer ==1:
                SHAPE_KNULLARRAYCONC = (SHAPE_IEQARRAYCONC/(n*F*CDecayArray))*np.exp(SHAPE_BETACONC*n*F*(OCPArray-EZero)/(R*T))
                SHAPE_KNULLCONC    = np.mean(SHAPE_KNULLARRAYCONC)
        
        if DesiReactxx == 1:
            SHAPE_IEQARRAYCONC = 2.1815*n*F*CDecayArray*(Delta_EP_EPHalb*n*F*D*0.001*Scanrate/(R*T))**0.5 * np.exp(-Delta_EP_EPHalb*n*F*(Potenzialmaximumsarray-OCPArray)/(R*T))
            SHAPE_IEQCONC      = np.mean(SHAPE_IEQARRAYCONC)
            SHAPE_IEQma    = SHAPE_IEQCONC*1000000
            

            if KZeroCalcer ==1:
                SHAPE_KNULLARRAYCONC = (SHAPE_IEQARRAYCONC/(n*F*CDecayArray))*np.exp(-SHAPE_BETACONC*n*F*(OCPArray-EZero)/(R*T))
                SHAPE_KNULLCONC    = np.mean(SHAPE_KNULLARRAYCONC)
            
        #MATSUDA AYABE 1 --> GEHT HIER NICHT
        #_____________________________________________________________________________________________________________
        
        
        bild3 = f.add_subplot(223)
        bild3.annotate('Matsuda-Ayabe1 NOT possible here' , xy=(0.1, 0.93), xycoords='axes fraction')
        for axis in ['top','bottom','left','right']:
            bild3.spines[axis].set_linewidth(2)
        bild3.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
        
        if SeparateMA1 ==1:
            MA1bild3 = MA1Plot.add_subplot(111)
            MA1bild3.annotate('Matsuda-Ayabe1 NOT possible here', xy=(0.1, 0.93), xycoords='axes fraction')
            for axis in ['top','bottom','left','right']:
                MA1bild3.spines[axis].set_linewidth(2)
            MA1bild3.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            
            canvas = FigureCanvasTkAgg(MA1Plot, master=MA1Fenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, MA1Fenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        #_______________________________________________________________________________________________________________    
        
        
                           
        #MATSUDA AYABE 2
        #_________________________________________________________________________________________________________________
                           
        
        bild4  = f.add_subplot(224)
        bild4.plot(ConcentrationsArray,Delta_EP_EPHalb,linestyle='-',marker='.')
        bild4.set_xlabel('c/''c$_{max}$' , fontsize=12)
        bild4.set_ylabel(u'\u03B2' , fontsize=12)
        for axis in ['top','bottom','left','right']:
            bild4.spines[axis].set_linewidth(2)
        bild4.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        bild4.set_ylim([0,1])
        
        if DesiReactxx == 0:
            bild4.annotate(u'\u03B2=%5.3f'    % betaAy2, xy=(0.4, 0.93), xycoords='axes fraction') 
            bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
        if DesiReactxx == 1:
            bild4.annotate(u'\u03B2=%5.3f'  % betaAy2, xy=(0.4, 0.13), xycoords='axes fraction')  
            bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
        
        if SeparateMA2 ==1:
            MA2bild4  = MA2Plot.add_subplot(111)
            MA2bild4.plot(ConcentrationsArray,Delta_EP_EPHalb,linestyle='-',marker='.')
            MA2bild4.set_xlabel('c/''c$_{max}$' , fontsize=12)
            MA2bild4.set_ylabel(u'\u03B2' , fontsize=12)
            for axis in ['top','bottom','left','right']:
                MA2bild4.spines[axis].set_linewidth(2)
            MA2bild4.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            MA2bild4.set_ylim([0,1])
            
            if DesiReactxx == 0:
                MA2bild4.annotate(u'\u03B2=%5.3f'    % betaAy2, xy=(0.4, 0.93), xycoords='axes fraction')
                MA2bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
            if DesiReactxx == 1:
                MA2bild4.annotate(u'\u03B2=%5.3f'  % betaAy2, xy=(0.4, 0.13), xycoords='axes fraction')  
                MA2bild4.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % SHAPE_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
            
            canvas = FigureCanvasTkAgg(MA2Plot, master=MA2Fenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, MA2Fenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        

        
        
        bild2 = f.add_subplot(222)
        
        def RSIRRCONCPLOTTER():
        
            bild2.plot(ConcentrationsArray,FITTED_I_Peak,color='r',linestyle='-',marker='')
            bild2.plot(ConcentrationsArray,I_Peak,linestyle='',marker='.')
            bild2.set_xlabel('c/''c$_{max}$' , fontsize=12)
            bild2.set_ylabel('I$_p$ [mA]', fontsize=12)
            for axis in ['top','bottom','left','right']:
                bild2.spines[axis].set_linewidth(2)
            bild2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
            #________________________________________________________
            #________________________________________________________
            #ANNOTOATIONS FOR PARAMETERS AND CALCULATIONS
            #________________________________________________________
            #________________________________________________________
            #__________________________________________________________________________________________________________________
         
            beta_RS = (np.absolute(RSIRRFitting[1]*0.000001*(R*T)**0.5 / (0.4967*n*F*A*(Scanrate*0.001)**0.5 *(n*F*D)**0.5)))**2
            

            if DesiReactxx == 0:
                bild2.annotate(u'\u03B2=%8.3f' % beta_RS, xy=(0.4, 0.93), xycoords='axes fraction')
                bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.85), xycoords='axes fraction')
                
            if DesiReactxx == 1:
                bild2.annotate(u'\u03B2=%8.3f' % beta_RS, xy=(0.4, 0.13), xycoords='axes fraction')
                bild2.annotate('I$_{eq,Av}$=%5.3f*10$^{-6} Acm^{-2}$'    % RS_IEQma, xy=(0.4, 0.25), xycoords='axes fraction')
            
            
           
        RSIRRCONCPLOTTER()

        f.tight_layout()
        
        canvas = FigureCanvasTkAgg(f, master=Fenster)
        canvas.show()
        canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

        toolbar = NavigationToolbar2TkAgg(canvas, Fenster)
        toolbar.update()
        canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        


        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #ONLY RSIRR PLOT OPTION
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        #________________________________________________________________________________________________________________
        
        
        if OnlyRSPlot == 1:
            
            RSIRRFenster = Toplevel()                                                         
            RSIRRFenster.title("Irreversible Randles-Sevcik Plot")                         
            RSIRRFenster.geometry("700x700")
        
            RSIRRPlot = Figure(figsize=(6, 6), dpi=100)
            
            bild2 = RSIRRPlot.add_subplot(111)
            
            RSIRRCONCPLOTTER()
                
            canvas = FigureCanvasTkAgg(RSIRRPlot, master=RSIRRFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, RSIRRFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
        
        
        #Das ist das Canvas für den Separaten LSV Plot, was vllt in der for schleife nicht geht???
        
        if SeparateLSV == 1:
            canvas = FigureCanvasTkAgg(SepLSVPlot, master=SepLSVFenster)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

            toolbar = NavigationToolbar2TkAgg(canvas, SepLSVFenster)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
            
        if AsTxtSaver ==1:
            RSIRRCONCASTXTsaver()
    
    


# In[21]:


def KLROTASTXTsaver():

 
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        if get_n == 1:
            f.write("n-KL")
            f.write("\t")
            f.write(str(np.asscalar(n_KL)))
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("c in mol/cm^3")
            f.write("\t")
            f.write(str(c))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base")
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
            if RefToCotti  ==1:
                f.write("D ref. to Cottrell")
                f.write("\n")
                
            
  

            
        if get_c == 1:
            f.write("c-KL")
            f.write("\t")
            f.write(str(0.000001*np.asscalar(c_KL)))
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
                
                
                
        if get_D == 1:
            f.write("D-KL")  
            f.write("\t")
            f.write(str(0.000001*np.asscalar(D_KL)))
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("c in mol/cm^3")
            f.write("\t")
            f.write(str(c))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
        
        
        
            
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("KL-Plot Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("RotRates in rpm")   
        f.write("\t")
        f.write("I-Lim in MicAmp")   
        f.write("\t")
        f.write("Inv root RotRates")   
        f.write("\t")
        f.write("Inv I Lim ")   
        f.write("\t")
        f.write("Fit Inv I Lim")   
        f.write("\t")
        f.write("\n")
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(RotRatesArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(ILimMittel[i])))
            f.write("\t")
            f.write(str(np.asscalar(InvWurzelRot[i])))
            f.write("\t")
            f.write(str(np.asscalar(InvILimMittel[i])))
            f.write("\t")
            f.write(str(np.asscalar(InvILimFITTED[i])))
            f.write("\t")
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Rot-Speed")
        f.write("\t")
    
        for jj in range(len(RotRatesArray)):
        
            f.write(str(np.asscalar(RotRatesArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(len(RotRatesArray)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(len(RotRatesArray)):
                b = str(np.asscalar(KORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(KORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")
#__________________________________________________________________________________________________-        
        for i in range(3):
            f.write("\n")
            
        f.write("Uncorrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Rot-Speed")
        f.write("\t")
    
        for jj in range(len(RotRatesArray)):
        
            f.write(str(np.asscalar(RotRatesArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(len(RotRatesArray)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(ShortLength)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(len(RotRatesArray)):
                b = str(np.asscalar(UNKORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(UNKORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")

    root.destroy()




# In[22]:


def KLCONASTXTsaver():

 
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        if get_n == 1:
            f.write("n-KL")
            f.write("\t")
            f.write(str(np.asscalar(n_KL)))
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("Rotation speed in rpm")
            f.write("\t")
            f.write(str(RotRate))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base")
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
            if RefToCotti  ==1:
                f.write("D ref. to Cottrell")
                f.write("\n")
                
            
  

            
        if get_c == 1:
            f.write("c-KL")
            f.write("\t")
            f.write("Getting c not possible here")
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("Rotation speed in rpm")
            f.write("\t")
            f.write(str(RotRate))
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
                
                
                
        if get_D == 1:
            f.write("D-KL")  
            f.write("\t")
            f.write(str(0.000001*np.asscalar(D_KL)))
            f.write("\n")
            f.write("KL_Offset")
            f.write("\t")
            f.write(str(np.asscalar(KL_Offset)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("Rotation speed in rpm")
            f.write("\t")
            f.write(str(RotRate))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Visc in cm^2/s")
            f.write("\t")
            f.write(str(Visc))
            f.write("\n")
            f.write("KL_StartPot in V")
            f.write("\t")
            f.write(str(KLStart))
            f.write("\n")
            f.write("KL_EndPot in V")
            f.write("\t")
            f.write(str(KLEnd))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")

            
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("KL-Plot Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("Concentrations in mol/cm^3")   
        f.write("\t")
        f.write("I-Lim in MicAmp")   
        f.write("\t")
        f.write("Inv Concentrations")   
        f.write("\t")
        f.write("Inv I Lim ")   
        f.write("\t")
        f.write("Fit Inv I Lim")   
        f.write("\t")
        f.write("\n")
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(c*ConcentrationsArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(ILimMittel[i])))
            f.write("\t")
            f.write(str(np.asscalar(1/c*ConcentrationsArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(InvILimMittel[i])))
            f.write("\t")
            f.write(str(np.asscalar(InvILimFITTED[i])))
            f.write("\t")
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Concentrations as Fract of highest")
        f.write("\t")
    
        for jj in range(len(ConcentrationsArray)):
        
            f.write(str(np.asscalar(ConcentrationsArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(len(ConcentrationsArray)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(len(ConcentrationsArray)):
                b = str(np.asscalar(KORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(KORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")
#__________________________________________________________________________________________________-        
        for i in range(3):
            f.write("\n")
            
        f.write("Uncorrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Concentrations as Fract of highest")
        f.write("\t")
    
        for jj in range(len(ConcentrationsArray)):
        
            f.write(str(np.asscalar(ConcentrationsArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(len(ConcentrationsArray)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(ShortLength)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(len(ConcentrationsArray)):
                b = str(np.asscalar(UNKORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(UNKORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")

    root.destroy()




# In[23]:


def TAFELASTXTsaver():

 
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        

        f.write("Beta-Tafel-Average")
        f.write("\t")
        f.write(str(np.asscalar(BETA_T)))
        f.write("\n")
        f.write("Ieq-Tafel-Average in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(IEQT)))
        f.write("\n")
        if KZeroCalcer ==1:
            f.write("k_zero Tafel in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(KNULLT)))
            f.write("\n")
        f.write("Beta-Shape-Average")
        f.write("\t")
        f.write(str(np.asscalar(BETA_F)))
        f.write("\n")
        f.write("Ieq-Shape-Average A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(IEQF)))
        f.write("\n")
        if KZeroCalcer ==1:
            f.write("k_zero Shape in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(KNULLF)))
            f.write("\n")
        f.write("\n")
  
        f.write("Parameters")
        f.write("\n")
        f.write("A in cm^2")
        f.write("\t")
        f.write(str(A))
        f.write("\n")
        f.write("T in K")
        f.write("\t")
        f.write(str(T))
        f.write("\n")
        if VarParCon ==1:
            f.write("Rotation speed in rpm")
            f.write("\t")
            f.write(str(RotRate))
            f.write("\n")
        

        f.write("Tafel region in V")
        f.write("\t")
        f.write(str(Tafelreg))
        f.write("\n")
        f.write("Pot of Lim-Current in V")
        f.write("\t")
        f.write(str(ManLimPot))
        f.write("\n")
        f.write("Shift of ideal Tafelstart in V")
        f.write("\t")
        f.write(str(TKSTP))
        f.write("\n")
        f.write("Refining Factor")
        f.write("\t")
        f.write(str(REFF))
        f.write("\n")
        f.write("Ohmic Resistance in Ohm")
        f.write("\t")
        f.write(str(Ru))
        f.write("\n")
        if VarParCon ==1:
            f.write("Rotation Rate")
            f.write("\t")
            f.write(str(RotRate))
            f.write("\n")
            
        f.write("Correction Type")
        f.write("\t")
        if AutoCorr == 1:
            f.write("Auto")
            f.write("\n")
        if BaseCorr == 1:
            f.write("Base")
            f.write("\n")
        if ManuCorr == 1:
            f.write("Manual by")
            f.write("\t")
            f.write(str(ManuCorrCurr))
            f.write("\n")
        if Smooth_Data == 1:
            f.write("Data got smoothed")
            f.write("\n")
            f.write("Smoothing factor")
            f.write("\t")
            f.write(str(SmoothFact))
            f.write("\n")
            f.write("Smoothing Startpotential in V")
            f.write("\t")
            f.write(str(SmoothFrom))
            f.write("\n")
            f.write("Smoothing Endpotential in V")
            f.write("\t")
            f.write(str(SmoothTo))
            f.write("\n")
            



            

#______________________________________________________________________________________________________________________
#
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("Individual Tafel-Plot Results")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        if VarParRot ==1:
            f.write("Rotation Rates in rpm")   
            f.write("\t")
        if VarParCon ==1:
            f.write("Concentrations in mol/cm^3")   
            f.write("\t")
        f.write("Individual Tafelbeta")   
        f.write("\t")
        f.write("Individual Tafel-Ieq in A/cm^2")   
        f.write("\t")
        if KZeroCalcer ==1:
            f.write("Individual k_zero Tafel in cm/s")   
            f.write("\t")
        f.write("Individual Shapebeta")   
        f.write("\t")
        f.write("Individual Shape-Ieq in A/cm^2")   
        f.write("\t")
        if KZeroCalcer ==1:
            f.write("Individual k_zero Shape in cm/s")   
            f.write("\t")
        f.write("Potential of 0.5max Current in V")   
        f.write("\t")
        f.write("Potential of 0.75max Current in V")   
        f.write("\t")
        f.write("Fit-Startpot. in Tafelplot in V")   
        f.write("\t")
        f.write("Fit-Endpot. in Tafelplot in V")   
        f.write("\t")
        f.write("Tafel-Startpot in LSV-Plot in V")   
        f.write("\t")
        f.write("Tafel-Endpot in LSV-Plot in V")   
        f.write("\t")
        f.write("\n")
        f.write("\n")

        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            if VarParRot ==1:
                f.write(str(np.asscalar(RotRatesArray[i])))
                f.write("\t")
            if VarParCon ==1:
                f.write(str(np.asscalar(c*ConcentrationsArray[i])))
                f.write("\t")   
            f.write(str(np.asscalar(TafelbetaArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(TafelIEQArray[i])))
            f.write("\t")
            if KZeroCalcer ==1:
                f.write(str(np.asscalar(KNULLTARRAY[i])))
                f.write("\t")
            f.write(str(np.asscalar(FormBetaArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(FormIEQArray[i])))
            f.write("\t")
            if KZeroCalcer ==1:
                f.write(str(np.asscalar(KNULLFARRAY[i])))
                f.write("\t")
            f.write(str(np.asscalar(PotVonHalbLim[i])))
            f.write("\t")
            f.write(str(np.asscalar(PotVonDreivLim[i])))
            f.write("\t")
            f.write(str(np.asscalar(EtaStartEchtARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(EtaEndEchtARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(TafelStartImLSVARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(TafelEndImLSVARRAY[i])))
            f.write("\t")
            f.write("\n")
            

            

        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        
        if VarParRot == 1:
            f.write("Rotation Rates in rpm")
            f.write("\t")
            for jj in range(int(NumMeas)):
                f.write(str(np.asscalar(RotRatesArray[jj])))
                f.write("\t")
                f.write("\t")
            f.write("\n")
            f.write("\t")
        
        if VarParCon == 1:
            f.write("Concentrations as Fract of highest")
            f.write("\t")
            for jj in range(int(NumMeas)):
                f.write(str(np.asscalar(ConcentrationsArray[jj])))
                f.write("\t")
                f.write("\t")
            f.write("\n")
            f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(KORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(KORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")
#__________________________________________________________________________________________________-        
        for i in range(3):
            f.write("\n")
            
        f.write("Uncorrected LSVs")
        f.write("\n")
        f.write("\n")
        if VarParRot ==1:
            f.write("Rotation rates in rpm")
            f.write("\t")
        if VarParCon ==1:
            f.write("Concentrations as Fract of highest")
            f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            if VarParRot ==1:
                f.write(str(np.asscalar(RotRatesArray[jj])))
                f.write("\t")
                f.write("\t")
            
            
            if VarParCon ==1:
                f.write(str(np.asscalar(ConcentrationsArray[jj])))
                f.write("\t")
                f.write("\t")
        
        f.write("\n")
        f.write("\t")
     
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        for i in range(3):
            f.write("\n")
    
        for i in range(int(ShortLength)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(UNKORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(UNKORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")

#_________________________________________                
                
        f.write("Tafel Plots Data")
        f.write("\n")
        f.write("\n")
        
        if VarParRot == 1:
            f.write("Rotation Rates in rpm")
            f.write("\t")
            for jj in range(int(NumMeas)):
                f.write(str(np.asscalar(RotRatesArray[jj])))
                f.write("\t")
                f.write("\t")
                f.write("\t")
            f.write("\n")
            f.write("\t")
        
        if VarParCon == 1:
            f.write("Concentrations as Fract of highest")
            f.write("\t")
            for jj in range(int(NumMeas)):
                f.write(str(np.asscalar(ConcentrationsArray[jj])))
                f.write("\t")
                f.write("\t")
                f.write("\t")
            f.write("\n")
            f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("Eta in V")
            f.write("\t")
            f.write("LN(I*Ilim/(I-ILim))")
            f.write("\t")
            f.write("Tafelfit LN(I*Ilim/(I-ILim))")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(EtaArray[j:j+1:,i:i+1:]))
                d = str(np.asscalar(LogKinStromArray[j:j+1:,i:i+1:]))
                e = str(np.asscalar(LogFitStromArray[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
                f.write(e)
                f.write("\t")
            f.write("\n")
    

    
    root.destroy()




# In[24]:


def RSREVSCASTXTsaver():
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        if get_n == 1:
            f.write("n-RS-Rev")
            f.write("\t")
            f.write(str(np.asscalar(n_RS)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("c in mol/cm^3")
            f.write("\t")
            f.write(str(c))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base")
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
            if RefToCotti  ==1:
                f.write("D ref. to Cottrell")
                f.write("\n")
                
            
  

            
        if get_c == 1:
            f.write("c-RS-REV")
            f.write("\t")
            f.write(str(0.000001*np.asscalar(c_RS)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
                
                
                
        if get_D == 1:
            f.write("D-RS-REV")  
            f.write("\t")
            f.write(str(0.000001*np.asscalar(D_RS)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("c in mol/cm^3")
            f.write("\t")
            f.write(str(c))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
        
        
        
            
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("RS-REV-Plot Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("Scanrates in mV/s")   
        f.write("\t")
        f.write("Root of Scanrates in (mV/s)^0.5")   
        f.write("\t")
        f.write("I-Peak in microampere")   
        f.write("\t")
        f.write("Fit I-Peak in microampere")   
        f.write("\t")
        f.write("\n")
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(ScanratesArray[i])))
            f.write("\t")
            f.write(str(np.asscalar((ScanratesArray[i])**0.5)))
            f.write("\t")
            f.write(str(np.asscalar(I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(FITTED_I_Peak[i])))
            f.write("\t")
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Scanrate in mV/s")
        f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            f.write(str(np.asscalar(ScanratesArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(KORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(KORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")
#__________________________________________________________________________________________________-        
        for i in range(3):
            f.write("\n")
            
        f.write("Uncorrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Scanrate in mV/s")
        f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            f.write(str(np.asscalar(ScanratesArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(ShortLength)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(UNKORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(UNKORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")

    root.destroy()



 


# In[25]:


def RSREVCONCASTXTsaver():
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        if get_n == 1:
            f.write("n-RS-Rev")
            f.write("\t")
            f.write(str(np.asscalar(n_RS)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("Scanrate in mV/s")
            f.write("\t")
            f.write(str(Scanrate))
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("c_max in mol/cm^3")
            f.write("\t")
            f.write(str(c))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base")
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
            if RefToCotti  ==1:
                f.write("D ref. to Cottrell")
                f.write("\n")
                
            
  

            
        if get_c == 1:
            f.write("c-RS-REV")
            f.write("\t")
            f.write("getting c not possible here")
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("Scanrate in mV/s")
            f.write("\t")
            f.write(str(Scanrate))
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("D in cm^2/s")
            f.write("\t")
            f.write(str(D))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
                
                
                
        if get_D == 1:
            f.write("D-RS-REV")  
            f.write("\t")
            f.write(str(0.000001*np.asscalar(D_RS)))
            f.write("\n")
            f.write("\n")
            f.write("Parameters")
            f.write("\n")
            f.write("Scanrate in mV/s")
            f.write("\t")
            f.write(str(Scanrate))
            f.write("\n")
            f.write("given n")
            f.write("\t")
            f.write(str(n))
            f.write("\n")
            f.write("A in cm^2")
            f.write("\t")
            f.write(str(A))
            f.write("\n")
            f.write("T in K")
            f.write("\t")
            f.write(str(T))
            f.write("\n")
            f.write("c in mol/cm^3")
            f.write("\t")
            f.write(str(c))
            f.write("\n")
            f.write("Refining Factor")
            f.write("\t")
            f.write(str(REFF))
            f.write("\n")
            f.write("Ohmic Resistance in Ohm")
            f.write("\t")
            f.write(str(Ru))
            f.write("\n")
            f.write("Correction Type")
            f.write("\t")
            if AutoCorr == 1:
                f.write("Auto")
                f.write("\n")
            if BaseCorr == 1:
                f.write("Base") 
                f.write("\n")
            if ManuCorr == 1:
                f.write("Manual by")
                f.write("\t")
                f.write(str(ManuCorrCurr))
                f.write("\n")
            if Smooth_Data == 1:
                f.write("Data got smoothed")
                f.write("\n")
                f.write("Smoothing factor")
                f.write("\t")
                f.write(str(SmoothFact))
                f.write("\n")
                f.write("Smoothing Startpotential in V")
                f.write("\t")
                f.write(str(SmoothFrom))
                f.write("\n")
                f.write("Smoothing Endpotential in V")
                f.write("\t")
                f.write(str(SmoothTo))
                f.write("\n")
        
        
        
            
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("RS-REV-Plot Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("Concentrations in mol/cm^3")   
        f.write("\t")
        f.write("I-Peak in microampere")   
        f.write("\t")
        f.write("Fit I-Peak in microampere")   
        f.write("\t")
        f.write("\n")
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(c*ConcentrationsArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(FITTED_I_Peak[i])))
            f.write("\t")
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Concentr. as Fract. of highest")
        f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            f.write(str(np.asscalar(ConcentrationsArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(KORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(KORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")
#__________________________________________________________________________________________________-        
        for i in range(3):
            f.write("\n")
            
        f.write("Uncorrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Concentr. as Fract. of highest")
        f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            f.write(str(np.asscalar(ConcentrationsArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(ShortLength)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(UNKORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(UNKORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")

    root.destroy()



 


# In[26]:


def RSIRRSCASTXTsaver():
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")

        f.write("Beta-RS")
        f.write("\t")
        f.write(str(np.asscalar(RS_BETA)))
        f.write("\n")
        if RefToCotti  ==1:
            f.write("nBeta-RS from Cottrellref.")
            f.write("\t")
            f.write(str(np.asscalar(nBETA_RS)))
            f.write("\n")
        
        f.write("Beta-LN-Plot")
        f.write("\t")
        f.write(str(np.asscalar(LNPLOT_BETA)))
        f.write("\n")
        f.write("Average Beta-Shape-Plot")
        f.write("\t")
        f.write(str(np.asscalar(SHAPE_BETA)))
        f.write("\n")
        f.write("Average IEQ-RS in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(RS_IEQ)))
        f.write("\n")
        f.write("Average IEQ-LN-Plot in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(LNPLOT_IEQ)))
        f.write("\n")
        f.write("Average IEQ-Shape-Plot in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(SHAPE_IEQ)))
        f.write("\n")
        if KZeroCalcer ==1:
            f.write("Average KZero-RS in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(RS_KNULL)))
            f.write("\n")
            f.write("Average KZero-LN-Plot in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(LNPLOT_KNULL)))
            f.write("\n")
            f.write("Average KZero-Shape-Plot in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_KNULL)))
            f.write("\n")
        f.write("\n")
       

        
        
        f.write("Parameters")
        f.write("\n")
        f.write("D in cm^2/s")
        f.write("\t")
        f.write(str(D))
        f.write("\n")
        f.write("A in cm^2")
        f.write("\t")
        f.write(str(A))
        f.write("\n")
        f.write("T in K")
        f.write("\t")
        f.write(str(T))
        f.write("\n")
        f.write("c in mol/cm^3")
        f.write("\t")
        f.write(str(c))
        f.write("\n")
        f.write("Refining Factor")
        f.write("\t")
        f.write(str(REFF))
        f.write("\n")
        f.write("Ohmic Resistance in Ohm")
        f.write("\t")
        f.write(str(Ru))
        f.write("\n")
        f.write("Correction Type")
        f.write("\t")
        if AutoCorr == 1:
            f.write("Auto")
            f.write("\n")
        if BaseCorr == 1:
            f.write("Base")
            f.write("\n")
        if ManuCorr == 1:
            f.write("Manual by")
            f.write("\t")
            f.write(str(ManuCorrCurr))
            f.write("\t")
            f.write("(microams*s)/mV")
            f.write("\n")
        if Smooth_Data == 1:
            f.write("Data got smoothed")
            f.write("\n")
            f.write("Smoothing factor")
            f.write("\t")
            f.write(str(SmoothFact))
            f.write("\n")
            f.write("Smoothing Startpotential in V")
            f.write("\t")
            f.write(str(SmoothFrom))
            f.write("\n")
            f.write("Smoothing Endpotential in V")
            f.write("\t")
            f.write(str(SmoothTo))
            f.write("\n")
        if RefToCotti  ==1:
            f.write("D ref. to Cottrell")
            f.write("\n")
                
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("RS-IRR-Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("Scanrates in mV/s")   
        f.write("\t")
        f.write("Root of Scanrates in (V/s)^0.5")   
        f.write("\t")
        f.write("I-Peak in microampere")   
        f.write("\t")
        f.write("Fit I-Peak in microampere")   
        f.write("\t")
        f.write("Nat. Logarithm of Scanrates")   
        f.write("\t")
        f.write("E-Peak in V")   
        f.write("\t")
        f.write("Fit E-Peak in V")   
        f.write("\t")
        f.write("E-HalfPeak in V")   
        f.write("\t")
        f.write("beta-Shape-Plot")   
        f.write("\t")
        f.write("Ieq-Randles-Plot")   
        f.write("\t")
        f.write("Ieq-LN-Plot")   
        f.write("\t")
        f.write("Ieq-Shape-Plot")   
        f.write("\t")
        if KZeroCalcer ==1:
            f.write("K_Zero-Randles-Plot")   
            f.write("\t")
            f.write("K_Zero-LN-Plot")   
            f.write("\t")
            f.write("K_Zero-Shape-Plot")   
            f.write("\t")
        f.write("\n")
        
        
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(ScanratesArray[i])))
            f.write("\t")
            f.write(str(np.asscalar((WurzelScanratesarray[i]))))
            f.write("\t")
            f.write(str(np.asscalar(I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(FITTED_I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(LnScanratesarray[i])))
            f.write("\t")
            f.write(str(np.asscalar(PeakPotentiale[i])))
            f.write("\t")
            f.write(str(np.asscalar(0.001*FIT_PeakPotentiale[i])))
            f.write("\t")
            f.write(str(np.asscalar(HalbPeakPotenziale[i])))
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_BETAARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(RS_IEQARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(LNPLOT_IEQARRAY[i])))
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_IEQARRAY[i])))
            f.write("\t")
            if KZeroCalcer ==1:
                f.write(str(np.asscalar(RS_KNULLARRAY[i])))
                f.write("\t")
                f.write(str(np.asscalar(LNPLOT_KNULLARRAY[i])))
                f.write("\t")
                f.write(str(np.asscalar(SHAPE_KNULLARRAY[i])))
                f.write("\t")
            
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Scanrate in mV/s")
        f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            f.write(str(np.asscalar(ScanratesArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(KORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(KORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")
#__________________________________________________________________________________________________-        
        for i in range(3):
            f.write("\n")
            
        f.write("Uncorrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Scanrate in mV/s")
        f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            f.write(str(np.asscalar(ScanratesArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(ShortLength)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(UNKORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(UNKORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")

    root.destroy()



 


# In[27]:


def RSIRRCONCASTXTsaver():
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")

        f.write("Beta-RS")
        f.write("\t")
        f.write(str(np.asscalar(RS_BETACONC)))
        f.write("\n")
        f.write("Average Beta-Shape-Plot")
        f.write("\t")
        f.write(str(np.asscalar(SHAPE_BETACONC)))
        f.write("\n")
        f.write("Average IEQ-RS in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(RS_IEQCONC)))
        f.write("\n")
        f.write("Average IEQ-Shape-Plot in A/cm^2")
        f.write("\t")
        f.write(str(np.asscalar(SHAPE_IEQCONC)))
        f.write("\n")
        if KZeroCalcer ==1:
            f.write("Average KZero-RS in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(RS_KNULLCONC)))
            f.write("\n")
            f.write("Average KZero-Shape-Plot in cm/s")
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_KNULLCONC)))
            f.write("\n")
        f.write("\n")
       

        
        
        f.write("Parameters")
        f.write("\n")
        f.write("D in cm^2/s")
        f.write("\t")
        f.write(str(D))
        f.write("\n")
        f.write("A in cm^2")
        f.write("\t")
        f.write(str(A))
        f.write("\n")
        f.write("T in K")
        f.write("\t")
        f.write(str(T))
        f.write("\n")
        f.write("Highest c in mol/cm^3")
        f.write("\t")
        f.write(str(c))
        f.write("\n")
        f.write("Refining Factor")
        f.write("\t")
        f.write(str(REFF))
        f.write("\n")
        f.write("Ohmic Resistance in Ohm")
        f.write("\t")
        f.write(str(Ru))
        f.write("\n")
        f.write("Correction Type")
        f.write("\t")
        if AutoCorr == 1:
            f.write("Auto")
            f.write("\n")
        if BaseCorr == 1:
            f.write("Base")
            f.write("\n")
        if ManuCorr == 1:
            f.write("Manual by")
            f.write("\t")
            f.write(str(ManuCorrCurr))
            f.write("\t")
            f.write("(microams*s)/mV")
            f.write("\n")
        if Smooth_Data == 1:
            f.write("Data got smoothed")
            f.write("\n")
            f.write("Smoothing factor")
            f.write("\t")
            f.write(str(SmoothFact))
            f.write("\n")
            f.write("Smoothing Startpotential in V")
            f.write("\t")
            f.write(str(SmoothFrom))
            f.write("\n")
            f.write("Smoothing Endpotential in V")
            f.write("\t")
            f.write(str(SmoothTo))
            f.write("\n")               
        
        
        for i in range(3):
            f.write("\n")
        
        f.write("RS-IRR-Data")
        f.write("\n")
        f.write("\n")
        f.write("OCPs in V")   
        f.write("\t")
        f.write("Concentrations in mol/cm^3")   
        f.write("\t")
        f.write("I-Peak in microampere")   
        f.write("\t")
        f.write("Fit I-Peak in microampere")   
        f.write("\t")
        f.write("E-Peak in V")    
        f.write("\t")
        f.write("E-HalfPeak in V")   
        f.write("\t")
        f.write("beta-Shape-Plot")   
        f.write("\t")
        f.write("Ieq-Randles-Plot")   
        f.write("\t")
        f.write("Ieq-Shape-Plot")   
        f.write("\t")
        if KZeroCalcer ==1:
            f.write("K_Zero-Randles-Plot")   
            f.write("\t")
            f.write("K_Zero-Shape-Plot")   
            f.write("\t")
        f.write("\n")
        
        
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(OCPArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(CDecayArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(FITTED_I_Peak[i])))
            f.write("\t")
            f.write(str(np.asscalar(PeakPotentiale[i])))
            f.write("\t")
            f.write(str(np.asscalar(HalbPeakPotenziale[i])))
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_BETAARRAYCONC[i])))
            f.write("\t")
            f.write(str(np.asscalar(RS_IEQARRAYCONC[i])))
            f.write("\t")
            f.write(str(np.asscalar(SHAPE_IEQARRAYCONC[i])))
            f.write("\t")
            if KZeroCalcer ==1:
                f.write(str(np.asscalar(RS_KNULLARRAYCONC[i])))
                f.write("\t")
                f.write(str(np.asscalar(SHAPE_KNULLARRAYCONC[i])))
                f.write("\t")
            
            f.write("\n")
            

        
        for i in range(3):
            f.write("\n")
#_________________________________________________________________________________            
        f.write("Korrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Concentrations as fract. of highest")
        f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            f.write(str(np.asscalar(ConcentrationsArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(Length)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(KORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(KORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")
#__________________________________________________________________________________________________-        
        for i in range(3):
            f.write("\n")
            
        f.write("Uncorrected LSVs")
        f.write("\n")
        f.write("\n")
        f.write("Scanrate in mV/s")
        f.write("\t")
    
        for jj in range(int(NumMeas)):
        
            f.write(str(np.asscalar(ConcentrationsArray[jj])))
            f.write("\t")
            f.write("\t")
        
        f.write("\n")
        f.write("\t")
    
    
    
        for jj in range(int(NumMeas)):
        
            f.write("E in V")
            f.write("\t")
            f.write("I in micAmp")
            f.write("\t")
        
        f.write("\n")
        f.write("\n")
        f.write("\n")
        
        
    
        for i in range(int(ShortLength)):
        
            f.write(str([i]))
            f.write("\t")
            for j in range(int(NumMeas)):
                b = str(np.asscalar(UNKORRPOTENZIALSUPPERARRAY[j:j+1:,i:i+1:]))
                d = str(np.asscalar(UNKORRSTROMSUPPERARRAY[j:j+1:,i:i+1:]))
                f.write(b)
                f.write("\t")
                f.write(d)
                f.write("\t")
            f.write("\n")

    root.destroy()



 


# In[28]:


def CottrellWindow():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Cottrell-Analysis")                         
    Fenster.geometry("150x200")
    

    var3 = IntVar()
    Checkbutton(Fenster, text="Baseline Correction", variable=var3).grid(row=1,column=0, sticky=W)
    
    c_Label = Label(Fenster,text="c in mol/cm^3")
    c_Label.grid(row=2, column=0)
    
    Concents = []

    for i in range(NumMeas):
        
        en = Entry(Fenster)
        en.grid(row=i+3, column=0)
        Concents.append(en)    
        
    
    def AcceptParams():

        global BaseCorr
        
        Concs= []
        for entry in Concents:
            Concs.append(float(entry.get()))
        global ConcentrationsArray
        ConcentrationsArray = np.array(Concs) 
        
        BaseCorr = var3.get()
        
    def Next():
        
        def quit():
            Fenster.destroy()
        quit() 
        
        if BaseCorr == 0:
            CottrellWindowLevel2()
        if BaseCorr == 1:
            BaseGetter()
            CottrellWindowLevel2()
        
        
    
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=30,column=0)
    button=Button(Fenster,text="Next",command=Next).grid(row=31,column=0)


# In[29]:


def CottrellWindowLevel2():
    Fenster = Toplevel()                                                         
    Fenster.title("Cottrell-Analysis")                         
    Fenster.geometry("1200x600")
    
    
    f = Figure(figsize=(12, 6), dpi=100)

    global ShortLength
    ShortLength = int(len(np.squeeze(np.transpose(Potenzial[0::,0:1:]))))
    
    global ZEITSUPERARRAY
    global INVWURZZEITSUPPERARRAY
    global UNKORRSTROMSUPPERARRAY
    global STROMSUPPERARRAY
    
    ZEITSUPERARRAY             = np.empty([NumMeas,ShortLength])
    INVWURZZEITSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    UNKORRSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    STROMSUPPERARRAY           = np.empty([NumMeas,ShortLength])
    
    
    for i in range(int(NumMeas)):
        StromarraysROH             = np.squeeze(np.transpose(Strom[0::,i:i+1:]))   
        ZeitArraysROH              = np.squeeze(np.transpose(Potenzial[0::,i:i+1:]))
        InvWurzZeit                = 1/((ZeitArraysROH**0.5))
        
        if DesiReactxx == 0:
            Stromarrays            = -StromarraysROH
            
        if DesiReactxx == 1:
            Stromarrays            = StromarraysROH
   
            
        if BaseCorr == 1:
            
            BaseStromArraysROH     = np.squeeze(np.transpose(BaseStrom[0::,i:i+1:]))
            if DesiReactxx == 0:
                BaseStromArrays    = -BaseStromArraysROH[::]
            if DesiReactxx == 1:
                BaseStromArrays    = BaseStromArraysROH[::]
            Stromarrays            = Stromarrays - BaseStromArrays
                                              

        for j in range(ShortLength):
                
            ZEITSUPERARRAY[i,j]          = ZeitArraysROH[j] 
            INVWURZZEITSUPPERARRAY[i,j]  = InvWurzZeit[j]
            UNKORRSTROMSUPPERARRAY[i,j]  = StromarraysROH[j]
            STROMSUPPERARRAY[i,j]        = Stromarrays[j]

        
        
    #Chronoplot
    #______________________________________________________
    Chronoplot = f.add_subplot(121)
    for i in range(NumMeas):
        Chronoplot.plot(ZEITSUPERARRAY[i],UNKORRSTROMSUPPERARRAY[i], color='k',linestyle='-',marker='')
        Chronoplot.plot(ZEITSUPERARRAY[i],STROMSUPPERARRAY[i],color='r',linestyle='-',marker='')
    Chronoplot.set_xlabel('t in [s]', fontsize=12)
    Chronoplot.set_ylabel('I [$\mu$A]', fontsize=12)
    for axis in ['top','bottom','left','right']:
        Chronoplot.spines[axis].set_linewidth(2)
    Chronoplot.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
    #Cottiplot für Fit
    #_______________________________________________________
        
    Cottiplot = f.add_subplot(122)
    for i in range(NumMeas):
        Cottiplot.plot(INVWURZZEITSUPPERARRAY[i],STROMSUPPERARRAY[i],color='r',linestyle='',marker='.')
    Cottiplot.set_xlabel(r'$\sqrt{t^{-1}}$' '[s $^{-0.5}$]''t in [s]', fontsize=12)
    Cottiplot.set_ylabel('I [$\mu$A]', fontsize=12)
    for axis in ['top','bottom','left','right']:
        Cottiplot.spines[axis].set_linewidth(2)
    Cottiplot.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
        
        
    f.tight_layout()
    canvas = FigureCanvasTkAgg(f, master=Fenster)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, Fenster)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    

    CottrellWindowLevel3()
    


# In[30]:


def CottrellWindowLevel3():
    
    Fenster = Toplevel()                                                         
    Fenster.title("Cottrell-Analysis")                         
    Fenster.geometry("150x200")
    
    FitFrom_Label = Label(Fenster, text="Fit From")
    FitFrom_Label.grid(row=1, column=0)
    FitFrom_Eingabe = Entry(Fenster)                                               
    FitFrom_Eingabe.grid(row=2, column=0)

    FitTo_Label = Label(Fenster,text="Fit To")
    FitTo_Label.grid(row=3, column=0)
    FitTo_Eingabe = Entry(Fenster)
    FitTo_Eingabe.grid(row=4, column=0)
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Save as txt", variable=var1).grid(row=5,column=0, sticky=W)
    
    def Next():
        
        global FitFrom
        global FitTo
        global AsTxtSaver
        
        FitFrom    = float(FitFrom_Eingabe.get())
        FitTo      = float(FitTo_Eingabe.get())
        AsTxtSaver = int(var1.get())
        
        def quit():
            Fenster.destroy()
        quit() 
        
        CottrellWindowLevel4()
        
        
        
    button=Button(Fenster,text="Next",command=Next).grid(row=6,column=0)
    
    
    


# In[31]:


def CottrellWindowLevel4():  
    
    Fenster = Toplevel()                                                         
    Fenster.title("Cottrell-Plot")                         
    Fenster.geometry("700x700")
    
    #DIESE ZWEI FUNKTIONEN MÜSSEN HIER DRIN DEFINIERT WERDEN!!!
        
    def GeradenFit(xWert,ySchnitt,Steig):
        return  ySchnitt + Steig*xWert

    def EchtesPotFinden(WoSieSuchenSoll,WasSieSuchenSoll):                         
        idxStart = (np.abs(WoSieSuchenSoll-WasSieSuchenSoll)).argmin()               
        return WoSieSuchenSoll[idxStart] 
    
    

    global FITSTROMSUPPERARRAY
    FITSTROMSUPPERARRAY     = np.empty([NumMeas,ShortLength])
    global CottiSlopesArray
    CottiSlopesArray        = np.empty(int(NumMeas))
    global CottiCutsArray
    CottiCutsArray          = np.empty(int(NumMeas))
    
    
    f = Figure(figsize=(6, 6), dpi=100)

          
    
    for i in range(int(NumMeas)):
    
        FitFromPotenzial       = EchtesPotFinden(INVWURZZEITSUPPERARRAY[i,::],FitFrom)
        FitToPotenzial         = EchtesPotFinden(INVWURZZEITSUPPERARRAY[i,::],FitTo)
        
        #print FitFromPotenzial
        #print FitToPotenzial
        
        IdxFitFrom = np.asscalar(np.where(INVWURZZEITSUPPERARRAY[i,::] == FitFromPotenzial) [0])
        IdxFitTo   = np.asscalar(np.where(INVWURZZEITSUPPERARRAY[i,::] == FitToPotenzial) [0])
        
        if FitFromPotenzial<FitToPotenzial:
            CottiFitting, pcov  = curve_fit(GeradenFit, INVWURZZEITSUPPERARRAY[i,IdxFitTo:IdxFitFrom], STROMSUPPERARRAY[i,IdxFitTo:IdxFitFrom])
        
        if FitFromPotenzial>FitToPotenzial:
            CottiFitting, pcov  = curve_fit(GeradenFit, INVWURZZEITSUPPERARRAY[i,IdxFitFrom:IdxFitTo], STROMSUPPERARRAY[i,IdxFitFrom:IdxFitTo])
        
        
        CottiSlopesArray[i] = CottiFitting[1]
        CottiCutsArray[i]   = CottiFitting[0]
        
    #Fitstromarrays berechnen
    #________________________________________________________________
    for i in range(len(CottiSlopesArray)):
 
        FITSTROMSUPPERARRAY[i] = INVWURZZEITSUPPERARRAY[i,::]*CottiSlopesArray[i] + CottiCutsArray[i]       
            

    #PLOTTEN
    #_________________________________________________________________
    #Cottiplot für Fit
    #_________________________________________________________________
        
    Cottiplot = f.add_subplot(111)
    for i in range(NumMeas):
        Cottiplot.plot(INVWURZZEITSUPPERARRAY[i],STROMSUPPERARRAY[i],color='r',linestyle='',marker='.')
        Cottiplot.plot(INVWURZZEITSUPPERARRAY[i],FITSTROMSUPPERARRAY[i],color='r',linestyle='-',marker='')
    
        Cottiplot.annotate('m=%5.3f'    % CottiSlopesArray[i], xy=(0.15, 0.93 - 0.05*i), xycoords='axes fraction')
    
    
    Cottiplot.set_xlabel(r'$\sqrt{t^{-1}}$' '[s $^{-0.5}$]''t in [s]', fontsize=12)
    Cottiplot.set_ylabel('I [$\mu$A]', fontsize=12)
    for axis in ['top','bottom','left','right']:
        Cottiplot.spines[axis].set_linewidth(2)
    Cottiplot.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
        
        
        
        
    f.tight_layout()
    canvas = FigureCanvasTkAgg(f, master=Fenster)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, Fenster)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)    
    
    
    if AsTxtSaver ==1:
        CottiAsTxtSaver()
       


# In[32]:


def CottiAsTxtSaver():
    
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as f:
        f.write("Results")
        f.write("\n")
        f.write("\n")
        f.write("c in mol/cm^3")
        f.write("\t")
        f.write("Cottrellslope in (µA cm^3)/(s^0.5 mol)")
        f.write("\n")
        f.write("\n")
        
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(ConcentrationsArray[i])))
            f.write("\t")
            f.write(str(np.asscalar(CottiSlopesArray[i])))
            f.write("\t")
            f.write("\n")
        
        f.write("\n")
        
        f.write("Correction type")
        f.write("\t")
        if BaseCorr == 0:
            f.write("none")
        if BaseCorr == 1:
            f.write("Base")
        
        f.write("\n")
        f.write("Fit from")
        f.write("\t")
        f.write(str((FitFrom)))
        f.write("\t")
        f.write("\n")
        f.write("Fit to")
        f.write("\t")
        f.write(str((FitTo)))
        f.write("\t")
        
        f.write("\n")
        f.write("\n")
        
        f.write("Data")
        f.write("\n")
        
        for i in range(int(NumMeas)):
            f.write(str(np.asscalar(ConcentrationsArray[i])))
            f.write("\t")
            f.write("\t")
            f.write("\t")
            f.write("\t")
        f.write("\n")
        
        for i in range(int(NumMeas)):
            f.write("t in s")
            f.write("\t")
            f.write(" 1/t^0.5 in 1/s^0.5")
            f.write("\t")
            f.write("I in microamps")
            f.write("\t")
            f.write("Fit-I in microamps")
            f.write("\t")
        f.write("\n")
        
        for j in range(int(ShortLength)):
            for i in range(int(NumMeas)):
                f.write(str(np.asscalar(ZEITSUPERARRAY[i,j])))
                f.write("\t")
                f.write(str(np.asscalar(INVWURZZEITSUPPERARRAY[i,j])))
                f.write("\t")
                f.write(str(np.asscalar(STROMSUPPERARRAY[i,j])))
                f.write("\t")
                f.write(str(np.asscalar(FITSTROMSUPPERARRAY[i,j])))
                f.write("\t")
            
            f.write("\n")
            

    root.destroy()


# In[33]:


def Model_Zeroer():
    global model_a
    global model_b
    global model_c
    global model_d
    global model_e
    global model_f
    global model_g
    global Statistical
    Statistical     = 0
    model_a         = 0
    model_b         = 0
    model_c         = 0
    model_d         = 0
    model_e         = 0
    model_f         = 0
    model_g         = 0
    global xxx
    global yyy
    global PsiArray
    global XiArray
    global CV_Interpolation_Hin
    global CV_Interpolation_Back


def Semi_Inf_Planar():
    Model_Zeroer()
    global model_a
    model_a = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
       
def Finit_Planar():
    Model_Zeroer()
    global model_b
    model_b = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Semi_Inf_Zyl_Ext():
    Model_Zeroer()
    global model_c
    model_c = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Finit_Zyl_Ext():
    Model_Zeroer()
    global model_d
    model_d = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Finit_Zyl_Int():
    Model_Zeroer()
    global model_e
    model_e = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Semi_Inf_Sphere_Ext():
    Model_Zeroer()
    global model_f
    model_f = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Finit_Sphere_Int():
    Model_Zeroer()
    global model_g
    model_g = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
    
    
#-------------------------------------------------------------------------------------
#Fitters
#-------------------------------------------------------------------------------------
    
    
def Semi_Inf_Planar_FITTER():
    Model_Zeroer()
    global model_a
    model_a = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def Finit_Planar_FITTER():
    Model_Zeroer()
    global model_b
    model_b = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def Semi_Inf_Zyl_Ext_FITTER():
    Model_Zeroer()
    global model_c
    model_c = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def Finit_Zyl_Ext_FITTER():
    Model_Zeroer()
    global model_d
    model_d = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def Finit_Zyl_Int_FITTER():
    Model_Zeroer()
    global model_e
    model_e = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def Semi_Inf_Sphere_Ext_FITTER():
    Model_Zeroer()
    global model_f
    model_f = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
def Finit_Sphere_Int_FITTER():
    Model_Zeroer()
    global model_g
    model_g = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
        
        
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#Statistical Simulators
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------

def Statistical_Finit_Planar():
    Model_Zeroer()
    global model_b
    global Statistical
    Statistical  = 1
    model_b = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Statistical_Finit_Zyl_Ext():
    Model_Zeroer()
    global model_d
    model_d = 1
    global Statistical
    Statistical  = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
def Statistical_Finit_Zyl_Int():
    Model_Zeroer()
    global model_e
    global Statistical
    Statistical  = 1
    model_e = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
    
def Statistical_Finit_Sphere_Int():
    Model_Zeroer()
    global model_g
    global Statistical
    Statistical  = 1
    model_g = 1
    global FIT 
    FIT = 0
    Eingabe_CV_Simulator()
    
    
    
#-------------------------------------------------------------------------------------
#Statistical Fitters
#-------------------------------------------------------------------------------------
    
    
def Statistical_Finit_Planar_FITTER():
    Model_Zeroer()
    global model_b
    global Statistical
    Statistical  = 1
    model_b = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    
    
def Statistical_Finit_Zyl_Ext_FITTER():
    Model_Zeroer()
    global model_d
    global Statistical
    Statistical  = 1
    model_d = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
def Statistical_Finit_Zyl_Int_FITTER():
    Model_Zeroer()
    global model_e
    global Statistical
    Statistical  = 1
    model_e = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
        
    
def Statistical_Finit_Sphere_Int_FITTER():
    Model_Zeroer()
    global model_g
    global Statistical
    Statistical  = 1
    model_g = 1
    global FIT 
    FIT = 1
    if File_Was_Loaded == 1:
        Eingabe_CV_Simulator()
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()


# In[34]:


def No_Loaded_File_Warner():
    Fenster = Toplevel()                                                         
    Fenster.title("No loaded file found")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Caution! You loaded no file that\nshould be fittet. Go to\nopen file first and select\na file that should be fitted. \n Nevertheless, CV might be calculated")


# In[35]:


def Parameters_do_not_fit_warner():
    Fenster = Toplevel()                                                         
    Fenster.title("Entered Parameters are out of limits!")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Caution! The entered Parameters\ncannot be used for the\ncalculation of a CV\nIt should be always\nD>10^-6 cm^2/s,  0.00001 cm < distance < 10 cm \nand r>0.1 micrometer!")


# In[36]:


def Eingabe_CV_Simulator():
    Fenster = Toplevel()  
    if FIT == 0:
        Fenster.geometry("700x630")
    if FIT == 1:
        Fenster.geometry("700x630")
        
    
    Fenster.title("Set Parameters")        
    
    
    n_Label = Label(Fenster, text="n*")
    n_Label.grid(row=0, column=0)
    n_Eingabe = Entry(Fenster)                                               
    n_Eingabe.grid(row=0, column=1)
    
    alpha_Label = Label(Fenster, text="alpha_forward*")
    alpha_Label.grid(row=1, column=0)
    alpha_Eingabe = Entry(Fenster)                                               
    alpha_Eingabe.grid(row=1, column=1)
    
    kzero_Label = Label(Fenster,text="k_zero [cm/s]*")
    kzero_Label.grid(row=2, column=0)
    kzero_Eingabe = Entry(Fenster)
    kzero_Eingabe.grid(row=2, column=1)
    
    var15 = IntVar()
    Checkbutton(Fenster, text="Fin. het. k [cm/s]", variable=var15).grid(row=3, column=0, sticky=W)
    kmax_Eingabe = Entry(Fenster)
    kmax_Eingabe.grid(row=3, column=1)
    
    
    
    
    if FIT ==1:
        Ezero_Label = Label(Fenster,text="E_zero vs Ref [V]*")
        Ezero_Label.grid(row=4, column=0)
        Ezero_Eingabe = Entry(Fenster)
        Ezero_Eingabe.grid(row=4, column=1)
    
    T_Label = Label(Fenster,text="T [°C]*")
    T_Label.grid(row=5, column=0)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.grid(row=5, column=1)
    
    Scanrate_Label = Label(Fenster,text="Scanrate [mV/s]*")
    Scanrate_Label.grid(row=6, column=0)
    Scanrate_Eingabe = Entry(Fenster)
    Scanrate_Eingabe.grid(row=6, column=1)
    
    if Statistical  == 1:
        
        DistFunc_Nodes_Label = Label(Fenster,text="Sub int. of Distribution*")
        DistFunc_Nodes_Label.grid(row=10, column=0)
        DistFunc_Nodes_Eingabe = Entry(Fenster)
        DistFunc_Nodes_Eingabe.grid(row=10, column=1)
        
        if model_b == 1:
            Number_N_Label = Label(Fenster,text="sheets per mm*")
        if model_d == 1 or model_e == 1:
            Number_N_Label = Label(Fenster,text="Cylinder bottoms per mm^2*")
        if model_g == 1:
            Number_N_Label = Label(Fenster,text="Spheres per mm^3*")
        Number_N_Label.grid(row=11, column=0)
        Number_N_Eingabe = Entry(Fenster)
        Number_N_Eingabe.grid(row=11, column=1)
        
        
        var16 = IntVar()
        Checkbutton(Fenster, text="Modify Distribution function", variable=var16).grid(row=12, column=0, sticky=W)
        
        DFP_a_Label = Label(Fenster,text="Dist.func. param. a")
        DFP_a_Label.grid(row=13, column=0)
        DFP_a_Eingabe = Entry(Fenster)
        DFP_a_Eingabe.grid(row=13, column=1)
            
        DFP_b_Label = Label(Fenster,text="Dist.func. param. b")
        DFP_b_Label.grid(row=14, column=0)
        DFP_b_Eingabe = Entry(Fenster)
        DFP_b_Eingabe.grid(row=14, column=1)
        
        DFP_c_Label = Label(Fenster,text="Dist.func. param. c")
        DFP_c_Label.grid(row=15, column=0)
        DFP_c_Eingabe = Entry(Fenster)
        DFP_c_Eingabe.grid(row=15, column=1)
        
        
        if model_d == 1:
            r_Label = Label(Fenster,text="Fiber radius 10^-6[m]*")
            r_Label.grid(row=16, column=0)
            r_Eingabe = Entry(Fenster)
            r_Eingabe.grid(row=16, column=1)
        
        if model_b == 1:
            sheet_thickness_Label = Label(Fenster,text="Sheet thickness 10^-6[m]*")
            sheet_thickness_Label.grid(row=16, column=0)
            sheet_thickness_Eingabe = Entry(Fenster)
            sheet_thickness_Eingabe.grid(row=16, column=1)
            
             
            
    if Statistical  == 0:
        if model_b == 1:
            d_Label = Label(Fenster,text="distance 10^-6[m]*")
            d_Label.grid(row=14, column=0)
            d_Eingabe = Entry(Fenster)
            d_Eingabe.grid(row=14, column=1)
    
        if model_c == 1 or model_d  == 1 or model_e == 1 or model_f == 1 or model_g == 1:
            r_Label = Label(Fenster,text="radius 10^-6[m]*")
            r_Label.grid(row=14, column=0)
            r_Eingabe = Entry(Fenster)
            r_Eingabe.grid(row=14, column=1)
    
        if model_d  == 1:
            d_Label = Label(Fenster,text="distance 10^-6[m]*")
            d_Label.grid(row=15, column=0)
            d_Eingabe = Entry(Fenster)
            d_Eingabe.grid(row=15, column=1)
        
    

    if FIT ==0:
        End_Label = Label(Fenster,text="Abs(E_fin) vs. Ezero [V]*")
        End_Label.grid(row=17, column=0)
        End_Eingabe = Entry(Fenster)
        End_Eingabe.grid(row=17, column=1)
        
    A_Label = Label(Fenster,text="A in cm^2*")
    A_Label.grid(row=18, column=0)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.grid(row=18, column=1)
        
    c_Label = Label(Fenster,text="c in mol/L*")
    c_Label.grid(row=19, column=0)
    c_Eingabe = Entry(Fenster)
    c_Eingabe.grid(row=19, column=1)    
           
    
    if FIT ==1:
        Ru_Label = Label(Fenster,text="Ru in Ohm*")
        Ru_Label.grid(row=20, column=0)
        Ru_Eingabe = Entry(Fenster)
        Ru_Eingabe.grid(row=20, column=1)
        
    
    D_Label = Label(Fenster,text="D 10^-6[cm^2/s]*")
    D_Label.grid(row=21, column=0)
    D_Eingabe = Entry(Fenster)
    D_Eingabe.grid(row=21, column=1)
    
    var17 = IntVar()
    Checkbutton(Fenster, text="Mod. Xi-Resol.", variable=var17).grid(row=22, column=0, sticky=W)
    
    ModDeltXi_Label = Label(Fenster,text="delta Xi")
    ModDeltXi_Label.grid(row=23, column=0)
    ModDeltXi_Eingabe = Entry(Fenster)
    ModDeltXi_Eingabe.grid(row=23, column=1)
    
    
    
    
    #-----------------------------------------------------------------------------------------------
    #Neue Spalte im Fenster
    #-----------------------------------------------------------------------------------------------
    
    
    if FIT ==0:
        var4 = IntVar()
        Checkbutton(Fenster, text="Anodic sweep", variable=var4).grid(row=0, column=3, sticky=W)
    
        var5 = IntVar()
        Checkbutton(Fenster, text="Cathodic sweep", variable=var5).grid(row=1, column=3, sticky=W)
     
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Unequal D", variable=var1).grid(row=2, column=3, sticky=W)
    
    Df_Label = Label(Fenster,text="D_f 10^-6[cm^2/s]")
    Df_Label.grid(row=3, column=4)
    Df_Eingabe = Entry(Fenster)
    Df_Eingabe.grid(row=3, column=3)
    
    Db_Label = Label(Fenster,text="D_b 10^-6[cm^2/s]")
    Db_Label.grid(row=4, column=4)
    Db_Eingabe = Entry(Fenster)
    Db_Eingabe.grid(row=4, column=3)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Preceding ch. eq.", variable=var2).grid(row=5, column=3, sticky=W)
    
    kp_Label = Label(Fenster,text="k_p [1/s]")
    kp_Label.grid(row=6, column=4)
    kp_Eingabe = Entry(Fenster)
    kp_Eingabe.grid(row=6, column=3)
    
    kmp_Label = Label(Fenster,text="k_-p [1/s]")
    kmp_Label.grid(row=7, column=4)
    kmp_Eingabe = Entry(Fenster)
    kmp_Eingabe.grid(row=7, column=3)
    
    var3 = IntVar()
    Checkbutton(Fenster, text="Following ch. eq.", variable=var3).grid(row=8, column=3, sticky=W)
    
    kf_Label = Label(Fenster,text="k_f [1/s]")
    kf_Label.grid(row=9, column=4)
    kf_Eingabe = Entry(Fenster)
    kf_Eingabe.grid(row=9, column=3)
    
    kmf_Label = Label(Fenster,text="k_-f [1/s]")
    kmf_Label.grid(row=10, column=4)
    kmf_Eingabe = Entry(Fenster)
    kmf_Eingabe.grid(row=10, column=3)
    
    
    if FIT ==1:
        var7 = IntVar()
        Checkbutton(Fenster, text="Lin. Basecorr", variable=var7).grid(row=11, column=3, sticky=W)
        
        LinStart_Label = Label(Fenster,text="Start vs. Ref. [V]")
        LinStart_Label.grid(row=12, column=4)
        LinStart_Eingabe = Entry(Fenster)
        LinStart_Eingabe.grid(row=12, column=3)
        
        LinEnd_Label = Label(Fenster,text="End vs. Ref. [V]")
        LinEnd_Label.grid(row=13, column=4)
        LinEnd_Eingabe = Entry(Fenster)
        LinEnd_Eingabe.grid(row=13, column=3)
        
    
    if FIT == 0:
        var8 = IntVar()
        Checkbutton(Fenster, text="Show Xi vs. Chi", variable=var8).grid(row=14, column=3, sticky=W)

        var9 = IntVar()
        Checkbutton(Fenster, text="Save Chi vs. Xi", variable=var9).grid(row=15, column=3, sticky=W)
        
        if Statistical  == 0:
            var10 = IntVar()
            Checkbutton(Fenster, text="Save convolution functions", variable=var10).grid(row=16, column=3, sticky=W)

    if Statistical  == 1:
        if FIT == 0:
            var11 = IntVar()
            Checkbutton(Fenster, text="Show Individuals", variable=var11).grid(row=17, column=3, sticky=W)
        
        var12 = IntVar()
        Checkbutton(Fenster, text="Save Individuals", variable=var12).grid(row=18, column=3, sticky=W)
        
        var13 = IntVar()
        Checkbutton(Fenster, text="Show Distribution function", variable=var13).grid(row=19, column=3, sticky=W)
        
        var14 = IntVar()
        Checkbutton(Fenster, text="Save Distribution function", variable=var14).grid(row=20, column=3, sticky=W)
        
        
    
    var6 = IntVar()
    Checkbutton(Fenster, text="Save as txt", variable=var6).grid(row=21, column=3, sticky=W)
    

    
    #-------------------------------------------------------
    #Akzeptierfunktion schreiben
    #-------------------------------------------------------
    
    def AcceptParams():
        global n
        global alpha
        global kzero
        global kp
        global kmp
        global kf
        global kmf
        global Ezero
        global T
        global distance
        global r
        global Scanrate
        global Precedeq
        global Succeedeq
        global Unequal_D
        global Cathsweeper
        global AstxtSaver
        global A
        global Area
        global c
        global Ru
        global LinBaseCorr
        global D_f
        global D_b
        global Show_Individuals
        global Save_Individuals
        global Show_Dist_Func
        global Save_Dist_Func
        global DistFunc_Nodes
        global DistFunc_Modifier
        global average_d_or_r
        global Number_N
        global DFP_a
        global DFP_b
        global DFP_c
        global sheet_thickness
        global FinKin
        global kfin
        global DeltaXiModder
        
#========================================================================
    #Read entries
#========================================================================
        
        n             = (float(n_Eingabe.get()))
        alpha         = (float(alpha_Eingabe.get()))
        kzero         = (float(kzero_Eingabe.get()))
        T             = (float(T_Eingabe.get())) + 273.15
        Scanrate      = (float(Scanrate_Eingabe.get()))*0.001
        FinKin        = var15.get()
        kfin          = 0
        c             = (float(c_Eingabe.get()))*0.001
        A             = (float(A_Eingabe.get()))
        AstxtSaver    = var6.get()
        DeltaXiModder = var17.get()
        
        if DeltaXiModder == 1:
            global ModDeltXi
            ModDeltXi = (float(ModDeltXi_Eingabe.get()))
            
        
        if FinKin == 1:
            kfin = kzero/(float(kmax_Eingabe.get()))
        
        
        
        if Statistical  == 0:
            if model_b == 1:
                distance     = (float(d_Eingabe.get()))*0.0001
            if model_c == 1  or model_f == 1 or model_e == 1  or model_g == 1:
                r            = (float(r_Eingabe.get()))*0.0001
            if model_d  == 1:
                distance     = (float(d_Eingabe.get()))*0.0001
                r            = (float(r_Eingabe.get()))*0.0001
            
        Unequal_D   = var1.get()
        if Unequal_D ==1:
            D_f     = (float(Df_Eingabe.get()))*0.000001
            D_b     = (float(Db_Eingabe.get()))*0.000001
        if Unequal_D ==0:
            global D
            D_f       = (float(D_Eingabe.get()))*0.000001
            D_b       = (float(D_Eingabe.get()))*0.000001
        
        Precedeq = var2.get()
        if Precedeq ==1:
            kp    = (float(kp_Eingabe.get()))
            kmp   = (float(kmp_Eingabe.get()))
           
        
        Succeedeq = var3.get()
        if Succeedeq ==1:
            kf   = (float(kf_Eingabe.get()))
            kmf  = (float(kmf_Eingabe.get()))
        
        
        if FIT ==0:
            Cathsweeper = var5.get()
            global EEE
            EEE         = (float(End_Eingabe.get()))
            global Just_Current_Template
            Just_Current_Template = var8.get()
            global Chi_vs_Xi_Saver
            Chi_vs_Xi_Saver = var9.get()
            
            if Statistical  == 0:
                global Convolution_Saver
                Convolution_Saver = var10.get()
        
        
        if FIT == 1:  
            Ru     = (float(Ru_Eingabe.get()))
            Ezero  = (float(Ezero_Eingabe.get()))
            LinBaseCorr = var7.get()
            if LinBaseCorr == 1:
                global LineStart
                global LineEnd
                LineStart = (float(LinStart_Eingabe.get()))
                LineEnd   = (float(LinEnd_Eingabe.get()))
        
       
         
        
#========================================================================
#Define distribution function
#========================================================================
        if Statistical  == 1:
            
            Number_N = (float(Number_N_Eingabe.get()))
            #there are Number_N sheets per mm, cylinders per mm^2 or spheres per mm^3  
            #--> the following calculates r_average or d_average in µm
            if model_b ==1:
                average_d_or_r = (1000*1.0/(Number_N))    #distance between two sheet centers
            if model_d == 1 or model_e == 1:
                average_d_or_r = 2*(1000000/(np.pi*Number_N))**0.5
            if model_g == 1:
                average_d_or_r = 2 * (1000000000*3.0/(4.0*np.pi*Number_N))**0.333333333
                
            
            
            DistFunc_Modifier = var16.get()
            
            if DistFunc_Modifier == 0:
                if model_b ==1:
                    #1D dist-Function
                    DFP_a =  2.0
                    DFP_b =  2.0
                    DFP_c =  1.0
    
                if model_d ==1 or model_e == 1:
                    #2D dist-Function
                    DFP_a =  3.3095
                    DFP_b =  3.0328
                    DFP_c =  1.0787

                if model_g ==1:
                    #3D dist-Function
                    DFP_a =  4.8065
                    DFP_b =  4.06342
                    DFP_c =  1.16391
            
            if DistFunc_Modifier == 1:
                DFP_a =  (float(DFP_a_Eingabe.get()))
                DFP_b =  (float(DFP_b_Eingabe.get()))
                DFP_c =  (float(DFP_c_Eingabe.get()))
                
            if FIT == 0:
                Show_Individuals  = var11.get()
            Save_Individuals  = var12.get()
            Show_Dist_Func    = var13.get()
            Save_Dist_Func    = var14.get()
            DistFunc_Nodes    = 1+(int(DistFunc_Nodes_Eingabe.get()))
            
            
            if model_d  == 1:
                r               = (float(r_Eingabe.get()))*0.0001                 #Cylinder radius in cm (given entry in µm)
            if model_b == 1:
                sheet_thickness = (float(sheet_thickness_Eingabe.get()))*0.0001   #Sheet thickness in cm (given entry in µm)
                
#========================================================================
            

        
            
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=24,column=2)
    button=Button(Fenster,text="Next",command=CV_Type_Chooser).grid(row=25,column=2)
    


# In[37]:


def CV_Type_Chooser():
    if model_a == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_a_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_b == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_b_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_c == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_c_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_d == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_d_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_e == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_e_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_f == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_f_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
    if model_g == 1:
        Start_End_Definer()
        if Statistical  == 0:
            model_g_ConvolutionInverter()
        if Statistical  == 1:
            DistFunc_Calcer()
   


# In[38]:


def Start_End_Definer():    
    
    if FIT == 1:
        
        global EEE
        global SSS
        global AsymmetrischesEnde
        
        StromarraysROH             = np.squeeze(np.transpose(Strom[0::,0:1:]))   
        UnkorrPotenzialarraysROH   = np.squeeze(np.transpose(Potenzial[0::,0:1:]))
        PotenzialarraysROH         = np.squeeze(np.transpose(Potenzial[0::,0:1:]))
         
        global Stromarrays
        global Potenzialarrays
                  
        Stromarrays = StromarraysROH
        
        UnkorrPotenzialarrays      = UnkorrPotenzialarraysROH[::] 
        Potenzialarrays            = (PotenzialarraysROH[::]   -  Ru*Stromarrays/1000000) - Ezero
        
        
        #SSS                = -np.abs(Potenzialarrays[0])
        #AsymmetrischesEnde = -np.abs(Potenzialarrays[-1])
        
        
        if DesiReactxx ==0:
            Stromarrays        = Stromarrays[::-1]
            Potenzialarrays    = Potenzialarrays[::-1]
            EEE                = np.abs(np.min(Potenzialarrays))    #Reduktion
            SSS                = -np.abs(Potenzialarrays[0])
            AsymmetrischesEnde = -np.abs(Potenzialarrays[-1])
            
        if DesiReactxx ==1:
            EEE                = np.abs(np.max(Potenzialarrays))    #Oxidation
            SSS                = -np.abs(Potenzialarrays[0])
            AsymmetrischesEnde = -np.abs(Potenzialarrays[-1])
            
        
    
        if LinBaseCorr == 1:
            global StartPotLinCorr
            global EndPotLinCorr
            global BaseSlope
            
            if DesiReactxx ==0:
                Potenzialarrays = Potenzialarrays[::-1]
                Stromarrays     = Stromarrays[::-1]
                HinScanEndindex = np.argmin(Potenzialarrays)
                LineStartIndex  = (np.abs(Potenzialarrays[-1:HinScanEndindex:-1]+Ezero-LineStart)).argmin()
                StartPotLinCorr = Potenzialarrays[LineStartIndex] #Startpot der Lin Korr vs. E_Zero
                LineEndIndex    = (np.abs(Potenzialarrays[-1:HinScanEndindex:-1]+Ezero-LineEnd)).argmin()
                EndPotLinCorr   = Potenzialarrays[LineEndIndex]   #Endpot der Lin Korr vs. E_Zero
                BaseSlope       = (Stromarrays[LineEndIndex]-Stromarrays[LineStartIndex])/(Potenzialarrays[LineEndIndex]-Potenzialarrays[LineStartIndex])
                Stromarrays     = (Stromarrays - Stromarrays[LineStartIndex])-BaseSlope*(Potenzialarrays-Potenzialarrays[LineStartIndex])
               
            
                
            if DesiReactxx ==1:
                HinScanEndindex = np.argmax(Potenzialarrays)
                LineStartIndex  = np.abs(Potenzialarrays[0:HinScanEndindex]+Ezero-LineStart).argmin()
                StartPotLinCorr = Potenzialarrays[LineStartIndex] #Startpot der Lin Korr vs. E_Zero
                LineEndIndex    = np.abs(Potenzialarrays[0:HinScanEndindex]+Ezero-LineEnd).argmin()
                EndPotLinCorr   = Potenzialarrays[LineEndIndex]   #Endpot der Lin Korr vs. E_Zero
                BaseSlope       = (Stromarrays[LineEndIndex]-Stromarrays[LineStartIndex])/(Potenzialarrays[LineEndIndex]-Potenzialarrays[LineStartIndex])
                Stromarrays     = (Stromarrays - Stromarrays[LineStartIndex])-BaseSlope*(Potenzialarrays-Potenzialarrays[LineStartIndex])
            
            
            
    


# In[39]:


def DistFunc_Calcer(): 

    global R_Array
    global P_Array
    
    
    y_Array = np.arange(0,100,0.01)
    P_Array = (DFP_c*DFP_b**(DFP_a/DFP_c) / gamma(DFP_a/DFP_c)) *(y_Array**(DFP_a-1) *np.exp(-DFP_b*y_Array**DFP_c))/average_d_or_r
    R_Array = average_d_or_r*y_Array
    
    
    IndexMax         =  np.asscalar(np.where(P_Array == np.max(P_Array))[0])
    Kurvenende       = (np.abs(P_Array[IndexMax::]-0.01*P_Array[IndexMax])).argmin()
    Array_Ab_Maximum = P_Array[IndexMax::]
    Endwert          = Array_Ab_Maximum[Kurvenende]
    EndIndex         = np.asscalar(np.where(P_Array == Endwert)[0])

    R_Array = R_Array[0:EndIndex] 
    P_Array = P_Array[0:EndIndex]
    
    
    global Probab_Array
    global R_Mean_Array  
    global RR_Mean_Array

    Probab_Array  = np.empty(DistFunc_Nodes) 
    R_Mean_Array  = np.empty(DistFunc_Nodes)
    RR_Mean_Array = np.empty(DistFunc_Nodes)
    
    DistFuncInterpolation = InterpolatedUnivariateSpline(R_Array,P_Array)
    R_Div = (R_Array[-1] - R_Array[0])/(1.0*float(DistFunc_Nodes))

    for i in range(len(Probab_Array)-1):
        Probab_Array[i] = quad(DistFuncInterpolation, R_Div*i,R_Div*(i+1) )[0]   #Wahscheinl. v. großem Segment
        R_Mean_Array[i] = R_Div*i                                                #Das sind die Intervallgrenzen
    
        R_Array_Resolved  = np.linspace(R_Div*i,R_Div*(i+1),50)
        PR_Indiv_Resolved = np.zeros(len(R_Array_Resolved))         #Indiv. Probability of further resolved Integral-Part
        for j in range(len(R_Array_Resolved)-1):
            #berechnen von Mittlerem R von jedem einzelnen großen Segment.
            #Jedes große Segment wird in 50 Teile geiteilt und darüber dann gewichtet integriert um mittleres R
            #von großem Segment zu finden. Das wird dann mit der Gewichtung des großen Segments am Gesamtanteil verrechnet.
            PR_Indiv_Resolved[j] = 0.5*(R_Array_Resolved[j]+R_Array_Resolved[j+1])*quad(DistFuncInterpolation, R_Array_Resolved[j],R_Array_Resolved[j+1] )[0]
        
        RR_Mean_Array[i] = (1/Probab_Array[i])*np.sum(PR_Indiv_Resolved)    #gewichteter Flächenanteil 
                                                                            #Wahrscheinlichkeit, dass mittleres R aus großem
                                                                            #Segment auftaucht.
    Probab_Array = Probab_Array[0:-1]
    #Das sind eigentlich die Abstände von zwei Fasermitten, ALSO EIGENTLCIH KEINE R!!!
    R_Array      = R_Array
    R_Mean_Array = R_Mean_Array[0:-1]
   

    DistFunc_User()
    
    if Show_Dist_Func == 1:
        Dist_Func_Shower()

    


# In[40]:


def DistFunc_User():
    global distance 
    global a
    
    if FIT == 1:
        Start_End_Definer()
        
    
    distance = 0.01    #to run calculation with default values before real calculation starts
    a        = 0.01    #to run calculation with default values before real calculation starts
    
    #HIER MÜSSEN ERST DIE CONVOLUTION INVERTERS GEZÜNDET WERDEN um len(yyy) zu bekommen
    if model_b == 1:
        model_b_ConvolutionInverter()
    if model_d == 1:
        model_d_ConvolutionInverter()
    if model_e == 1:
        model_e_ConvolutionInverter()
    if model_g == 1:
        model_g_ConvolutionInverter()
    
    
    global ChiSuperarray
    ChiSuperarray = np.zeros([len(yyy),len(R_Mean_Array)])
    
    
    #internal spherical finite diffusion --> calculation of area-weightning
    #----------------------------------------------------------------------------------------------------
    if model_g == 1:
        SUMMANDS    = np.zeros(len(Probab_Array))
        for i in range(len(Probab_Array)):
            SUMMANDS[i] = Probab_Array[i]*(RR_Mean_Array[i]/2.0)**2
        N_spheres = A/(4*np.pi*np.sum(SUMMANDS)) 
        
        INDIV_AREAS = np.zeros(len(Probab_Array))
        for i in range(len(Probab_Array)):
            INDIV_AREAS[i] = (N_spheres*Probab_Array[i])*(4*np.pi*(RR_Mean_Array[i]/2.0)**2)
    #----------------------------------------------------------------------------------------------------  
    
    #internal cylindrical finite diffusion --> calculation of area-weightning
    #----------------------------------------------------------------------------------------------------
    if model_e == 1:
        lunit = 0.001
        SUMMANDS    = np.zeros(len(Probab_Array))
        for i in range(len(Probab_Array)):
            SUMMANDS[i] = Probab_Array[i]*(RR_Mean_Array/2.0)[i]
        N_cyls = A/(2*np.pi*lunit*np.sum(SUMMANDS)) 
        
        INDIV_AREAS = np.zeros(len(Probab_Array))
        for i in range(len(Probab_Array)):
            INDIV_AREAS[i] = (N_cyls*Probab_Array[i])*(2*np.pi*lunit*(RR_Mean_Array/2.0)[i])
    #----------------------------------------------------------------------------------------------------  
     
        
    for i in range(len(Probab_Array)):   #is one shorter than RR_Mean_Array. Takes only the "useful" RR_Mean (not the last)
        
        #finite planar
        if model_b == 1:  
            distance = 0.5*RR_Mean_Array[i]*0.0001 - 0.5*sheet_thickness    #sheet distance/2 
        
        #internal cylindrical finite
        if model_e == 1:
            a = 0.5*RR_Mean_Array[i]*0.0001                                 #pore radii in µm
                
        #internal spherical finite
        if model_g == 1:
            a = 0.5*RR_Mean_Array[i]*0.0001                                 #pore radii in µm
        
        #external cylindrical finite
        if model_d == 1:
            #not necssary to define "a" here, because it will be re-defined in corresponding modell 
            distance = 0.5*RR_Mean_Array[i]*0.0001  - r                 #center to wall distances in µm

            
        INTERRUPT = 0
        
        if model_b == 1:
            if distance < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            model_b_ConvolutionInverter()
        if model_d == 1:
            if distance < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            model_d_ConvolutionInverter()
        if model_e == 1:
            if a < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            model_e_ConvolutionInverter()
        if model_g == 1:
            if a < 0.00001:
                Parameters_do_not_fit_warner()
                INTERRUPT = 1
                break
            model_g_ConvolutionInverter()
         
        for j in range(len(yyy)):
            if model_b == 1 or model_d == 1:
                ChiSuperarray[j,i] = yyy[j] *np.asscalar(Probab_Array[i])   #weightning by d-probability (planar fin and cylextfin)
            if model_e == 1 or model_g == 1:
                ChiSuperarray[j,i] = yyy[j] *INDIV_AREAS[i]/A       #weightning by area for internal cyl and internal sph

    
    global Chi_ges_Array
    Chi_ges_Array = np.empty(len(yyy))
    for i in range(len(yyy)):
        for j in range(len(Probab_Array)):
            Chi_ges_Array[i] = np.sum(ChiSuperarray[i,::])
            
    #====================================================================================================
    #Plotten von simuliertem
    #====================================================================================================
    if INTERRUPT == 0:
        
        root = Toplevel()
        root.title("Simulated CV")

        if Xi_i > -6.9:
            Incorrect_Range_Warner()

        if FIT == 0:
            global PotCalc
            global CurrCalc     
            PotCalc  = R*T/(n*F)*xxx 
            CurrCalc = (n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000*ChiSuperarray

            Abbildung = Figure(figsize=(5, 4), dpi=100)
            b = Abbildung.add_subplot(111)

            if Just_Current_Template == 0:
                if Show_Individuals == 0:
                    b.plot(PotCalc[0:-1],Chi_ges_Array[0:-1]*(n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000, color = 'r')     
                if Show_Individuals == 1:
                    b.plot(PotCalc[0:-1],CurrCalc[0:-1,::])
                b.axvline(0,color='k')
                b.axhline(0,color='k')
                to = np.abs(EEE)
                b.set_xlim([-(to+0.05*to),(to+0.05*to)])    
                b.set_xlabel('E vs. E_zero. / V', fontsize=12)
                b.set_ylabel('I / mA', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b.spines[axis].set_linewidth(2)
                    b.spines[axis].set_color('k')
                b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)

            if Just_Current_Template == 1:
                if Show_Individuals == 0:
                    b.plot(xxx[0:-1],Chi_ges_Array[0:-1], color = 'r') 
                    if xxx[0]<xxx[1]:
                        b.annotate(u'\u03A8' '=%5.3f'  % max(Chi_ges_Array[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                    if xxx[0]>xxx[1]:
                        b.annotate(u'\u03A8' '=%5.3f'  % min(Chi_ges_Array[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                if Show_Individuals == 1:
                    b.plot(xxx[0:-1],ChiSuperarray[0:-1,::])
                b.axvline(0,color='k')
                b.axhline(0,color='k')
                to = n*F*np.abs(EEE)/(R*T)
                b.set_xlim([-(to+0.1*to),(to+0.1*to)])    
                b.set_xlabel(u'\u03BE', fontsize=12)
                b.set_ylabel(u'\u03C7', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b.spines[axis].set_linewidth(2)
                    b.spines[axis].set_color('k')
                b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)

            canvas = FigureCanvasTkAgg(Abbildung, master=root)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            toolbar = NavigationToolbar2TkAgg(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)


        if FIT == 1:
            PotCorr  = R*T/(n*F) 
            CurrCorr = (n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000000
            global Pot
            global Curr
            Pot              = Potenzialarrays
            Curr             = Stromarrays
            PotCalc          = xxx*PotCorr
            Chi_ges_Array    = CurrCorr*Chi_ges_Array       #Unterschied zu Vorher :)
            NormLengthPot    = (1.0/float(len(Pot)))
            NormLengthxxx    = (1.0/float(len(PotCalc)))
            NormArrayPot     = np.arange(0,1,NormLengthPot)
            NormArrayxxx     = np.arange(0,1,NormLengthxxx)
            Fitinterpolation = InterpolatedUnivariateSpline(NormArrayxxx,Chi_ges_Array,k=3)
            InterpolatedCurr = np.empty(len(NormArrayPot))
            for i in range(len(NormArrayPot)):
                InterpolatedCurr[i] = Fitinterpolation(NormArrayPot[i])
            Abweichungsarray = np.empty(len(Curr)-1)
            for i in range(len(Curr)-1):
                Abweichungsarray[i] = (((Curr[i] - InterpolatedCurr[i])/np.max(Curr))**2)/(float(len(Curr))-1)
            Varianz            = np.sum(Abweichungsarray)
            global Standardabweichung
            Standardabweichung = Varianz**0.5
            Abbildung = Figure(figsize=(5, 4), dpi=100)
            b = Abbildung.add_subplot(111)
            b.plot(PotCalc[1:-1:5],0.001*Chi_ges_Array[1:-1:5],color='r', marker = '.', linestyle = '')
            b.plot(PotCalc[:-1:5],0.001*Chi_ges_Array[:-1:5],color='r')
            b.plot(Pot[:-1:],0.001*Curr[:-1:],color='k')       
            b.annotate(u'\u03C3' '=%5.3f'  % Standardabweichung , xy=(0.05, 0.83), xycoords='axes fraction')
            b.axvline(0,color='k')
            b.axhline(0,color='k')
            b.set_xlabel('E vs. E_zero. / V', fontsize=12)
            b.set_ylabel('I / mA', fontsize=12)
            for axis in ['top','bottom','left','right']:
                b.spines[axis].set_linewidth(2)
                b.spines[axis].set_color('k')
            b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            b.set_xlim([-(np.abs(SSS+0.05*SSS)),(np.abs(EEE+0.05*EEE))])
            canvas = FigureCanvasTkAgg(Abbildung, master=root)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            toolbar = NavigationToolbar2TkAgg(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)



        if AstxtSaver == 1:
            Simulated_CV_as_txt_saver()

    


# In[41]:


def Dist_Func_Shower():   
    root = Toplevel()
    root.title("Distance distribution-function")
    
    Abbildung = Figure(figsize=(5, 4), dpi=100)
    b = Abbildung.add_subplot(111)
    for i in range(len(R_Mean_Array)):
        #mal 0.5. weil R_Mean_Array quasi durchmesser sind
        b.axvline(R_Mean_Array[i], color = 'k', linewidth = 0.5, linestyle ='-')   #das sind nicht die Mittelwerte, sondern die Stützstellen
        b.axvline(RR_Mean_Array[i], color = 'k', linewidth = 0.5, linestyle ='--') #Das sind die eigentlichen Mittelwerte
        
    #mal 0.5 weil R_Mean_Array quasi urchmesser sind und mal 2, damit Integral passt
    b.plot(R_Array,P_Array, color = 'b')
    b.set_xlabel('midpoint dist. x / micrometer', fontsize=12)
    b.set_ylabel('P(midpoint dist. x) / micrometer^-1', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b.spines[axis].set_linewidth(2)
        b.spines[axis].set_color('k')
    b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    


# In[42]:


#define all Talbot-Inverters

def TalbotPlanarSemiHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanarSemiHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotPlanarSemiBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanarSemiBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotPlanarFinitHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanarFinitHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotPlanarFinitBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanarFinitBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotZylSemiHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylSemiHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotZylSemiBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylSemiBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   


def TalbotPlanFinHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanFinHin(z)/z)*dz;
            
    return ((h/(2j*pi))*ans).real 

def TalbotPlanFinBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(PlanFinBack(z)/z)*dz;
            
    return ((h/(2j*pi))*ans).real 

def TalbotZylFinHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylFinHin(z)/z)*dz;
        
    return ((h/(2j*pi))*ans).real   

def TalbotZylFinBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylFinBack(z)/z)*dz;
        
    return ((h/(2j*pi))*ans).real   

def TalbotZylIntFinHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylIntFinHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotZylIntFinBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(ZylIntFinBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotSpherExtSemiHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(SpherExtSemiHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotSpherExtSemiBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(SpherExtSemiBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   


def TalbotSpherIntFinHin(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(SpherIntFinHin(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   

def TalbotSpherIntFinBack(t,N):
    h = 2*pi/N;
    shift = 0.0;
    ans   = 0.0;
    if t == 0:
        print "ERROR:   Inverse transform can not be calculated for t=0"
        return ("Error");
    for k in range(0,N):
        theta  = -pi + (k+1./2)*h;
        z      = shift + N/t*(0.5017*theta*cot(0.6407*theta) - 0.6122 + 0.2645j*theta); 
        dz     = N/t*(-0.5017*0.6407*theta*(csc(0.6407*theta)**2)+0.5017*cot(0.6407*theta)+0.2645j);
        ans    = ans + exp(z*t)*(SpherIntFinBack(z))/z*dz;
            
    return ((h/(2j*pi))*ans).real   


# In[43]:


#Define Convolution functions in laplace Domain

def PlanarSemiHin(s):
    return 1/s**0.5

def PlanarSemiBack(s):
    return 1/s**0.5

def PlanarFinitHin(s):
    return coth(distance*(s/D_f)**0.5)/s**0.5

def PlanarFinitBack(s):
    return coth(distance*(s/D_b)**0.5)/s**0.5

def ZylSemiHin(s):
    if np.isnan(kv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_f**0.5)*s**0.5) )) == False:
        return kv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_f**0.5)*s**0.5) )
    if np.isnan(kv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_f**0.5)*s**0.5) )) == True:
        return 1/s**0.5

def ZylSemiBack(s):
    if np.isnan(kv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_b**0.5)*s**0.5) )) == False:
        return kv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_b**0.5)*s**0.5) )
    if np.isnan(kv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_b**0.5)*s**0.5) )) == True:
        return 1/s**0.5
def PlanFinHin(s):
    return coth((distance/D_f**0.5)*s**0.5)/s**0.5

def PlanFinBack(s):
    return coth((distance/D_b**0.5)*s**0.5)/s**0.5

def ZylSemiHin(s):
    if np.isnan(kv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_f**0.5)*s**0.5) )) == False:
        return kv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_f**0.5)*s**0.5) )
    if np.isnan(kv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_f**0.5)*s**0.5) )) == True:
        return 1/s**0.5

def ZylSemiBack(s):
    if np.isnan(kv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_b**0.5)*s**0.5) )) == False:
        return kv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_b**0.5)*s**0.5) )
    if np.isnan(kv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * kv(1,(a/D_b**0.5)*s**0.5) )) == True:
        return 1/s**0.5

def ZylFinHin(s):
    return (kv(0,(a/D_f**0.5)*s**0.5)*iv(1,(d/D_f**0.5)*s**0.5) + iv(0,(a/D_f**0.5)*s**0.5)*kv(1,(d/D_f**0.5)*s**0.5))/(s**0.5 *(kv(1,(a/D_f**0.5)*s**0.5)*iv(1,(d/D_f**0.5)*s**0.5) - iv(1,(a/D_f**0.5)*s**0.5)*kv(1,(d/D_f**0.5)*s**0.5)))

def ZylFinBack(s):
    return (kv(0,(a/D_b**0.5)*s**0.5)*iv(1,(d/D_b**0.5)*s**0.5) + iv(0,(a/D_b**0.5)*s**0.5)*kv(1,(d/D_b**0.5)*s**0.5))/(s**0.5 *(kv(1,(a/D_b**0.5)*s**0.5)*iv(1,(d/D_b**0.5)*s**0.5) - iv(1,(a/D_b**0.5)*s**0.5)*kv(1,(d/D_b**0.5)*s**0.5)))

def ZylIntFinHin(s):
     if np.isnan(iv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_f**0.5)*s**0.5) )) == False:
        return iv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_f**0.5)*s**0.5) )
     if np.isnan(iv(0,(a/D_f**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_f**0.5)*s**0.5) )) == True:
        return 1/s**0.5

def ZylIntFinBack(s):
     if np.isnan(iv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_b**0.5)*s**0.5) )) == False:
        return iv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_b**0.5)*s**0.5) )
     if np.isnan(iv(0,(a/D_b**0.5)*s**0.5)/(s**0.5 * iv(1,(a/D_b**0.5)*s**0.5) )) == True:
        return 1/s**0.5

def SpherExtSemiHin(s):
    return 1/(s**0.5 + D_f**0.5/a)

def SpherExtSemiBack(s):
    return 1/(s**0.5 + D_b**0.5/a)

def SpherIntFinHin(s):
    return 1/(s**0.5 *coth(a*(s/D_f)**0.5) - D_f**0.5/a)

def SpherIntFinBack(s):
    return 1/(s**0.5 *coth(a*(s/D_b)**0.5) - D_b**0.5/a)




# In[44]:


#Invert the functions of Interest

def model_a_ConvolutionInverter():
    if D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    
    if D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                      = np.logspace(-8,5,300)
        ft_Array_PlanarSemiHin       = np.empty(len(t_Array)) 
        ft_Array_PlanarSemiBack      = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_PlanarSemiHin[i]     = TalbotPlanarSemiHin(float(t_Array[i]),24)
            ft_Array_PlanarSemiBack[i]    = TalbotPlanarSemiBack(float(t_Array[i]),24)
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_PlanarSemiHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_PlanarSemiBack, k=3)
        
        CV_Calculator()
    

    
def model_b_ConvolutionInverter():
    #also in statistical case distance will be passed as distance
    if distance < 0.00001 or distance > 10.0 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if distance >= 0.00001 and distance <= 10.0 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                       = np.logspace(-8,5,300)
        ft_Array_PlanarFinitHin       = np.empty(len(t_Array)) 
        ft_Array_PlanarFinitBack      = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_PlanarFinitHin[i]     = TalbotPlanarFinitHin(float(t_Array[i]),24)
            ft_Array_PlanarFinitBack[i]    = TalbotPlanarFinitBack(float(t_Array[i]),24)   
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_PlanarFinitHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_PlanarFinitBack, k=3)
        CV_Calculator()
    
    
    
def model_c_ConvolutionInverter():
    global a
    a = r
    
    if a>10 or a<0.00001 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
        
    if a<=10 and a>=0.00001 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array              = np.logspace(-8,5,300)
        ft_Array_ZylSemiHin  = np.empty(len(t_Array)) 
        ft_Array_ZylSemiBack = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_ZylSemiHin[i]     = TalbotZylSemiHin(float(t_Array[i]),24)
            ft_Array_ZylSemiBack[i]    = TalbotZylSemiBack(float(t_Array[i]),24)
        #----------------------------------------------------------------------------------------
        #Interpolationen finden
        #----------------------------------------------------------------------------------------
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_ZylSemiHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_ZylSemiBack, k=3)
        CV_Calculator()
    

    
def model_d_ConvolutionInverter():
    global a
    a = r

    global d 
    d = a + distance
    
    if a>10 or a<0.00001 or distance>10 or distance<0.00001 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if a<=10 and a>=0.00001 and distance<=10 and distance>=0.00001 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array              = np.logspace(-8,5,300)
        ft_Array_ZylSemiHin  = np.empty(len(t_Array)) 
        ft_Array_ZylFinHin   = np.empty(len(t_Array)) 
        ft_Array_PlanFinHin  = np.empty(len(t_Array)) 
        ft_Array_ZylSemiBack = np.empty(len(t_Array)) 
        ft_Array_ZylFinBack  = np.empty(len(t_Array)) 
        ft_Array_PlanFinBack = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_ZylSemiHin[i]     = TalbotZylSemiHin(float(t_Array[i]),24)
            ft_Array_PlanFinHin[i]     = TalbotPlanFinHin(float(t_Array[i]),24)
            ft_Array_ZylFinHin[i]      = TalbotZylFinHin(float(t_Array[i]),24)
            ft_Array_ZylSemiBack[i]    = TalbotZylSemiBack(float(t_Array[i]),24)
            ft_Array_PlanFinBack[i]    = TalbotPlanFinBack(float(t_Array[i]),24)
            ft_Array_ZylFinBack[i]     = TalbotZylFinBack(float(t_Array[i]),24)
        #--------------------------------------------------------------------------------------------    
        ft_Array_ZylSemi_revHin  = ft_Array_ZylSemiHin[::-1]
        ft_Array_ZylFin_revHin   = ft_Array_ZylFinHin[::-1] 
        ft_Array_PlanFin_revHin  = ft_Array_PlanFinHin[::-1]
        ft_Array_ZylSemi_revBack = ft_Array_ZylSemiBack[::-1]
        ft_Array_ZylFin_revBack  = ft_Array_ZylFinBack[::-1] 
        ft_Array_PlanFin_revBack = ft_Array_PlanFinBack[::-1]
        Schalter1Hin  = 0
        Schalter2Hin  = 0
        Schalter1Back = 0
        Schalter2Back = 0
        GesamtHin  = np.empty(len(t_Array))
        GesamtBack = np.empty(len(t_Array))
        #---------------------------------------------------------------------------------------
        #Hin
        #---------------------------------------------------------------------------------------
        for i in range(len(t_Array)):
            if Schalter1Hin == 0:
                if a <=0.01:
                    if distance >= 0.0001:
                        GesamtHin[i] = ft_Array_ZylFin_revHin[i]
                        if np.abs(ft_Array_PlanFin_revHin[i]- ft_Array_ZylFin_revHin[i]) <= 10**(-5):
                            Schalter1Hin = 1
                    if distance < 0.0001:
                        GesamtHin[i] = ft_Array_ZylFin_revHin[i]
                        if np.abs(ft_Array_PlanFin_revHin[i]- ft_Array_ZylFin_revHin[i]) <= 10**(-4):
                            Schalter1Hin = 1
                if a > 0.01:
                    GesamtHin[i] = ft_Array_ZylFin_revHin[i]
                    if np.abs(ft_Array_PlanFin_revHin[i]- ft_Array_ZylFin_revHin[i]) <= 10**(-3):
                        Schalter1Hin = 1
                if a > 0.1:
                    GesamtHin[i] = ft_Array_ZylFin_revHin[i]
                    if np.abs(ft_Array_PlanFin_revHin[i]- ft_Array_ZylFin_revHin[i]) <= 10**(-2):
                        Schalter1Hin = 1
                if a > 1:
                    GesamtHin[i] = ft_Array_ZylFin_revHin[i]
                    if np.abs(ft_Array_PlanFin_revHin[i]- ft_Array_ZylFin_revHin[i]) <= 10**(-1):
                        Schalter1Hin = 1
            if Schalter1Hin == 1:
                GesamtHin[i] = ft_Array_PlanFin_revHin[i]
        for i in range(len(t_Array)):
            if Schalter2Hin == 0:
                if np.abs(GesamtHin[i] - ft_Array_ZylSemi_revHin[i]) <= 10**(-5):
                    Schalter2Hin = 1
            if Schalter2Hin == 1:
                GesamtHin[i] = ft_Array_ZylSemi_revHin[i]  
        #---------------------------------------------------------------------------------------
        #Back
        #---------------------------------------------------------------------------------------       
            if Schalter1Back == 0:
                if a <=0.01:
                    if distance >= 0.0001:
                        GesamtBack[i] = ft_Array_ZylFin_revBack[i]
                        if np.abs(ft_Array_PlanFin_revBack[i]- ft_Array_ZylFin_revBack[i]) <= 10**(-5):
                            Schalter1Back = 1
                    if distance < 0.0001:
                        GesamtBack[i] = ft_Array_ZylFin_revBack[i]
                        if np.abs(ft_Array_PlanFin_revBack[i]- ft_Array_ZylFin_revBack[i]) <= 10**(-4):
                            Schalter1Back = 1
                if a > 0.01:
                    GesamtBack[i] = ft_Array_ZylFin_revBack[i]
                    if np.abs(ft_Array_PlanFin_revBack[i]- ft_Array_ZylFin_revBack[i]) <= 10**(-3):
                        Schalter1Back = 1
                if a > 0.1:
                    GesamtBack[i] = ft_Array_ZylFin_revBack[i]
                    if np.abs(ft_Array_PlanFin_revBack[i]- ft_Array_ZylFin_revBack[i]) <= 10**(-2):
                        Schalter1Back = 1
                if a > 1:
                    GesamtBack[i] = ft_Array_ZylFin_revBack[i]
                    if np.abs(ft_Array_PlanFin_revBack[i]- ft_Array_ZylFin_revBack[i]) <= 10**(-1):
                        Schalter1Back = 1
            if Schalter1Back == 1:
                GesamtBack[i] = ft_Array_PlanFin_revBack[i]
        for i in range(len(t_Array)):
            if Schalter2Back == 0:
                if np.abs(GesamtBack[i] - ft_Array_ZylSemi_revBack[i]) <= 10**(-5):
                    Schalter2Back = 1
            if Schalter2Back == 1:
                GesamtBack[i] = ft_Array_ZylSemi_revBack[i]     
        #----------------------------------------------------------------------------------------
        #Interpolationen finden
        #----------------------------------------------------------------------------------------
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,GesamtHin[::-1], k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,GesamtBack[::-1], k=3)
        CV_Calculator()
    
    
    
def model_e_ConvolutionInverter():
    global a
    if Statistical == 0:
        a = r
    
    if a>10 or a<0.00001 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if a<=10 and a>=0.00001 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                = np.logspace(-8,5,300)
        ft_Array_ZylIntFinHin  = np.empty(len(t_Array)) 
        ft_Array_ZylIntFinBack = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_ZylIntFinHin[i]     = TalbotZylIntFinHin(float(t_Array[i]),24)
            ft_Array_ZylIntFinBack[i]    = TalbotZylIntFinBack(float(t_Array[i]),24)
        #----------------------------------------------------------------------------------------
        #Interpolationen finden
        #----------------------------------------------------------------------------------------
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_ZylIntFinHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_ZylIntFinBack, k=3)
        CV_Calculator()
    


def model_f_ConvolutionInverter():
    global a
    a = r
    
    if a>10 and a<0.00001 and D_f <0.000001 and D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if a<=10 and a>=0.00001 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                      = np.logspace(-8,5,300)
        ft_Array_SphereExtSemiHin    = np.empty(len(t_Array)) 
        ft_Array_SphereExtSemiBack   = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_SphereExtSemiHin[i]     = TalbotSpherExtSemiHin(float(t_Array[i]),24)
            ft_Array_SphereExtSemiBack[i]    = TalbotSpherExtSemiBack(float(t_Array[i]),24)          
        #----------------------------------------------------------------------------------------
        #Interpolationen finden
        #----------------------------------------------------------------------------------------
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_SphereExtSemiHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_SphereExtSemiBack, k=3)
        CV_Calculator()
    
    

def model_g_ConvolutionInverter():
    global a
    if Statistical == 0:
        a = r
    
    if a>10 or a<0.00001 or D_f <0.000001 or D_b <0.000001:
        Parameters_do_not_fit_warner()
    
    if a<=10 and a>=0.00001 and D_f >=0.000001 and D_b >=0.000001:
        #--------------------------------------------------------------------------------------------
        t_Array                     = np.logspace(-8,5,300)
        ft_Array_SphereIntFinHin    = np.empty(len(t_Array)) 
        ft_Array_SphereIntFinBack   = np.empty(len(t_Array)) 
        for i in range(len(t_Array)):
            ft_Array_SphereIntFinHin[i]     = TalbotSpherIntFinHin(float(t_Array[i]),24)
            ft_Array_SphereIntFinBack[i]    = TalbotSpherIntFinBack(float(t_Array[i]),24)  
        global CV_InterpolationHin
        CV_InterpolationHin = InterpolatedUnivariateSpline(t_Array,ft_Array_SphereIntFinHin, k=3)
        global CV_InterpolationBack
        CV_InterpolationBack = InterpolatedUnivariateSpline(t_Array,ft_Array_SphereIntFinBack, k=3)
        CV_Calculator()
    
    


# In[45]:


def CV_Calculator():
    
    if Statistical == 0:
        root = Toplevel()
        root.title("Simulated CV")
    
    global Kp
    global Kf
    global p
    global f
    p               = 0.0
    Kp              = 1000000000000.0
    Kf              = 0.00000000001
    f               = 0.0
    if Precedeq == 1:
        p               = kp + kmp
        Kp              = kp/kmp
    if Succeedeq == 1:
        Kf              = kf/kmf
        f               = kf + kmf
    global t_Array_cont
    global Xi_i
    if FIT == 0:
        ScaRa           = Scanrate
        delta_Xi        = 0.1
        if DeltaXiModder == 1:
            delta_Xi = ModDeltXi
        E_f             = EEE
        delta_E         = delta_Xi*R*T/(n*F)
        Num_of_E_forw   = int(2*E_f/delta_E)
        E_i             = E_f - (Num_of_E_forw+1)*delta_E
        delta_t         = delta_Xi*R*T/(n*F*ScaRa) 
        tges            = (np.abs(E_i)+E_f+(E_f-(E_i)))/ScaRa
        t_Array_cont    = np.arange(0,tges+delta_t,delta_t)
        sigma           = n*F*ScaRa/(R*T)
        Xi_i            = n*F*E_i/(R*T)   #initial Xi
        Xi_f            = n*F*(E_f)/(R*T)  #final Xi
        Xi_Asym         = n*F*(E_i)/(R*T)
        Xi_Array_hin    = np.arange(Xi_i,Xi_f,delta_Xi) 
        Xi_Array_back   = np.arange(Xi_f,Xi_Asym-delta_Xi,-delta_Xi) 
        Xi_Array        = np.concatenate([Xi_Array_hin,Xi_Array_back])
        Exp_Array       = np.exp(-Xi_Array)
        Preced_Array    = np.exp(-p*t_Array_cont)
        Follow_Array    = np.exp(-f*t_Array_cont)
        Kinetik_Array   = D_f**0.5/(kzero*Exp_Array**(-alpha))
        
    if FIT == 1:
        ScaRa           = Scanrate
        delta_Xi        = 0.1
        if DeltaXiModder == 1:
            delta_Xi = ModDeltXi
        E_f             = EEE
        delta_E         = delta_Xi*R*T/(n*F)
        Num_of_E_forw   = int((E_f-SSS)/delta_E)
        E_i             = E_f - (Num_of_E_forw+1)*delta_E
        delta_t         = delta_Xi*R*T/(n*F*ScaRa) 
        tges            = (np.abs(E_i)+E_f+(E_f-(E_i)))/ScaRa
        t_Array_cont    = np.arange(0,tges+delta_t,delta_t)
        sigma           = n*F*ScaRa/(R*T)
        Xi_i            = n*F*E_i/(R*T)                   #initial Xi
        Xi_f            = n*F*(E_f)/(R*T)                 #final Xi
        Xi_Asym         = n*F*(AsymmetrischesEnde)/(R*T)  #Asymmetrisch endendes Xi
        Xi_Array_hin    = np.arange(Xi_i,Xi_f,delta_Xi) 
        Xi_Array_back   = np.arange(Xi_f,Xi_Asym-delta_Xi,-delta_Xi) 
        Xi_Array        = np.concatenate([Xi_Array_hin,Xi_Array_back])
        Exp_Array       = np.exp(-Xi_Array)
        Preced_Array    = np.exp(-p*t_Array_cont)
        Follow_Array    = np.exp(-f*t_Array_cont)
        Kinetik_Array   = D_f**0.5/(kzero*Exp_Array**(-alpha))
        
    
    Fin_Kin_p_Array = 1/(1 + kfin*Exp_Array**(-alpha))      #hier mit minus, weil Exp-Array schon umgedreht ist!
    Fin_Kin_s_Array = 1/(1 + kfin*Exp_Array**((1-alpha)))
     
    
    

    #-------------------------------------------------------------------------------------------------- 
    #Jetzt CV berechnen
    #--------------------------------------------------------------------------------------------------
    global WuFuInt_p
    WuFuInt_p  = np.empty(len(t_Array_cont))
    for i in range(len(t_Array_cont)):
        WuFuInt_p[i] = np.pi**0.5 *CV_InterpolationHin(t_Array_cont[i])
    global WuFuInt_f
    WuFuInt_f  = np.empty(len(t_Array_cont))
    for i in range(len(t_Array_cont)):
        WuFuInt_f[i] = np.pi**0.5 *CV_InterpolationBack(t_Array_cont[i])
    
    Chi_Array  = np.zeros(len(Xi_Array))

    for i in range(len(Xi_Array)-1):
        Summandenarray_p_GG = np.zeros(i)
        Summandenarray_f_GG = np.zeros(i)
        Summandenarray_p    = np.zeros(i)
        Summandenarray_f    = np.zeros(i)
        for j in range(i):
            Summandenarray_p_GG[j]   = (WuFuInt_p[i-j+1]-WuFuInt_p[i-j])*Chi_Array[j]*Preced_Array[i-j]
            Summandenarray_f_GG[j]   = (WuFuInt_f[i-j+1]-WuFuInt_f[i-j])*Chi_Array[j]*Follow_Array[i-j]
            Summandenarray_p[j]      = (WuFuInt_p[i-j+1]-WuFuInt_p[i-j])*Chi_Array[j]
            Summandenarray_f[j]      = (WuFuInt_f[i-j+1]-WuFuInt_f[i-j])*Chi_Array[j]

        Chi_Array[i] = (Fin_Kin_p_Array[i]*(1/(1+Kp))*(Kp - Kp*np.sum(Summandenarray_p) - np.sum(Summandenarray_p_GG)) -(Fin_Kin_s_Array[i]*(D_f/D_b)**0.5)*(Exp_Array[i]/(1+Kf))*(np.sum(Summandenarray_f) + Kf*np.sum(Summandenarray_f_GG))  )/(  Kinetik_Array[i] +  Fin_Kin_p_Array[i]*WuFuInt_p[1]+ (Fin_Kin_s_Array[i]*(D_f/D_b)**0.5)*WuFuInt_f[1]*Exp_Array[i])

    #-----------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------------
    global xxx
    global yyy
    xxx = (Xi_Array)
    yyy   = np.pi**0.5*(Chi_Array[::])/sigma**0.5
     
    
    
    if FIT == 0:
        if Cathsweeper == 1:
            xxx = -(Xi_Array)
            yyy = -(np.pi**0.5*(Chi_Array)/sigma**0.5)
    if FIT ==1:
        if DesiReactxx == 0:
            xxx = -(Xi_Array)
            yyy = -(np.pi**0.5*(Chi_Array)/sigma**0.5)
            
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON , nur wenn nicht statistisch--> der hat seinen eigenen Plotter
    #-------------------------------------------------------------------------------------------------------------------------
    
    if Statistical == 0:
        if Xi_i > -6.9:
            Incorrect_Range_Warner()
            
        if FIT == 0:
            global PotCalc
            global CurrCalc     
            PotCalc  = R*T/(n*F)*xxx 
            CurrCalc = (n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000*yyy
            Abbildung = Figure(figsize=(5, 4), dpi=100)
            b = Abbildung.add_subplot(111)

            if Just_Current_Template == 0:
                b.plot(PotCalc[:-1:],CurrCalc[:-1:],color='r')
                if Cathsweeper == 0:
                    b.annotate('I_p / mA' '=%5.3f'  % max(CurrCalc[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                if Cathsweeper == 1:
                    b.annotate('I_p / mA' '=%5.3f'  % min(CurrCalc[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                b.axvline(0,color='k')
                b.axhline(0,color='k')
                to = np.abs(EEE)
                b.set_xlim([-(to+0.05*to),(to+0.05*to)])    
                b.set_xlabel('E vs. E_zero. / V', fontsize=12)
                b.set_ylabel('I / mA', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b.spines[axis].set_linewidth(2)
                    b.spines[axis].set_color('k')
                b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)

            if Just_Current_Template == 1:
                b.plot(xxx[:-1:],yyy[:-1:],color='r')
                if Cathsweeper == 0:
                    b.annotate(u'\u03A8' '=%5.3f'  % max(yyy[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                if Cathsweeper == 1:
                    b.annotate(u'\u03A8' '=%5.3f'  % min(yyy[:-1:]) , xy=(0.05, 0.83), xycoords='axes fraction')
                b.axvline(0,color='k')
                b.axhline(0,color='k')
                to = n*F*np.abs(EEE)/(R*T)
                b.set_xlim([-(to+0.1*to),(to+0.1*to)])    
                b.set_xlabel(u'\u03BE', fontsize=12)
                b.set_ylabel(u'\u03C7', fontsize=12)
                for axis in ['top','bottom','left','right']:
                    b.spines[axis].set_linewidth(2)
                    b.spines[axis].set_color('k')
                b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)

            canvas = FigureCanvasTkAgg(Abbildung, master=root)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            toolbar = NavigationToolbar2TkAgg(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)


        if FIT == 1:
            PotCorr  = R*T/(n*F) 
            CurrCorr = (n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000000
            global Pot
            global Curr
            Pot  = Potenzialarrays
            Curr = Stromarrays
            xxx  = xxx*PotCorr
            yyy  = yyy*CurrCorr
            NormLengthPot    = (1.0/len(Pot))
            NormLengthxxx    = (1.0/len(xxx))
            NormArrayPot     = np.arange(0,1,NormLengthPot)
            NormArrayxxx     = np.arange(0,1,NormLengthxxx)
            Fitinterpolation = InterpolatedUnivariateSpline(NormArrayxxx,yyy,k=3)
            InterpolatedCurr = np.empty(len(NormArrayPot))
            for i in range(len(NormArrayPot)):
                InterpolatedCurr[i] = Fitinterpolation(NormArrayPot[i])
            Abweichungsarray = np.empty(len(Curr)-1)
            for i in range(len(Curr)-1):
                Abweichungsarray[i] = (((Curr[i] - InterpolatedCurr[i])/np.max(Curr))**2)/(float(len(Curr))-1)
            Varianz            = np.sum(Abweichungsarray)
            global Standardabweichung
            Standardabweichung = Varianz**0.5
            Abbildung = Figure(figsize=(5, 4), dpi=100)
            b = Abbildung.add_subplot(111)
            b.plot(xxx[1:-1:5],0.001*yyy[1:-1:5],color='r', marker ='.', linestyle = '')   
            b.plot(xxx[:-1:],0.001*yyy[:-1:],color='r')   
            b.plot(Pot[:-1:],0.001*Curr[:-1:],color='k')  
            b.annotate(u'\u03C3' '=%5.3f'  % Standardabweichung , xy=(0.05, 0.83), xycoords='axes fraction')
            b.axvline(0,color='k')
            b.axhline(0,color='k')
            b.set_xlabel('E vs. E_zero. / V', fontsize=12)
            b.set_ylabel('I / mA', fontsize=12)
            for axis in ['top','bottom','left','right']:
                b.spines[axis].set_linewidth(2)
                b.spines[axis].set_color('k')
            b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
            b.set_xlim([-(np.abs(SSS+0.05*SSS)),(np.abs(EEE+0.05*EEE))])
            canvas = FigureCanvasTkAgg(Abbildung, master=root)
            canvas.show()
            canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
            toolbar = NavigationToolbar2TkAgg(canvas, root)
            toolbar.update()
            canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)

        if AstxtSaver == 1:
            Simulated_CV_as_txt_saver()

      


# In[46]:


def Incorrect_Range_Warner():
    Fenster = Toplevel()                                                         
    Fenster.title("Warning! Incorrect numerical range.")                         
    Fenster.geometry("400x200")
    Te = Text(Fenster, height=7, width=50)
    Te.pack()
    Te.insert(END, "Potential limits are not sufficiently\nfar away from E_zero. It should be\nat least\nXi_initial = -7 = nF(E_i-E_zero)/RT \nNevertheless, CV got calculated.")


# In[47]:


def Simulated_CV_as_txt_saver():
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")

        fi.write("n")
        fi.write("\t")
        fi.write(str(n))
        fi.write("\n")
        fi.write("T in [K]")
        fi.write("\t")
        fi.write(str(T))
        fi.write("\n")
        if Unequal_D ==0:
            fi.write("Common D in [cm^2/s]")
            fi.write("\t")
            fi.write(str(D_f))
            fi.write("\n")
        if Unequal_D ==1:
            fi.write("D_foreward in [cm^2/s]")
            fi.write("\t")
            fi.write(str(D_f))
            fi.write("\n")
            fi.write("D_backward in [cm^2/s]")
            fi.write("\t")
            fi.write(str(D_b))
            fi.write("\n")
        fi.write("Scanrate in [V/s]")
        fi.write("\t")
        fi.write(str(Scanrate))
        fi.write("\n")
        fi.write("alpha_forward")
        fi.write("\t")
        fi.write(str(alpha))
        fi.write("\n")
        fi.write("k_zero in [cm/s]")
        fi.write("\t")
        fi.write(str(kzero))
        fi.write("\n")
        if FinKin ==1:
            fi.write("k_het. max. in [cm/s]")
            fi.write("\t")
            fi.write(str(kzero/kfin))
            fi.write("\n")
        
        fi.write("A in [cm^2]")
        fi.write("\t")
        fi.write(str(A))
        fi.write("\n")
        fi.write("c in [mol/cm^3]")
        fi.write("\t")
        fi.write(str(c))
        fi.write("\n")
        
        
       
        if Precedeq == 1:
            fi.write("k_p [1/s]")
            fi.write("\t")
            fi.write(str(kp))
            fi.write("\n")
            fi.write("k_-p [1/s]")
            fi.write("\t")
            fi.write(str(kmp))
            fi.write("\n")
            fi.write("Kp")
            fi.write("\t")
            fi.write(str(Kp))
            fi.write("\n")
            fi.write("p")
            fi.write("\t")
            fi.write(str(p))
            fi.write("\n")
            
        if Succeedeq == 1:
            fi.write("k_f [1/s]")
            fi.write("\t")
            fi.write(str(kf))
            fi.write("\n")
            fi.write("k_-f [1/s]")
            fi.write("\t")
            fi.write(str(kmf))
            fi.write("\n")
            fi.write("Kf")
            fi.write("\t")
            fi.write(str(Kf))
            fi.write("\n")
            fi.write("f")
            fi.write("\t")
            fi.write(str(f))
            fi.write("\n")
        
        if DeltaXiModder == 1:
            fi.write("Xi-Resol. got modified")
            fi.write("\t")
            fi.write(str(ModDeltXi))
            fi.write("\n")
            
        
        if Statistical  == 0:
    
            if model_b == 1:
                fi.write("d in [cm]")
                fi.write("\t")
                fi.write(str(distance))
                fi.write("\n")

            if model_c == 1 or model_e == 1 or model_f == 1 or model_g == 1:
                fi.write("radius in [cm]")
                fi.write("\t")
                fi.write(str(r))
                fi.write("\n")

            if model_d   == 1:
                fi.write("radius in [cm]")
                fi.write("\t")
                fi.write(str(r))
                fi.write("\n")
                fi.write("d in [cm]")
                fi.write("\t")
                fi.write(str(distance))
                fi.write("\n")
        
        if Statistical  == 1:
            if DistFunc_Modifier == 1:
                fi.write("Dist. Func. was modified")
                fi.write("\n")
            if DistFunc_Modifier == 0:
                fi.write("Default Dist. Func. was used")
                fi.write("\n")
            
            fi.write("Dist. Func. param. a")
            fi.write("\t")
            fi.write(str(DFP_a))
            fi.write("\n")
            
            fi.write("Dist. Func. param. b")
            fi.write("\t")
            fi.write(str(DFP_b))
            fi.write("\n")
            
            fi.write("Dist. Func. param. c")
            fi.write("\t")
            fi.write(str(DFP_c))
            fi.write("\n")
            
            if model_b == 1:
                fi.write("sheets per mm")   
            if model_d == 1 or model_e == 1:
                fi.write("cylinders per mm^2")
            if model_g == 1:
                fi.write("spheres per mm^3")
            fi.write("\t")
            fi.write(str(Number_N))
            fi.write("\n")
                
            
            if model_d   == 1:
                fi.write("radius of Fiber in 10^-6[m]")
                fi.write("\t")
                fi.write(str(r*10000))
                fi.write("\n")
            if model_b   == 1:
                fi.write("sheet thickness 10^-6[m]")
                fi.write("\t")
                fi.write(str(sheet_thickness*10000))
                fi.write("\n")
                
                
        
            
        if FIT == 1:
            
            fi.write("E_zero vs. E_ref in [V]")
            fi.write("\t")
            fi.write(str(Ezero))
            fi.write("\n")
            fi.write("R_u in [Ohm]")
            fi.write("\t")
            fi.write(str(Ru))
            fi.write("\n")
            fi.write("Standard deviation")
            fi.write("\t")
            fi.write(str(Standardabweichung))
            fi.write("\n")
            if LinBaseCorr == 1:
                fi.write("Base was corrected")
                fi.write("\n")
                fi.write("Slope Start Pot vs. E_zero [V]")
                fi.write("\t")
                fi.write(str(StartPotLinCorr))
                fi.write("\n")
                fi.write("Slope End Pot vs. E_zero [V]")
                fi.write("\t")
                fi.write(str(EndPotLinCorr))
                fi.write("\n")
                fi.write("Slope of Baseline [micA/V]")
                fi.write("\t")
                fi.write(str(BaseSlope))
                fi.write("\n")
            
        
        for i in range(3):
            fi.write("\n")
        
        fi.write("CURVE DATA")
        
        
#---------------------------------------------------------------------------------------------------------------------------        
#Wenn der Fit null ist... also wenn nur berechnet wird...  

        if FIT == 0:  
            fi.write("\n")
            fi.write("E_vs_E_zero [V]")
            fi.write("\t")
            fi.write("I[mA]")
            fi.write("\t")
            if Chi_vs_Xi_Saver == 1:
                fi.write("Xi")
                fi.write("\t")
                fi.write("Chi")
                fi.write("\t")
            if Statistical  == 0:
                if Convolution_Saver == 1:
                    fi.write("t")
                    fi.write("\t")
                    fi.write("Convolution_integral_forward")
                    fi.write("\t")
                    fi.write("Convolution_integral_backward")
                    fi.write("\t")              
            fi.write("\n")

        if FIT == 0:
            for i in range(len(xxx)-1):
                fi.write(str(np.asscalar(PotCalc[i])))
                fi.write("\t")
                if Statistical  == 0:
                    fi.write(str(np.asscalar(CurrCalc[i])))
                    fi.write("\t")
                if Statistical  == 1:
                    fi.write(str(np.asscalar(Chi_ges_Array[i]*(n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000)))
                    fi.write("\t")
                if Statistical  == 0:
                    if Chi_vs_Xi_Saver == 1:
                        fi.write(str(np.asscalar(xxx[i])))
                        fi.write("\t")
                        fi.write(str(np.asscalar(yyy[i])))
                        fi.write("\t")
                    if Convolution_Saver == 1:
                        fi.write(str(np.asscalar(t_Array_cont[i])))
                        fi.write("\t")
                        fi.write(str(np.asscalar(WuFuInt_p[i])))
                        fi.write("\t")
                        fi.write(str(np.asscalar(WuFuInt_f[i])))
                        fi.write("\t")
                if Statistical  == 1:
                    if Chi_vs_Xi_Saver == 1:
                        fi.write(str(np.asscalar(xxx[i])))
                        fi.write("\t")
                        fi.write(str(np.asscalar(Chi_ges_Array[i])))
                        fi.write("\t")
                fi.write("\n")
            fi.write("\n")
            fi.write("\n")
            
        if FIT == 0:  
            if Statistical  == 1:
                if Save_Individuals == 1:
                    fi.write("Mean center to wall in 10^-6[m]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write(str(np.asscalar(0.5*RR_Mean_Array[i])))                        
                        fi.write("\t")
                    fi.write("\n")
                    
                    fi.write("Center to wall Integ-Limits in 10^-6[m]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write(str(np.asscalar(0.5*R_Mean_Array[i])))                        
                        fi.write("\t")
                    fi.write("\n")
                    
                    
                    fi.write("\n")
                    fi.write("E_vs_E_zero [V]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write("weighted I[mA]")
                        fi.write("\t")
                    fi.write("\n")
                    for i in range(len(PotCalc)-1):
                        fi.write(str(np.asscalar(PotCalc[i])))
                        fi.write("\t")
                        for j in range(len(R_Mean_Array)):
                            fi.write(str(np.asscalar(CurrCalc[i,j])))
                            fi.write("\t")       
                        fi.write("\n")
                    fi.write("\n")
                    fi.write("\n")
                
                    if Chi_vs_Xi_Saver == 1: 
                        fi.write("Mean center to wall in 10^-6[m]")
                        fi.write("\t")
                        for i in range(len(R_Mean_Array)):
                            fi.write(str(np.asscalar(0.5*RR_Mean_Array[i])))                        
                            fi.write("\t")
                        fi.write("\n")
                        
                        fi.write("Center to wall Integ-Limits in 10^-6[m]")
                        fi.write("\t")
                        for i in range(len(R_Mean_Array)):
                            fi.write(str(np.asscalar(0.5*R_Mean_Array[i])))                        
                            fi.write("\t")
                        fi.write("\n")
                    
                        fi.write("\n")
                        fi.write("Xi")
                        fi.write("\t")
                        for i in range(len(R_Mean_Array)):
                            fi.write("weighted Chi")
                            fi.write("\t")
                        fi.write("\n")
                        for i in range(len(PotCalc)-1):
                            fi.write(str(np.asscalar(xxx[i])))
                            fi.write("\t")
                            for j in range(len(R_Mean_Array)):
                                fi.write(str(np.asscalar(CurrCalc[i,j]/((n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5)*1000))))
                                fi.write("\t")       
                            fi.write("\n")
                        fi.write("\n")
                        fi.write("\n")
                                 
        if FIT == 0:  
            if Statistical  == 1:
                if Save_Dist_Func == 1:
                    fi.write("center to wall / 10^-6 m")
                    fi.write("\t")
                    fi.write("P(center to wall) /(10^-6 m)^-1")
                    fi.write("\t")
                    fi.write("\n")
                    for i in range(len(R_Array)):
                        #weil es ja eigentlich Abstände von Fasermittelpunkten sind!!! durch 2
                        fi.write(str(np.asscalar(0.5*R_Array[i])))
                        fi.write("\t")
                        #Mal 2 damit Normierung stimmt!
                        fi.write(str(np.asscalar(2*P_Array[i])))
                        fi.write("\t")
                        fi.write("\n")
                    fi.write("\n")
                    fi.write("\n")
                    fi.write("\n")                       

#---------------------------------------------------------------------------------------------------------------------------        
#Wenn der Fit gleich eins ist, also gefittet wird
#---------------------------------------------------------------------------------------------------------------------------        

        if FIT == 1:
            fi.write("\n")
            fi.write("Calculated Potential vs. E_zero in [V]")
            fi.write("\t")
            fi.write("Calculated Current in microampere")
            fi.write("\t")
            fi.write("\n")
            
            if Statistical  == 0:
                for i in range(len(xxx)-1):
                    fi.write(str(np.asscalar(xxx[i])))
                    fi.write("\t")
                    fi.write(str(np.asscalar(yyy[i])))
                    fi.write("\t")
                    fi.write("\n")
                fi.write("\n")
                fi.write("\n")
                fi.write("\n")
                
            if Statistical  == 1:
                for i in range(len(xxx)-1):
                    fi.write(str(np.asscalar(PotCalc[i])))
                    fi.write("\t")
                    fi.write(str(np.asscalar(Chi_ges_Array[i])))
                    fi.write("\t")
                    fi.write("\n")
                fi.write("\n")
                fi.write("\n")
                fi.write("\n")
                
                
            if Statistical  == 1:
                if Save_Individuals == 1:
                    fi.write("r to wall in 10^-6[m]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write(str(np.asscalar(0.5*RR_Mean_Array[i])))                        
                        fi.write("\t")
                    fi.write("\n")
                    
                    fi.write("Center to wall Integ-Limits in 10^-6[m]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write(str(np.asscalar(0.5*R_Mean_Array[i])))                        
                        fi.write("\t")
                    fi.write("\n")
                    
                    fi.write("\n")
                    fi.write("E_vs_E_zero [V]")
                    fi.write("\t")
                    for i in range(len(R_Mean_Array)):
                        fi.write("weighted I[mA]")
                        fi.write("\t")
                    fi.write("\n")
                    for i in range(len(PotCalc)-1):
                        fi.write(str(np.asscalar(PotCalc[i])))
                        fi.write("\t")
                        for j in range(len(R_Mean_Array)): #soll in milliamps sein, deshalb mal 1000
                            fi.write(str(np.asscalar(1000*ChiSuperarray[i,j]*(n*F*A*c*(D_f*n*F*Scanrate/(R*T))**0.5))))
                            fi.write("\t")       
                        fi.write("\n")
                    fi.write("\n")
                    fi.write("\n")
                
                if Save_Dist_Func == 1:
                    fi.write("center to wall / 10^-6 m")
                    fi.write("\t")
                    fi.write("P(center to wall) /(10^-6 m)^-1")
                    fi.write("\t")
                    fi.write("\n")
                    for i in range(len(R_Array)):
                        #weil es ja eigentlich Abstände von Fasermittelpunkten sind!!! durch 2
                        fi.write(str(np.asscalar(0.5*R_Array[i])))
                        fi.write("\t")
                        #Mal 2 damit Normierung stimmt!
                        fi.write(str(np.asscalar(2*P_Array[i])))
                        fi.write("\t")
                        fi.write("\n")
                    
                fi.write("\n")
                fi.write("\n")
                fi.write("\n")


#---------------------------------------------------------------------------------------------------------------------------        
     
            fi.write("R_Corrected measured Potential vs. E_zero in [V]")
            fi.write("\t")
            fi.write("Measured Current in microampere")
            fi.write("\t")
            fi.write("\n")
            
            for i in range(len(Pot)):
                fi.write(str(np.asscalar(Pot[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Curr[i])))
                fi.write("\t")
                fi.write("\n")            
                        
                
    root.destroy()




# In[48]:


#---------------------------------------------------------------------------------------
#IMPEDANCE-STUFF
#---------------------------------------------------------------------------------------

def Semi_Inf_Plan_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 1
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()
def Fin_Plan_Trans_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 1
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()
def Fin_Plan_Ref_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 1
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()
def Semi_Inf_Zyl_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 1
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()
def Semi_Inf_Sphe_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 1
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()
def Finit_Int_Sphe_Imp_Sim():
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 1
    global IMPFIT
    global IMPSIM
    IMPFIT  = 0
    IMPSIM  = 1
    Imp_Entry_Window()


#-----------------------------------------------------

def Semi_Inf_Plan_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 1
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()
def Fin_Plan_Trans_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 1
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()
def Fin_Plan_Ref_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 1
    Imp4 = 0
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()
def Semi_Inf_Zyl_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 1
    Imp5 = 0
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()
def Semi_Inf_Sphe_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 1
    Imp6 = 0
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()
def Finit_Int_Sphe_Imp_Fit():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    global Imp1
    global Imp2
    global Imp3
    global Imp4
    global Imp5
    global Imp6
    Imp1 = 0
    Imp2 = 0
    Imp3 = 0
    Imp4 = 0
    Imp5 = 0
    Imp6 = 1
    global IMPFIT
    global IMPSIM
    IMPFIT  = 1
    IMPSIM  = 0
    if File_Was_Loaded ==1:
        Imp_Entry_Window()




# In[49]:


def Imp_Entry_Window():
    
    Fenster = Toplevel()  
    if IMPFIT == 0:
        Fenster.geometry("400x600")
    if IMPFIT == 1:
        Fenster.geometry("400x600")
        
    if Imp1 == 1:
        Fenster.title("Semi-Infinite Planar Diffusion Impedance")
    if Imp2 == 1:
        Fenster.title("Finite Planar Transmissive Diffusion Impedance")
    if Imp3 == 1:
        Fenster.title("Finite Planar Refelctive Diffusion Impedance")
    if Imp4 == 1:
        Fenster.title("Semi-Infinite Cylindrical Diffusion Impedance")
    if Imp5 == 1:
        Fenster.title("Semi-Infinite Spherical Diffusion Impedance")
    if Imp6 == 1:
        Fenster.title("Finite Spherical internal Diffusion Impedance")
        
  
    
    n_Label = Label(Fenster, text="n*")
    n_Label.grid(row=0, column=0)
    n_Eingabe = Entry(Fenster)                                               
    n_Eingabe.grid(row=0, column=1)
        
    
    alpha_Label = Label(Fenster, text="alpha*")
    alpha_Label.grid(row=1, column=0)
    alpha_Eingabe = Entry(Fenster)                                               
    alpha_Eingabe.grid(row=1, column=1)
    
    
    kzero_Label = Label(Fenster,text="k_zero [cm/s]*")
    kzero_Label.grid(row=2, column=0)
    kzero_Eingabe = Entry(Fenster)
    kzero_Eingabe.grid(row=2, column=1)
    
    
    Ezero_Label = Label(Fenster,text="E_zero vs Ref [V]*")
    Ezero_Label.grid(row=3, column=0)
    Ezero_Eingabe = Entry(Fenster)
    Ezero_Eingabe.grid(row=3, column=1)
    
    Eeq_Label = Label(Fenster,text="E_eq vs Ref [V]*")
    Eeq_Label.grid(row=4, column=0)
    Eeq_Eingabe = Entry(Fenster)
    Eeq_Eingabe.grid(row=4, column=1)
    
    T_Label = Label(Fenster,text="T [°C]*")
    T_Label.grid(row=5, column=0)
    T_Eingabe = Entry(Fenster)
    T_Eingabe.grid(row=5, column=1)
    
    if IMPFIT == 0:
    
        Freq_High_Label = Label(Fenster,text="Highest Frequency [Hz]*")
        Freq_High_Label.grid(row=6, column=0)
        Freq_High_Eingabe = Entry(Fenster)
        Freq_High_Eingabe.grid(row=6, column=1)

        Freq_Low_Label = Label(Fenster,text="Lowest Frequency [Hz]*")
        Freq_Low_Label.grid(row=7, column=0)
        Freq_Low_Eingabe = Entry(Fenster)
        Freq_Low_Eingabe.grid(row=7, column=1)

        Freq_Num_Label = Label(Fenster,text="Number of Freq's*")
        Freq_Num_Label.grid(row=8, column=0)
        Freq_Num_Eingabe = Entry(Fenster)
        Freq_Num_Eingabe.grid(row=8, column=1)

    
    if Imp2 == 1:
        d_Label = Label(Fenster,text="distance 10^-6[m]*")
        d_Label.grid(row=9, column=0)
        d_Eingabe = Entry(Fenster)
        d_Eingabe.grid(row=9, column=1)
    
    if Imp3 == 1:
        d_Label = Label(Fenster,text="distance 10^-6[m]*")
        d_Label.grid(row=9, column=0)
        d_Eingabe = Entry(Fenster)
        d_Eingabe.grid(row=9, column=1)
    
    if Imp4 == 1:
        r_Label = Label(Fenster,text="radius 10^-6[m]*")
        r_Label.grid(row=9, column=0)
        r_Eingabe = Entry(Fenster)
        r_Eingabe.grid(row=9, column=1)
    
    if Imp5 == 1:
        r_Label = Label(Fenster,text="radius 10^-6[m]*")
        r_Label.grid(row=9, column=0)
        r_Eingabe = Entry(Fenster)
        r_Eingabe.grid(row=9, column=1)
        
    if Imp6 == 1:
        r_Label = Label(Fenster,text="radius 10^-6[m]*")
        r_Label.grid(row=9, column=0)
        r_Eingabe = Entry(Fenster)
        r_Eingabe.grid(row=9, column=1)
        
    
   
    Ru_Label = Label(Fenster,text="Ru in Ohm*")
    Ru_Label.grid(row=10, column=0)
    Ru_Eingabe = Entry(Fenster)
    Ru_Eingabe.grid(row=10, column=1)
    
    Capacity_Label = Label(Fenster,text="C in Microfarad*")
    Capacity_Label.grid(row=11, column=0)
    Capacity_Eingabe = Entry(Fenster)
    Capacity_Eingabe.grid(row=11, column=1)
    
    A_Label = Label(Fenster,text="A in cm^2*")
    A_Label.grid(row=12, column=0)
    A_Eingabe = Entry(Fenster)
    A_Eingabe.grid(row=12, column=1)
        
    c_Label = Label(Fenster,text="c_total in mol/L*")
    c_Label.grid(row=13, column=0)
    c_Eingabe = Entry(Fenster)
    c_Eingabe.grid(row=13, column=1)
        
    D_Label = Label(Fenster,text="D 10^-6[cm^2/s]*")
    D_Label.grid(row=14, column=0)
    D_Eingabe = Entry(Fenster)
    D_Eingabe.grid(row=14, column=1)
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Unequal D", variable=var1).grid(row=15, column=0, sticky=W)
    
    Dred_Label = Label(Fenster,text="D_red 10^-6[cm^2/s]")
    Dred_Label.grid(row=16, column=0)
    Dred_Eingabe = Entry(Fenster)
    Dred_Eingabe.grid(row=16, column=1)
    
    Dox_Label = Label(Fenster,text="D_ox 10^-6[cm^2/s]")
    Dox_Label.grid(row=17, column=0)
    Dox_Eingabe = Entry(Fenster)
    Dox_Eingabe.grid(row=17, column=1)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Preceeding ch eq", variable=var2).grid(row=18, column=0, sticky=W)
    
    k_prec_f_Label = Label(Fenster,text="kf_chem_preceding [1/s]")
    k_prec_f_Label.grid(row=19, column=0)
    k_prec_f_Eingabe = Entry(Fenster)
    k_prec_f_Eingabe.grid(row=19, column=1)
    
    k_prec_b_Label = Label(Fenster,text="kb_chem_preceding [1/s]")
    k_prec_b_Label.grid(row=20, column=0)
    k_prec_b_Eingabe = Entry(Fenster)
    k_prec_b_Eingabe.grid(row=20, column=1)
    
    var3 = IntVar()
    Checkbutton(Fenster, text="Succeding ch eq", variable=var3).grid(row=21, column=0, sticky=W)
    
    k_succ_f_Label = Label(Fenster,text="kf_chem_succeding [1/s]")
    k_succ_f_Label.grid(row=22, column=0)
    k_succ_f_Eingabe = Entry(Fenster)
    k_succ_f_Eingabe.grid(row=22, column=1)
    
    k_succ_b_Label = Label(Fenster,text="kb_chem_succeding [1/s]")
    k_succ_b_Label.grid(row=23, column=0)
    k_succ_b_Eingabe = Entry(Fenster)
    k_succ_b_Eingabe.grid(row=23, column=1)
    
    
    var4 = IntVar()
    Checkbutton(Fenster, text="Save as txt", variable=var4).grid(row=24, column=0, sticky=W)
    
    
    
    #-------------------------------------------------------
    #Akzeptierfunktion schreiben
    #-------------------------------------------------------
    
    def AcceptParams():
        global n
        global alpha
        global kzero
        global Ezero
        global Eeq
        global T
        global distance
        global r
        global precedeq
        global succedeq
        global Unequal_D
        global AstxtSaver
        global A
        global c
        global Ru
        global FreqUp
        global FreqLow
        global FreqNum
        global Cap
    
        
    #------------------------------------------------------------
    #Felder lesen
    #------------------------------------------------------------
        
        n           = (float(n_Eingabe.get()))
        alpha       = (float(alpha_Eingabe.get()))
        kzero       = (float(kzero_Eingabe.get()))
        T           = (float(T_Eingabe.get())) + 273.15
        Cap         = (float(Capacity_Eingabe.get()))*0.000001
        
        if IMPFIT == 0:
            FreqUp      = np.log10(float(Freq_High_Eingabe.get()))
            FreqLow     = np.log10(float(Freq_Low_Eingabe.get()))
            FreqNum     = (float(Freq_Num_Eingabe.get()))
        
    
        if Imp2 == 1:
            distance     = (float(d_Eingabe.get()))*0.0001
        if Imp3 == 1:
            distance     = (float(d_Eingabe.get()))*0.0001
        if Imp4 == 1:
            r            = (float(r_Eingabe.get()))*0.0001
        if Imp5 == 1:
            r            = (float(r_Eingabe.get()))*0.0001
        if Imp6 == 1:
            r            = (float(r_Eingabe.get()))*0.0001
    
         
        global Dox
        global Dred
        
        Unequal_D   = var1.get()
        if Unequal_D ==0:
            Dox       = (float(D_Eingabe.get()))*0.000001
            Dred      = (float(D_Eingabe.get()))*0.000001
        
        if Unequal_D ==1:
            Dox      = (float(Dox_Eingabe.get()))*0.000001
            Dred     = (float(Dred_Eingabe.get()))*0.000001
        
    
        Ru     = (float(Ru_Eingabe.get()))
        c      = (float(c_Eingabe.get()))*0.001
        A      = (float(A_Eingabe.get()))
        Ezero  = (float(Ezero_Eingabe.get()))
        Eeq    = (float(Eeq_Eingabe.get()))
        
        precedeq = var2.get()
        if precedeq == 1:
            global preced_kf
            global preced_kb
            preced_kf = (float(k_prec_f_Eingabe.get()))
            preced_kb = (float(k_prec_b_Eingabe.get()))
            
        succedeq = var3.get()
        if succedeq == 1:
            global succed_kf
            global succed_kb
            succed_kf = (float(k_succ_f_Eingabe.get()))  
            succed_kb = (float(k_succ_b_Eingabe.get())) 
            
        AstxtSaver  = var4.get()
        
        
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=25,column=1)
    button=Button(Fenster,text="Next",command=Imp_Type_Chooser).grid(row=26,column=1)
    
    
    


# In[50]:


def Imp_Type_Chooser():
    if Imp1 == 1:
        if IMPFIT == 0:
            PlanarSemiinfImpCalculator()
        if IMPFIT == 1:
            ImpStartEndDefiner()
            PlanarSemiinfImpCalculator()
            
    if Imp2 == 1:
        if IMPFIT == 0:
            PlanarTransmImpCalculator()
        if IMPFIT == 1:
            ImpStartEndDefiner()
            PlanarTransmImpCalculator()
        
    if Imp3 == 1:
        if IMPFIT == 0:
            PlanarReflImpCalculator()
        if IMPFIT == 1:
            ImpStartEndDefiner()
            PlanarReflImpCalculator()
        
    if Imp4 == 1:
        if IMPFIT == 0:
            ZylSemiinfImpCalculator()        
        if IMPFIT == 1:
            ImpStartEndDefiner()
            ZylSemiinfImpCalculator()
        
    if Imp5 == 1:
        if IMPFIT == 0:
            SpherSemiinfImpCalculator()
        if IMPFIT == 1:
            ImpStartEndDefiner()
            SpherSemiinfImpCalculator()
        
    if Imp6 == 1:
        if IMPFIT == 0:
            SpherIntFinImpCalculator()
        if IMPFIT == 1:
            ImpStartEndDefiner()
            SpherIntFinImpCalculator()
    
    


# In[51]:


def ImpStartEndDefiner():
    global Z_im_Meas
    global Z_re_Meas
    Z_im_Meas = np.squeeze(Strom)        #Nach dem Laden heißen die Arrays nur so wie früher
    Z_re_Meas = np.squeeze(Potenzial)    #deshalb hier die Umbenennung
    


# In[52]:


def PlanarSemiinfImpCalculator():
    
    root = Toplevel()
    root.title("Simulated Nyquist-Plot")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]

    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(-alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp((1-alpha)*n*F*(Eeq-EZero)/(R*T))


    ZFar = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +kf_el/((1+Kp)*((Freq+lp)*Dred)**0.5)   +kf_el*Kp/((1+Kp)*((Freq)*Dred)**0.5)       +kb_el*Ks/((1+Ks)*((Freq+ls)*Dox)**0.5)   +kb_el/((1+Ks)*((Freq)*Dox)**0.5)      )                              

    #----------------------------------------------------------------------------------------------
   
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
        #==============================================
        #calculate Standard deviation
        #==============================================
        #umdrehen zum Interpolieren
        #Z_imag_Calc = Z_imag_Calc[::-1]
        #Z_real_Calc = Z_real_Calc[::-1]
        #global ReStDev
        #global ImStDev
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
        #Zurück umdrehen nach Interpolieren
        #==============================================
        #Z_imag_Calc = Z_imag_Calc[::-1]
        #Z_real_Calc = Z_real_Calc[::-1]
        
        #==============================================
        
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm  ', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    
    


# In[53]:


def PlanarTransmImpCalculator():
    root = Toplevel()
    root.title("Simulated Nyquist-Plot")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]
    
    
    
    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    global d
    
    d       = distance
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(-alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp((1-alpha)*n*F*(Eeq-EZero)/(R*T))


    ZFar = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +kf_el*np.tanh(d*((Freq+lp)/Dred)**0.5)/((1+Kp)*((Freq+lp)*Dred)**0.5)   +kf_el*Kp*np.tanh(d*(Freq/Dred)**0.5)/((1+Kp)*((Freq)*Dred)**0.5)       +kb_el*Ks*np.tanh(d*((Freq+ls)/Dox)**0.5)/((1+Ks)*((Freq+ls)*Dox)**0.5)   +kb_el*np.tanh(d*(Freq/Dox)**0.5)/((1+Ks)*((Freq)*Dox)**0.5)      )                              
    
    #----------------------------------------------------------------------------------------------
   
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
        #==============================================
        #calculate Standard deviation
        #==============================================
        #umdrehen zum Interpolieren
        #Z_imag_Calc = Z_imag_Calc[::-1]
        #Z_real_Calc = Z_real_Calc[::-1]
        #global ReStDev
        #global ImStDev
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
        #Zurück umdrehen nach Interpolieren
        #==============================================
        #Z_imag_Calc = Z_imag_Calc[::-1]
        #Z_real_Calc = Z_real_Calc[::-1]
        
        #==============================================
        
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    
    


# In[54]:


def PlanarReflImpCalculator():
    root = Toplevel()
    root.title("Simulated Nyquist-Plot")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]
    
    
    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    global d
    
    d       = distance
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(-alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp((1-alpha)*n*F*(Eeq-EZero)/(R*T))


    def coth(x):
        return (np.tanh(x))**(-1)
    ZFar = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +kf_el*coth(d*((Freq+lp)/Dred)**0.5)/((1+Kp)*((Freq+lp)*Dred)**0.5)   +kf_el*Kp*coth(d*(Freq/Dred)**0.5)/((1+Kp)*((Freq)*Dred)**0.5)       +kb_el*Ks*coth(d*((Freq+ls)/Dox)**0.5)/((1+Ks)*((Freq+ls)*Dox)**0.5)   +kb_el*coth(d*(Freq/Dox)**0.5)/((1+Ks)*((Freq)*Dox)**0.5)      )                              

    #----------------------------------------------------------------------------------------------
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
        #==============================================
        #calculate Standard deviation
        #==============================================
        
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
          
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    
    


# In[55]:


def ZylSemiinfImpCalculator():
    root = Toplevel()
    root.title("Simulated Nyquist-Plot")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]
    
    
    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(-alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp((1-alpha)*n*F*(Eeq-EZero)/(R*T))


    
    def Bessel_P_eq(x):
        return (kv(0,(r*((x+lp)/Dred)**0.5)))/(kv(1,(r*((x+lp)/Dred)**0.5))*((x+lp)*Dred)**0.5)
    def Bessel_P_n(x):
        return (kv(0,(r*((x)/Dred)**0.5)))/(kv(1,(r*((x)/Dred)**0.5))*((x)*Dred)**0.5)
    def Bessel_S_eq(x):
        return (kv(0,(r*((x+ls)/Dox)**0.5)))/(kv(1,(r*((x+ls)/Dox)**0.5))*((x+ls)*Dox)**0.5)
    def Bessel_S_n(x):
        return (kv(0,(r*((x)/Dox)**0.5)))/(kv(1,(r*((x)/Dox)**0.5))*((x)*Dox)**0.5)


    ZFar = np.empty(len(Freq),dtype = complex)
    for i in range(len(ZFar)):
        if np.abs(kv(1,(r*(Freq[i]/Dox)**0.5)).imag) != 0 :
            if np.abs(kv(1,(r*(Freq[i]/Dred)**0.5)).imag) != 0 :
                ZFar[i] = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +(kf_el/(1+Kp))*Bessel_P_eq(Freq[i])   +(kf_el*Kp/(1+Kp))*Bessel_P_n(Freq[i])  +(kb_el*Ks/(1+Ks))*Bessel_S_eq(Freq[i])   +(kb_el/(1+Ks))*Bessel_S_n(Freq[i]) )                             
        if np.abs(kv(1,(r*(Freq[i]/Dox)**0.5))) == 0: 
            ZFar[i] = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +(kf_el/((1+Kp)*((Freq[i]+lp)*Dred)**0.5)) +(kf_el*Kp/((1+Kp)*((Freq[i])*Dred)**0.5))  +(kb_el*Ks/((1+Ks)*((Freq[i]+lp)*Dox)**0.5))   +(kb_el/((1+Ks)*((Freq[i])*Dox)**0.5)) )                             
        if np.abs(kv(1,(r*(Freq[i]/Dred)**0.5))) == 0:  
            ZFar[i] = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +(kf_el/((1+Kp)*((Freq[i]+lp)*Dred)**0.5)) +(kf_el*Kp/((1+Kp)*((Freq[i])*Dred)**0.5))  +(kb_el*Ks/((1+Ks)*((Freq[i]+lp)*Dox)**0.5))   +(kb_el/((1+Ks)*((Freq[i])*Dox)**0.5)) )                             
    
    #----------------------------------------------------------------------------------------------
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
        
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
       
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    
    


# In[56]:


def SpherSemiinfImpCalculator():
    root = Toplevel()
    root.title("Simulated Nyquist-Plot")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]
    
    
    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(-alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp((1-alpha)*n*F*(Eeq-EZero)/(R*T))


    ZFar = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +kf_el/((1+Kp)*(Dred/r + ((Freq+lp)*Dred)**0.5))   +kf_el*Kp/((1+Kp)*(Dred/r + ((Freq)*Dred)**0.5))       +kb_el*Ks/((1+Ks)*((Dox/r + ((Freq+ls)*Dox)**0.5) ))   +kb_el/((1+Ks)*((Dox/r + ((Freq)*Dox)**0.5) ))      )                              

    #----------------------------------------------------------------------------------------------
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
        
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
        
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    
    
   


# In[57]:


def SpherIntFinImpCalculator():
    root = Toplevel()
    root.title("Simulated Nyquist-Plot")
    
    #-------------------------------------------------------------------------------------------------- 
    #Jetzt Zuerst Parameter holen und setzen
    #--------------------------------------------------------------------------------------------------       
    #--------------------------------------------------------------------------------
    #Warburg Planar Semiinf Vorgelagertes UND Nachgelagertes homogenes chemisches GG
    #--------------------------------------------------------------------------------
    
    global Freq
    if IMPFIT == 0:
        Freq    = (1j)*2*np.pi*np.logspace(FreqLow,FreqUp,FreqNum)
    if IMPFIT == 1:
        global Z_imag_Meas
        global Z_real_Meas
        Freq        = (1j)*2*np.pi*np.squeeze(Frequency_Array[::-1])
        Z_imag_Meas = Z_im_Meas[::-1]
        Z_real_Meas = Z_re_Meas[::-1]
        
        if Freq[0].imag > Freq[1].imag:
            Freq        = Freq[::-1]
            Z_imag_Meas = Z_imag_Meas[::-1]
            Z_real_Meas = Z_real_Meas[::-1]
            
            
    global Kp
    global lp
    global Ks
    global ls
    global kfp
    global kbp
    global kfs
    global kbs
    c_tot   = c
    R_u     = Ru 
    EZero   = Ezero
    
    kfp     = 10000.0
    kbp     = 0.0001
    Kp      = kfp/kbp
    lp       = kfp+kbp
    
    kfs     = 0.0001
    kbs     = 10000.0
    Ks      = kfs/kbs
    ls       = kfs+kbs
    
    if precedeq == 1:
        #preceding
        kfp     = preced_kf
        kbp     = preced_kb
        Kp      = kfp/kbp
        lp       = kfp+kbp
    if succedeq == 1:
        #succeding
        kfs     = succed_kf
        kbs     = succed_kb
        Ks      = kfs/kbs
        ls       = kfs+kbs
    
    #other stuff
    global cox
    global cred
    global kf_ef
    global kb_el
    cox     = c_tot /((1 + 1/Kp)*np.exp(-(n*F/(R*T))*(Eeq-EZero))+(1+Ks))        #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    cred    = c_tot /((1 + 1/Kp) + (1+Ks)*np.exp((n*F/(R*T))*(Eeq-EZero)))       #Bewirkt, dass Immer die Gleichgewkonz genopmmen
    kf_el   = kzero*np.exp(-alpha*n*F*(Eeq-EZero)/(R*T))
    kb_el   = kzero*np.exp((1-alpha)*n*F*(Eeq-EZero)/(R*T))


    def coth(x):
        return (np.tanh(x))**(-1)

    ZFar = (R*T/(A*(n*F)**2))*(1/(alpha*kf_el*cred+(1-alpha)*kb_el*cox)) * (1  +kf_el/(coth(r*((Freq+lp)/Dred)**0.5)*(1+Kp)*(Dred/r + ((Freq+lp)*Dred)**0.5))   +kf_el*Kp/(coth(r*(Freq/Dred)**0.5)*(1+Kp)*(Dred/r + ((Freq)*Dred)**0.5))       +kb_el*Ks/(coth(r*((Freq+lp)/Dox)**0.5)*(1+Ks)*((Dox/r + ((Freq+ls)*Dox)**0.5) ))   +kb_el/(coth(r*(Freq/Dox)**0.5)*(1+Ks)*((Dox/r + ((Freq)*Dox)**0.5) ))      )                              

    #----------------------------------------------------------------------------------------------
    Z  = R_u + 1/(1.0/ZFar + Freq*Cap)    
    
    
    global xxx
    global f
    
    
    xxx = Z.real
    f   = Z.imag
    
        
    #-------------------------------------------------------------------------------------------------------------------------
    #PLOTTEN VON SIMULIERTEM
    #-------------------------------------------------------------------------------------------------------------------------
    
    
        
    Abbildung = Figure(figsize=(10, 4), dpi=100)
    b1 = Abbildung.add_subplot(121)
    b2 = Abbildung.add_subplot(122)
    
    
       
    if IMPFIT == 0:
        b1.plot(xxx[::],-f[::],color='k', linestyle ='-', marker ='.')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -f[::], linestyle='-',marker='.',color='b', label = '-Z.imag')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), xxx[::], linestyle='-',marker='.',color='r',label = 'Z.real')
    if IMPFIT == 1:
        global Z_imag_Calc
        Z_imag_Calc = f
        global Z_real_Calc
        Z_real_Calc = xxx
        
    
        ReStDev_Array = np.zeros(len(f))
        ImStDev_Array = np.zeros(len(f))
        for i in range(len(f)):
            ReStDev_Array[i] = ((Z_real_Calc[i]- Z_real_Meas[i])/np.max(Z_real_Calc))**2
            ImStDev_Array[i] = ((Z_imag_Calc[i] - Z_imag_Meas[i])/np.min(Z_imag_Calc))**2
        ReStDev = (np.sum(ReStDev_Array)/float(len(ReStDev_Array)))**0.5  
        ImStDev = (np.sum(ImStDev_Array)/float(len(ImStDev_Array)))**0.5  
        
        
        b1.plot(Z_real_Calc,-Z_imag_Calc,color='k', linestyle ='-', marker ='.')
        b1.plot(Z_real_Meas,-Z_imag_Meas,color='k', linestyle ='-', marker ='')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Calc, linestyle='-',marker='.',color='b', label = '-Z.imag_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), -Z_imag_Meas, linestyle='-',marker='',color='b', label = '-Z.imag_m')
        
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Calc, linestyle='-',marker='.',color='r',label = 'Z.real_c')
        b2.plot(np.log10(Freq.imag/(2*np.pi)), Z_real_Meas, linestyle='-',marker='',color='r',label = 'Z.real_m')
        
        b2.annotate(u'\u03C3''_re' '=%8.3f' % ReStDev, xy=(0.05, 0.93), xycoords='axes fraction')
        b2.annotate(u'\u03C3''_im' '=%8.3f' % ImStDev, xy=(0.05, 0.83), xycoords='axes fraction')
        
        
    b1.set_title("Nyquist-Plot")    
    b1.set_xlabel('Z real / Ohm', fontsize=12)
    b1.set_ylabel('-Z imag / Ohm', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b1.spines[axis].set_linewidth(2)
        b1.spines[axis].set_color('k')
    b1.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    
    
    b2.set_title("Bode-Plot")
    for axis in ['top','bottom','left','right']:
        b2.spines[axis].set_linewidth(2)
        b2.spines[axis].set_color('k')
    b2.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    b2.set_xlabel('log10(freq / s^-1)', fontsize=12)
    b2.set_ylabel('-Z_imag and Z_real / Ohm', fontsize=12)
    b2.plot()
    b2.legend()
    
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    
    if AstxtSaver == 1:
        Simulated_IMPEDANZ_as_txt_saver()
   
    
    
   


# In[58]:


def Simulated_IMPEDANZ_as_txt_saver():
    
    root = Toplevel() 

    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))


    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")
        
        
        if Imp1 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Planar Semi-infinit")
            fi.write("\n")
        if Imp2 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Planar Finit Transmissive")
            fi.write("\n")
        if Imp3 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Planar Finit Reflective")
            fi.write("\n")
        if Imp4 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Cylindrical semi infinte")
            fi.write("\n")
        if Imp5 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Spherical Semi-infinite")
            fi.write("\n")
        if Imp6 == 1:
            fi.write("Model:")
            fi.write("\t")
            fi.write("Spherical Finite Internal")
            fi.write("\n")

        fi.write("n")
        fi.write("\t")
        fi.write(str(n))
        fi.write("\n")
        fi.write("T in [K]")
        fi.write("\t")
        fi.write(str(T))
        fi.write("\n")
        
        fi.write("D_red in [cm^2/s]")
        fi.write("\t")
        fi.write(str(Dred))
        fi.write("\n")
        fi.write("Dox in [cm^2/s]")
        fi.write("\t")
        fi.write(str(Dox))
        fi.write("\n")
        if precedeq == 1:
            #preceding
            fi.write("k_f_preceeding [1/s]")
            fi.write("\t")
            fi.write(str(kfp))
            fi.write("\n")
            fi.write("k_b_preceeding [1/s]")
            fi.write("\t")
            fi.write(str(kbp))
            fi.write("\n")
        if succedeq == 1:
            #succeding
            fi.write("k_f_succeeding [1/s]")
            fi.write("\t")
            fi.write(str(kfs))
            fi.write("\n")
            fi.write("k_b_succeeding [1/s]]")
            fi.write("\t")
            fi.write(str(kbs))
            fi.write("\n")      
        
        
        fi.write("alpha")
        fi.write("\t")
        fi.write(str(alpha))
        fi.write("\n")
        fi.write("k_zero in [cm/s]")
        fi.write("\t")
        fi.write(str(kzero))
        fi.write("\n")
            

        if Imp2 == 1:
            fi.write("d in [cm]")
            fi.write("\t")
            fi.write(str(distance))
            fi.write("\n")
        if Imp3 == 1:
            fi.write("d in [cm]")
            fi.write("\t")
            fi.write(str(distance))
            fi.write("\n")
        if Imp4 == 1:
            fi.write("radius in [cm]")
            fi.write("\t")
            fi.write(str(r))
            fi.write("\n")
        if Imp5 == 1:
            fi.write("radius in [cm]")
            fi.write("\t")
            fi.write(str(r))
            fi.write("\n")
        if Imp6 == 1:
            fi.write("radius in [cm]")
            fi.write("\t")
            fi.write(str(r))
            fi.write("\n")

        fi.write("E_zero vs. E_ref in [V]")
        fi.write("\t")
        fi.write(str(Ezero))
        fi.write("\n")
        fi.write("Eeq vs. E_ref in [V]")
        fi.write("\t")
        fi.write(str(Eeq))
        fi.write("\n")
        fi.write("Ru in [Ohm]")        
        fi.write("\t")
        fi.write(str(Ru))
        fi.write("\n")
        fi.write("A in [cm^2]")
        fi.write("\t")
        fi.write(str(A))
        fi.write("\n")
        fi.write("c in [mol/L]")
        fi.write("\t")
        fi.write(str(c*1000.0))
        
        
        for i in range(3):
            fi.write("\n")
        
        fi.write("Nyquist DATA")
        
        if IMPFIT == 0:
            fi.write("\n")
            fi.write("Frequency [Hz]")
            fi.write("\t")
            fi.write("Z_real_Calc [Ohm]")
            fi.write("\t")
            fi.write("Z_imag_Calc [Ohm]")
            fi.write("\t")
            fi.write("\n")
            for i in range(len(xxx)):
                fi.write(str(np.asscalar((Freq[i].imag)/(2*np.pi))))
                fi.write("\t")
                fi.write(str(np.asscalar(xxx[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(f[i])))
                fi.write("\t")
                fi.write("\n")

        
        if IMPFIT == 1:
            fi.write("\n")
            fi.write("\n")
            fi.write("\n")
            fi.write("Standard dev real part")
            fi.write("\t")
            fi.write(str(np.asscalar(ReStDev)))
            fi.write("\n")
            fi.write("Standard dev imag part")
            fi.write("\t")
            fi.write(str(np.asscalar(ImStDev)))
            fi.write("\n")
            fi.write("\n")
            fi.write("\n")
            
            fi.write("Frequency [Hz]")
            fi.write("\t")
            fi.write("Z_real_Calc [Ohm]")
            fi.write("\t")
            fi.write("Z_imag_Calc [Ohm]")
            fi.write("\t")
            fi.write("Z_real_Meas [Ohm]")
            fi.write("\t")
            fi.write("Z_imag_Meas [Ohm]")
            fi.write("\t")
            fi.write("\n")
            
            for i in range(len(Z_real_Meas)):
                fi.write(str(np.asscalar(Frequency_Array[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_real_Calc[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_imag_Calc[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_real_Meas[i])))
                fi.write("\t")
                fi.write(str(np.asscalar(Z_imag_Meas[i])))
                fi.write("\t")
                fi.write("\n")
    
    root.destroy()


# In[59]:


#============================================================================================================   
#                 -------> JSTT_NNLS_DRT <-----------
#============================================================================================================   
def JSTTNNLS_DRT():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    if File_Was_Loaded ==1:
        JSTTNNLS_DRT_Entry_Window()
        
    
      
def JSTTNNLS_DRT_Entry_Window():         
    Fenster = Toplevel()  
    Fenster.geometry("350x250")
    Fenster.title("JSTT-NNLS-DRT")
    
    Resol_Label = Label(Fenster, text="Increase Data res. by fact.*")
    Resol_Label.grid(row=0, column=0)
    Resol_Eingabe = Entry(Fenster)                                               
    Resol_Eingabe.grid(row=0, column=1)
        
    reg_par_Label = Label(Fenster, text="Tikhonov regul. parameter*")
    reg_par_Label.grid(row=1, column=0)
    reg_par_Eingabe = Entry(Fenster)                                               
    reg_par_Eingabe.grid(row=1, column=1)
    
    Gauss_decay_Label = Label(Fenster,text="Gaussian-decay-factor*")
    Gauss_decay_Label.grid(row=2, column=0)
    Gauss_decay_Eingabe = Entry(Fenster)
    Gauss_decay_Eingabe.grid(row=2, column=1)
    
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Use imaginary part", variable=var1).grid(row=3, column=0, sticky=W)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Use real part", variable=var2).grid(row=4, column=0, sticky=W)
    
    var3 = IntVar()
    Checkbutton(Fenster, text="Save_Data", variable=var3).grid(row=5, column=0, sticky=W)
    
    def AcceptParams():
        global Resol
        global reg_par
        global Gauss_decay
        global Imag_User
        global Real_User
        global As_txt_saver
    

        
        Resol             = (float(Resol_Eingabe.get()))
        reg_par           = (float(reg_par_Eingabe.get()))
        Gauss_decay       = (float(Gauss_decay_Eingabe.get()))
       
        Imag_User         = var1.get()
        Real_User         = var2.get()
        As_txt_saver      = var3.get()
        
        
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=6,column=0)
    button=Button(Fenster,text="Next",command=JSTTNNLS_DRT_Transformation).grid(row=7,column=0)
    
    

    
#========================================================================================================
#========================================================================================================
def JSTTNNLS_DRT_Transformation():
    
    root = Toplevel()
    root.title("DRT-Function")
    
    
    global Z_imag_Meas
    global Z_real_Meas
    global frequencies
    Z_imag_Meas     = np.squeeze(Strom)        #Nach dem Laden heißen die Arrays nur so wie früher
    Z_real_Meas     = np.squeeze(Potenzial)    #deshalb hier die Umbenennung
    frequencies     = Frequency_Array
    
    frequencies     = Frequency_Array  #data[:,0]
    Z               = Z_real_Meas + 1j*Z_imag_Meas
    frequencies     = frequencies[np.imag(Z) < 0]
    Z               = Z[np.imag(Z) < 0]
    
    
    
    #==============================================
    #Interpolate higher resolution
    #==============================================
    
    Interpolation_Imag = InterpolatedUnivariateSpline(np.log10(frequencies[::-1]),Z.imag[::-1])
    Interpolation_Real = InterpolatedUnivariateSpline(np.log10(frequencies[::-1]),Z.real[::-1])
    Resolutor          = Resol     #calculate n-times resolution
    Freq_Resol         = np.logspace(np.log10(frequencies[-1]),np.log10(frequencies[0]),Resolutor*len(frequencies))
    Z_Resol_Imag       = Interpolation_Imag(np.log10(Freq_Resol))
    Z_Resol_Real       = Interpolation_Real(np.log10(Freq_Resol))
    frequencies        = Freq_Resol[::-1]
    Z                  = Z_Resol_Real[::-1] + 1j*Z_Resol_Imag[::-1]
    #==============================================
    
    def make_tau(f_max,f_min,M):
        tau = np.zeros(M)
        tau[0] = 1/(2*np.pi*f_max) #first value in array tau is tau_min
        tau[M-1] = 1/(2*np.pi*f_min) #last value is tau_max
        
        for k in range(1,M-1):
            tau[k] =10**(np.log10(tau[0])+(k)/float(M-1)*np.log10(tau[-1]/tau[0]))
        return tau
    
    
    
    global tau
    tau =make_tau(10**8,10**(-8),len(frequencies))

    
    #==========================================================================================

    def build_matrix(freq,tau):
        '''
        set up the data matrix X needed to compute vector w
        which fullfills Xw=y with y being the experimental Z-data
        '''    
        X = np.zeros((len(freq),len(tau)),dtype=np.complex)
        for m in range(len(freq)):
            for n in range(len(tau)):
                X[m][n] = 1/(1+2j*np.pi*freq[m]*tau[n])
        return X

    X = build_matrix(frequencies,tau)
    
    def compute_DRT(X,data,damp):
        
        A = np.vstack((X,damp*np.ones(X.shape)))
        c = np.squeeze(np.concatenate([data,np.zeros(data.shape)]))

        b = nnls(A,c)[0]
        return b

    
    if Imag_User == 1:
        Im_DRT = compute_DRT(-X.imag,np.abs(Z.imag),reg_par)
    
    if Real_User == 1:
        Re_DRT = compute_DRT(X.real,np.abs((Z.real-Z.real[0])),reg_par) 
    
    
    #==========================================================================================
    #Gaussige Verschmierung
    #==========================================================================================
    Knechtarray_Re = np.zeros((len(tau),len(tau)))
    Knechtarray_Im = np.zeros((len(tau),len(tau)))
    
    for i in range(len(tau)):
        for j in range(len(tau)):
            epsilon = Gauss_decay
            if Real_User == 1:
                Knechtarray_Re[i,j] = Re_DRT[i]*np.exp(-epsilon*(np.abs(np.log10(tau[i])-np.log10(tau[j]))**2)) 
            if Imag_User == 1:
                Knechtarray_Im[i,j] = Im_DRT[i]*np.exp(-epsilon*(np.abs(np.log10(tau[i])-np.log10(tau[j]))**2)) 
    
    global broadened_Re
    global broadened_Im
    broadened_Re = np.zeros(len(tau))
    broadened_Im = np.zeros(len(tau))
    
    for i in range(len(tau)):
        if Real_User == 1:
            broadened_Re[i] = np.sum(Knechtarray_Re[::,i])
        if Imag_User == 1:
            broadened_Im[i] = np.sum(Knechtarray_Im[::,i])

    
    
    Abbildung = Figure(figsize=(5, 4), dpi=100)
    b = Abbildung.add_subplot(111)
    if Imag_User == 1:
        b.plot(np.log10(tau[::]),broadened_Im,color='b', linestyle ='-',label = 'Im_DRT')
    if Real_User == 1:
        b.plot(np.log10(tau[::]),broadened_Re,color='r', linestyle ='-', label = 'Re_DRT')
    b.legend()    
    b.set_xlabel('log10(tau / s)', fontsize=12)
    b.set_ylabel('gamma(tau)', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b.spines[axis].set_linewidth(2)
        b.spines[axis].set_color('k')
    b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    if As_txt_saver == 1:
        JSTTNNLS_DRT_as_Txt_Saver()
    

def JSTTNNLS_DRT_as_Txt_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))

    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")
        
        fi.write("Resol")
        fi.write("\t")
        fi.write(str(Resol))
        fi.write("\n")
        fi.write("Reg_Par")
        fi.write("\t")
        fi.write(str(reg_par))
        fi.write("\n")
        fi.write("Gauss_Decay")
        fi.write("\t")
        fi.write(str(Gauss_decay))
        fi.write("\n")
        
        fi.write("\n")
        fi.write("\n")
        fi.write("log10(tau / s)")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("gamma(tau)_Imag")
            fi.write("\t")
        if Real_User == 1:
            fi.write("gamma(tau)_Real")
            fi.write("\n")
        
        for i in range(len(tau)):
            fi.write(str(np.asscalar(np.log10(tau[i]))))
            fi.write("\t")
            if Imag_User == 1:
                fi.write(str(np.asscalar(broadened_Im[i])))
                fi.write("\t")
            if Real_User == 1:
                fi.write(str(np.asscalar(broadened_Re[i])))
                fi.write("\t")
            fi.write("\n")
    
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency [Hz]")
        fi.write("\t")
        fi.write("Z_real_Meas [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_Meas [Ohm]")
        fi.write("\t")
        fi.write("\n")
            
        for i in range(len(Z_real_Meas)):
            fi.write(str(np.asscalar(Frequency_Array[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_real_Meas[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_imag_Meas[i])))
            fi.write("\t")
            fi.write("\n")
    
    root.destroy()


# In[60]:


#============================================================================================================   
#                 -------> DRT_TOOLS_NNLS_DRT <-----------
#============================================================================================================ 
def DRT_Tools_NNLS_DRT():
    if File_Was_Loaded ==0:
        No_Loaded_File_Warner()
    if File_Was_Loaded ==1:
        DRT_Tools_NNLS_DRT_Entry_Window()

#============================================================================================================   
#    -------> define all DRT Tools functions in python 2.7 language except quadprog <-----------
#       as there is not quadprog in python this transfom uses nnls algorthm in the end
#============================================================================================================ 

        
def quad_format_combined(A_re,A_im,b_re,b_im,M_re,M_im,damp): 
    H = 2*(0.5*(np.dot(A_re.T,A_re)+np.dot(A_im.T,A_im))+damp*M_re)
    c = -2*0.5*(np.dot(b_im.T,A_im)+np.dot(b_re.T,A_re))
    return H , c

def quad_format(A,b,M,damp):
    H = 2*(np.dot(A.T,A)+damp*M)
    c = -2*np.dot(b.T,A)
    return H , c    

##############################################################################
#map array to gamma
##############################################################################

def map_array_to_gamma(freq_map, freq_coll, x, epsilon, rbf_type):
    if rbf_type == "gaussian":
        rbf = lambda y,y0: np.exp(-(epsilon*(y-y0))**2)
    elif rbf_type == "C0_matern":
        rbf = lambda y,y0: np.exp(-abs(epsilon*(y-y0)))
    elif rbf_type == "C2_matern":
        rbf = lambda y,y0: np.exp(-abs(epsilon*(y-y0)))*(1+abs(epsilon*(y-y0)))
    elif rbf_type == "C4_matern":
        rbf = lambda y,y0: 1/3.0*np.exp(-abs(epsilon*(y-y0)))*(3+3*abs(epsilon*(y-y0))+abs(epsilon*(y-y0))**2)
    elif rbf_type == "C6_matern":
        rbf = lambda y,y0: 1/15.0*np.exp(-abs(epsilon*(y-y0)))*(15+15*abs(epsilon*(y-y0))+6*abs(epsilon*(y-y0))**2+abs(epsilon*(y-y0))**3)          
    elif rbf_type == "inverse_quadratic":
        rbf = lambda y,y0: 1.0/(1+(epsilon*(y-y0))**2)
    else:
        print('ERROR - Unexpected RBF input at map_array_to_gamma')
    y0  = -np.log(freq_coll)
    out_gamma = np.zeros(len(freq_map))   
    
    
    for iter_freq_map in range(len(freq_map)):
        freq_map_loc = freq_map[iter_freq_map]
        y = -np.log(freq_map_loc)
        rbf_temp = rbf(y,y0)
        out_gamma[iter_freq_map] = np.dot(x.T,rbf_temp)
    
    return out_gamma

##############################################################################
#inner products
##############################################################################
def inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type):
    a = epsilon*np.log(freq_n/freq_m)
    if rbf_type == "gaussian":
        out_IP = epsilon**3*(3-6*a**2+a**4)*np.exp(-(a**2/2))*(np.pi/2)**0.5
        return out_IP
    elif rbf_type == "C0_matern":
        out_IP = epsilon**3*(1+abs(a))*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C2_matern":
        out_IP = epsilon**3/6.0*(3+3*abs(a)-6*abs(a)**2+abs(a)**3)*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C4_matern":
        out_IP = epsilon**3/30.0*(45+45*abs(a)-15*abs(a)**3-5*abs(a)**4+abs(a)**5)*np.exp(-abs(a))                     
        return out_IP
    elif rbf_type == "C6_matern":
        out_IP = epsilon**3/140.0*(2835+2835*abs(a)+630*abs(a)**2-315*abs(a)**3-210*abs(a)**4-42*abs(a)**5+abs(a)**7)*np.exp(-abs(a))         
        return out_IP
    elif rbf_type == "inverse_quadratic":
        out_IP = 48*(16+5*a**2*(-8 + a**2))*np.pi*epsilon**3/((4 + a**2)**5)
        return out_IP
    else:
        print('ERROR - Unexpected RBF input at inner_prod_rbf_2')

def inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type):
    a = epsilon*np.log(freq_n/freq_m)

    if rbf_type == "gaussian":
        out_IP = -epsilon*(-1+a**2)*np.exp(-(a**2/2))*(np.pi/2)**0.5
        return out_IP
    elif rbf_type == "C0_matern":
        out_IP = epsilon*(1-abs(a))*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C2_matern":
        out_IP = epsilon/6.0*(3+3*abs(a)-abs(a)**3)*np.exp(-abs(a))
        return out_IP
    elif rbf_type == "C4_matern":
        out_IP = epsilon/30.0*(105+105*abs(a)+30*abs(a)**2-5*abs(a)**3-5*abs(a)**4-abs(a)**5)*np.exp(-abs(a))                     
        return out_IP
    elif rbf_type == "C6_matern":
        out_IP = epsilon/140.0*(10395 +10395*abs(a)+3780*abs(a)**2+315*abs(a)**3-210*abs(a)**4-84*abs(a)**5-14*abs(a)**6-abs(a)**7)*np.exp(-abs(a))         
        return out_IP
    elif rbf_type == "inverse_quadratic":
        out_IP = 4*epsilon*(4-3*a**2)*np.pi/((4+a**2)**3)
        return out_IP
    else:
        print('ERROR - Unexpected RBF input at inner_prod_rbf.')
    
#############################################################################
#g_i
#############################################################################
def g_i(freq_n, freq_m, epsilon, rbf_type):
    alpha = 2*np.pi*freq_n/freq_m

    if rbf_type == "gaussian":
        rbf = lambda x: np.exp(-(epsilon*x)**2)
    elif rbf_type == "C0_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))
    elif rbf_type == "C2_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))*(1+abs(epsilon*x))
    elif rbf_type == "C4_matern":
        rbf = lambda x: (1/3.0)*np.exp(-abs(epsilon*x))*(3+3*abs(epsilon*x)+abs(epsilon*x)**2)
    elif rbf_type == "C6_matern":
        rbf = lambda x: (1/15.0)*np.exp(-abs(epsilon*x))*(15+15*abs(epsilon*x)+6*abs(epsilon*x)**2+abs(epsilon*x)**3)          
    elif rbf_type == "inverse_quadratic":
        rbf = lambda x: 1.0/(1+(epsilon*x)**2)
    else:
        print('ERROR - Unexpected RBF input at g_i')
    
    integrand_g_i = lambda x: 1.0/(1+alpha**2 *np.exp(2*x))*rbf(x)
    
    out_val  = quad(integrand_g_i, -100, 100,epsabs=1.0e-06, epsrel=1.0e-06)[0]#,'RelTol',1E-9,'AbsTol',1e-9);       
    return out_val 


#############################################################################
#g_ii
#############################################################################

def g_ii(freq_n, freq_m, epsilon, rbf_type):
    alpha = 2*np.pi*freq_n/freq_m

    if rbf_type == "gaussian":
        rbf = lambda x: np.exp(-(epsilon*x)**2)
    elif rbf_type == "C0_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))
    elif rbf_type == "C2_matern":
        rbf = lambda x: np.exp(-abs(epsilon*x))*(1+abs(epsilon*x))
    elif rbf_type == "C4_matern":
        rbf = lambda x: (1/3.0)*np.exp(-abs(epsilon*x))*(3+3*abs(epsilon*x)+abs(epsilon*x)**2)
    elif rbf_type == "C6_matern":
        rbf = lambda x: (1/15.0)*np.exp(-abs(epsilon*x))*(15+15*abs(epsilon*x)+6*abs(epsilon*x)**2+abs(epsilon*x)**3)          
    elif rbf_type == "inverse_quadratic":
        rbf = lambda x: 1.0/(1+(epsilon*x)**2)

    else:
        print('ERROR - Unexpected RBF input at g_ii')
    
    integrand_g_ii = lambda x: alpha/(1.0/np.exp(x)+alpha**2 *np.exp(x))*rbf(x)

    out_val  = quad(integrand_g_ii, -100, 100,epsabs=1.0e-06, epsrel=1.0e-06)[0]       
    return out_val 
    
#############################################################################
#compute_L_re
#############################################################################
def compute_L_re(freq):
    tau    = 1.0/freq
    N_freq = len(freq)
    out_L      = np.zeros((N_freq-1,N_freq+2))
    out_L_temp = np.zeros((N_freq-1,N_freq+1))
    
    for p in range(N_freq-1):
        delta_loc          = np.log(tau[p+1]/tau[p])
        out_L_temp[p,p+1]  = -1/delta_loc
        out_L_temp[p,p+2]  =  1/delta_loc
        
    out_L[::,1::] = out_L_temp   
    return out_L

#############################################################################
#compute_L_im
############################################################################
def compute_L_im(freq):
    tau    = 1.0/freq
    N_freq = len(freq)
    
    out_L      = np.zeros((N_freq-1,N_freq+2))
    out_L_temp = np.zeros((N_freq-1,N_freq))
    
    for p in range(N_freq-1):
        delta_loc        = np.log(tau[p+1]/tau[p])
        out_L_temp[p,p]  = -1/delta_loc
        out_L_temp[p,p+1]=  1/delta_loc
        
    out_L[::,2::] = out_L_temp   
    return out_L

#############################################################################
#assemble A
############################################################################

def assemble_A_im(freq, epsilon, rbf_type,L=0):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_A_im_temp = np.zeros(len(freq))
    out_A_im      = np.zeros((len(freq), len(freq)+2))
    if std_freq/mean_freq < 1:  
        for iter_freq_n in range(len(freq)): 
            freq_n     = freq[iter_freq_n]
            freq_m = freq[0]
            C[iter_freq_n, 0] = g_ii(freq_n, freq_m, epsilon, rbf_type) 

        for iter_freq_m in range(len(freq)):
            freq_n = freq[0]
            freq_m = freq[iter_freq_m]
            R[0, iter_freq_m] = g_ii(freq_n, freq_m, epsilon, rbf_type)

        out_A_im_temp = toeplitz(C,R)
    
    else:
        for iter_freq_n in range(len(freq)):
            for iter_freq_m in range(len(freq)):
                freq_n = freq[iter_freq_n] 
                freq_m = freq[iter_freq_m]
                out_A_im_temp[iter_freq_n, iter_freq_m] = g_ii(freq_n, freq_m, epsilon, rbf_type) 
    out_A_im[:, 2::] = out_A_im_temp
    if L==1:
        out_A_im[:,0] = -2*np.pi*(freq[:])         
    return out_A_im


def assemble_A_re(freq, epsilon, rbf_type):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_A_re_temp = np.zeros(len(freq))
    out_A_re      = np.zeros((len(freq), len(freq)+2))
    
    if std_freq/mean_freq < 1:  #(error in frequency difference <1% make sure that the terms are evenly distributed)
        for iter_freq_n in range(len(freq)): 
            freq_n     = freq[iter_freq_n]
            freq_m = freq[0]
            C[iter_freq_n, 0] = g_i(freq_n, freq_m, epsilon, rbf_type) 

        for iter_freq_m in range(len(freq)):
            freq_n = freq[0]
            freq_m = freq[iter_freq_m]
            R[0, iter_freq_m] = g_i(freq_n, freq_m, epsilon, rbf_type)

        out_A_re_temp = toeplitz(C,R)
    
    else:
        for iter_freq_n in range(len(freq)):
            for iter_freq_m in range(len(freq)):
                freq_n = freq[iter_freq_n] 
                freq_m = freq[iter_freq_m]
                out_A_re_temp[iter_freq_n, iter_freq_m] = g_i(freq_n, freq_m, epsilon, rbf_type)
    
    out_A_re[:, 2::] = out_A_re_temp
    out_A_re[:,1] = 1
    return out_A_re
    
##############################################################################
#assemble M
##############################################################################
def assemble_M_re(freq, epsilon, rbf_type, der_used):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_M_re_temp = np.zeros(len(freq))
    out_M_re      = np.zeros((len(freq)+2, len(freq)+2))
    
    if der_used == "1st-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
    
            for iter_freq_m in range(len(freq)):
                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
                
            out_M_re_temp = toeplitz(C,R)

        else:
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_re_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
                    
    #-------------------------------------------------------------------------------------------------------------------           
       
    if der_used == "2nd-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1  #% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                    
            for iter_freq_m in range(len(freq)):

                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                
            out_M_re_temp = toeplitz(C,R)

        else: #%if log of tau is not evenly distributed
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_re_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)

    out_M_re[2::, 2::] = out_M_re_temp
    
    return out_M_re


def assemble_M_im(freq, epsilon, rbf_type, der_used):
    std_freq      = np.std(np.diff(np.log(1.0/freq)))
    mean_freq     = np.mean(np.diff(np.log(1.0/freq)))
    R             = np.zeros((1,len(freq)))
    C             = np.zeros((len(freq),1))
    out_M_im_temp = np.zeros(len(freq))
    out_M_im      = np.zeros((len(freq)+2, len(freq)+2))

    if der_used == "1st-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
    
            for iter_freq_m in range(len(freq)):
                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
            
            out_M_im_temp = toeplitz(C,R)

        else:
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_im_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf(freq_n, freq_m, epsilon, rbf_type)
                    
    #-------------------------------------------------------------------------------------------------------------------           
       
    if der_used == "2nd-order":
        if std_freq/mean_freq < 1:  #%(error in frequency difference <1  #% make sure that the terms are evenly distributed)
            for iter_freq_n in range(len(freq)):
                freq_n = freq[iter_freq_n]
                freq_m = freq[0]
                C[iter_freq_n, 0] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                    
            for iter_freq_m in range(len(freq)):

                freq_n = freq[0]
                freq_m = freq[iter_freq_m]
                R[0, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)
                
            out_M_im_temp = toeplitz(C,R)

        else: #%if log of tau is not evenly distributed
            for iter_freq_n in range(len(freq)):
                for iter_freq_m in range(len(freq)):
                    freq_n = freq[iter_freq_n]
                    freq_m = freq[iter_freq_m]
                    out_M_im_temp[iter_freq_n, iter_freq_m] = inner_prod_rbf_2(freq_n, freq_m, epsilon, rbf_type)

    out_M_im[2::, 2::] = out_M_im_temp
    
    return out_M_im

##############################################################################
#----------------------------------------------------------------------------#
##############################################################################
    
def compute_A_im(freq,L):
    omega  = 2*np.pi*freq
    tau    = 1.0/freq
    N_freq = len(freq)
    
    out_A_im = np.zeros((N_freq,N_freq+2))
    if L == 1:
        out_A_im[::,0] = -2*np.pi*freq[::]
    
    for p in range(N_freq):
        for q in range(N_freq):
            if  q ==0:
                out_A_im[p, q+2] = 0.5*((omega[p]*tau[q])/(1+(omega[p]*tau[q])**2)) *np.log(tau[q+1]/tau[q])
            elif q == N_freq-1:
                out_A_im[p, q+2] = 0.5*((omega[p]*tau[q])/(1+(omega[p]*tau[q])**2))*np.log(tau[q]/tau[q-1])
            else:
                out_A_im[p, q+2] = 0.5*((omega[p]*tau[q])/(1+(omega[p]*tau[q])**2))*np.log(tau[q+1]/tau[q-1])
    
    return out_A_im

def compute_A_re(freq):
    omega  = 2*np.pi*freq
    tau    = 1.0/freq
    N_freq = len(freq)
    
    out_A_re = np.zeros((N_freq,N_freq+2))
    out_A_re[::,1] = 1
    
    for p in range(N_freq):
        for q in range(N_freq):
            if  q ==0:
                out_A_re[p, q+2] = 0.5*((1.0)/(1+(omega[p]*tau[q])**2)) *np.log(tau[q+1]/tau[q])
            elif q == N_freq-1:
                out_A_re[p, q+2] = 0.5*((1.0)/(1+(omega[p]*tau[q])**2))*np.log(tau[q]/tau[q-1])
            else:
                out_A_re[p, q+2] = 0.5*((1.0)/(1+(omega[p]*tau[q])**2))*np.log(tau[q+1]/tau[q-1])
    return out_A_re

#=====================================================================================================================
#translation of all DRT functions from matlab finished. The main-window will be defined as separate function
#=====================================================================================================================
        
    
def DRT_Tools_NNLS_DRT_Entry_Window():
    Fenster = Toplevel()  
    Fenster.geometry("350x450")
    Fenster.title("DRT_Tools_NNLS_DRT")
    
    Resol_Label = Label(Fenster, text="Increase Data res. by fact.*")
    Resol_Label.grid(row=0, column=0)
    Resol_Eingabe = Entry(Fenster)                                               
    Resol_Eingabe.grid(row=0, column=1)
        
    reg_par_Label = Label(Fenster, text="Tikhonov regul. parameter*")
    reg_par_Label.grid(row=1, column=0)
    reg_par_Eingabe = Entry(Fenster)                                               
    reg_par_Eingabe.grid(row=1, column=1)
    
    epsilon_Label = Label(Fenster,text="epsilon*")
    epsilon_Label.grid(row=2, column=0)
    epsilon_Eingabe = Entry(Fenster)
    epsilon_Eingabe.grid(row=2, column=1)
    
    var1 = IntVar()
    Checkbutton(Fenster, text="Use imaginary part", variable=var1).grid(row=3, column=0, sticky=W)
    
    var2 = IntVar()
    Checkbutton(Fenster, text="Use real part", variable=var2).grid(row=4, column=0, sticky=W)
    
    var3 = IntVar()
    Checkbutton(Fenster, text="Combined re_im", variable=var3).grid(row=5, column=0, sticky=W)
    
    Sep_Label = Label(Fenster,text="Choose ONE type of radial basis function*")
    Sep_Label.grid(row=6, column=0)
    
    var4 = IntVar()
    Checkbutton(Fenster, text="gaussian", variable=var4).grid(row=7, column=0, sticky=W)
    
    var5 = IntVar()
    Checkbutton(Fenster, text="C2_Matern", variable=var5).grid(row=8, column=0, sticky=W)
    
    var6 = IntVar()
    Checkbutton(Fenster, text="C4_Matern", variable=var6).grid(row=9, column=0, sticky=W)
    
    var7 = IntVar()
    Checkbutton(Fenster, text="C6_Matern", variable=var7).grid(row=10, column=0, sticky=W)
    
    var8 = IntVar()
    Checkbutton(Fenster, text="Inv._quadratic", variable=var8).grid(row=11, column=0, sticky=W)
    
    Sep2_Label = Label(Fenster,text="Choose ONE type of derivative order*")
    Sep2_Label.grid(row=12, column=0)
    
    var9 = IntVar()
    Checkbutton(Fenster, text="1st_order", variable=var9).grid(row=13, column=0, sticky=W)
    
    var10 = IntVar()
    Checkbutton(Fenster, text="2nd_order", variable=var10).grid(row=14, column=0, sticky=W)
    
    var11 = IntVar()
    Checkbutton(Fenster, text="Save_Data", variable=var11).grid(row=15, column=0, sticky=W)
    
    
    def AcceptParams():
        global Resol
        global Tikh_Par
        global eps
        global Imag_User
        global Real_User
        global Comb_User
        global deru
        global RBF_Type
        global As_txt_saver
    

        Resol             = (float(Resol_Eingabe.get()))
        Tikh_Par          = (float(reg_par_Eingabe.get()))
        eps               = (float(epsilon_Eingabe.get()))
       
        Imag_User         = var1.get()
        Real_User         = var2.get()
        Comb_User         = var3.get()
        gaussrbf          = var4.get()
        c2rbf             = var5.get()
        c4rbf             = var6.get()
        c6rbf             = var7.get()
        invqurbf          = var8.get()
        deru1             = var9.get()
        deru2             = var10.get()
        
        #rbf_type
        if gaussrbf == 1:
            RBF_Type = 'gaussian'
        if c2rbf == 1:
            RBF_Type = 'C2_matern'
        if c4rbf == 1:
            RBF_Type = 'C4_matern'
        if c6rbf == 1:
            RBF_Type = 'C6_matern'
        if invqurbf  == 1:
            RBF_Type = 'inverse_quadratic'
          
        #der_used
        if deru1  == 1:
            deru = '1st-order'
        if deru2  == 1:
            deru = '2nd-order'
            
    
        As_txt_saver      = var11.get()
        
        
    
    button=Button(Fenster,text="Accept Parameters",command=AcceptParams).grid(row=16,column=0)
    button=Button(Fenster,text="Next",command=DRT_Tools_NNLS_DRT_Transformation).grid(row=17,column=0)

#====================================================================================================================    
#here is the GUI translation in python 2.7
#====================================================================================================================    


def DRT_Tools_NNLS_DRT_Transformation():
    ##############################################################################
    
    root = Toplevel()
    root.title("DRT_Tools_DRT-Function")
    
    
    global Z_imag_Meas
    global Z_real_Meas
    global freq
    Z_imag_Meas     = np.squeeze(Strom)        #Nach dem Laden heißen die Arrays nur so wie früher
    Z_real_Meas     = np.squeeze(Potenzial)    #deshalb hier die Umbenennung
    freq            = np.squeeze(Frequency_Array)
                       
    
    Z               = Z_real_Meas + 1j*Z_imag_Meas
    freq            = freq[np.imag(Z) < 0]
    Z               = Z[np.imag(Z) < 0]
  
    
    ##############################################################################
    #initialise parameters and bounds
    epsilon       = eps
    rbf_type      = RBF_Type
    der_used      = deru
    damp          = Tikh_Par

    taumax        = ceil(np.max(np.log10(1./freq)))+1
    taumin        = floor(np.min(np.log10(1./freq)))-1
    freq_out      = np.logspace(-taumin,-taumax,Resol*len(freq))

    ##############################################################################
    A_im          = assemble_A_im(freq,epsilon,rbf_type)
    M_im          = assemble_M_im(freq,epsilon,rbf_type,der_used)
    H_im, f_im    = quad_format(A_im,-Z.imag,M_im,damp)

    A_re          = assemble_A_re(freq,epsilon,rbf_type)
    M_re          = assemble_M_re(freq,epsilon,rbf_type,der_used)
    H_re, f_re    = quad_format(A_re,(Z.real),M_re,damp)

    H_co, f_co    = quad_format_combined(A_re, A_im,(Z.real), -Z.imag, M_re, M_im,damp)

    x_imag_2      = nnls(H_im,np.abs(f_im))[0]
    x_real_2      = nnls(H_re,np.abs(f_re))[0]
    x_comb_2      = nnls(H_co,np.abs(f_co))[0]
    
    global gamma_imag_fine
    global gamma_real_fine
    global gamma_comb_fine
    global TAUARRAY
    
    TAUARRAY         = 1./freq_out 

    gamma_imag_fine  = map_array_to_gamma(freq_out,freq,x_imag_2[2:],epsilon,rbf_type)
    gamma_real_fine  = map_array_to_gamma(freq_out,freq,x_real_2[2:],epsilon,rbf_type)
    gamma_comb_fine  = map_array_to_gamma(freq_out,freq,x_comb_2[2:],epsilon,rbf_type)
    
    
    Abbildung = Figure(figsize=(5, 4), dpi=100)
    b = Abbildung.add_subplot(111)
    if Imag_User == 1:
        b.plot(np.log10(TAUARRAY),gamma_imag_fine,color='b', linestyle ='-',label = 'Im_DRT')
    if Real_User == 1:
        b.plot(np.log10(TAUARRAY),gamma_real_fine,color='r', linestyle ='-', label = 'Re_DRT')
    if Comb_User == 1:
        b.plot(np.log10(TAUARRAY),gamma_comb_fine,color='k', linestyle ='-', label = 'Comb_DRT')
    
    b.legend()    
    b.set_xlabel('log10(tau / s)', fontsize=12)
    b.set_ylabel('gamma(tau)', fontsize=12)
    for axis in ['top','bottom','left','right']:
        b.spines[axis].set_linewidth(2)
        b.spines[axis].set_color('k')
    b.tick_params(tickdir='in', width = 2, length=6, labelsize=12)
    
    canvas = FigureCanvasTkAgg(Abbildung, master=root)
    canvas.show()
    canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

    toolbar = NavigationToolbar2TkAgg(canvas, root)
    toolbar.update()
    canvas._tkcanvas.pack(side=TOP, fill=BOTH, expand=1)
    
    if As_txt_saver == 1:
        DRT_Tools_NNLS_DRT_as_Txt_Saver()
        
        
def DRT_Tools_NNLS_DRT_as_Txt_Saver():
    root = Toplevel() 
    root.filename = asksaveasfilename(title = "Save Your Data",filetypes = (("save as","*.txt"),("all files","*.*")))

    with open(root.filename,"w") as fi:
        fi.write("Parameters")
        fi.write("\n")
        fi.write("\n")
        
        fi.write("DRT calc. from")
        fi.write("\n") 
        fi.write("Imaginary-Part:")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("yes")
        if Imag_User != 1:
            fi.write("no")
        fi.write("\n") 
        fi.write("Real-Part:")
        fi.write("\t")
        if Real_User == 1:
            fi.write("yes")
        if Real_User != 1:
            fi.write("no")
        fi.write("\n")
        fi.write("Combined Real_Im:")
        fi.write("\t")
        if Comb_User == 1:
            fi.write("yes")
        if Comb_User != 1:
            fi.write("no")
        fi.write("\n")
        
        fi.write("Resolution increase")
        fi.write("\t")
        fi.write(str(Resol))
        fi.write("\n")
        fi.write("Reg_Par")
        fi.write("\t")
        fi.write(str(Tikh_Par))
        fi.write("\n")
        fi.write("epsilon")
        fi.write("\t")
        fi.write(str(eps))
        fi.write("\n")
        fi.write("RBF-Type:")
        fi.write("\t")
        if RBF_Type == 'gaussian':
            fi.write("gaussian")
        if RBF_Type == 'C2_matern':
            fi.write("C2_matern")
        if RBF_Type == 'C4_matern':
            fi.write("C4_matern")
        if RBF_Type == 'C6_matern':
            fi.write("C6_matern")
        if RBF_Type == 'inverse_quadratic':
            fi.write("inverse_quadratic")
        fi.write("\n")
        fi.write("Derivative-Order:")
        fi.write("\t")
        if deru == '1st-order':
            fi.write("1st-order")
        if deru == '2nd-order':
            fi.write("2nd-order")
    
        
        fi.write("\n")
        fi.write("\n")
        fi.write("log10(tau / s)")
        fi.write("\t")
        if Imag_User == 1:
            fi.write("gamma(tau)_Imag")
            fi.write("\t")
        if Real_User == 1:
            fi.write("gamma(tau)_Real")
            fi.write("\t")
        if Comb_User == 1:
            fi.write("gamma(tau)_Combined")
        fi.write("\n")
        
        for i in range(len(TAUARRAY)):
            fi.write(str(np.asscalar(np.log10(TAUARRAY[i]))))
            fi.write("\t")
            if Imag_User == 1:
                fi.write(str(np.asscalar(gamma_imag_fine[i])))
                fi.write("\t")
            if Real_User == 1:
                fi.write(str(np.asscalar(gamma_real_fine[i])))
                fi.write("\t")
            if Comb_User == 1:
                fi.write(str(np.asscalar(gamma_comb_fine[i])))
                fi.write("\t")
            fi.write("\n")
    
        fi.write("\n")
        fi.write("\n")
        fi.write("\n")
        fi.write("Frequency [Hz]")
        fi.write("\t")
        fi.write("Z_real_Meas [Ohm]")
        fi.write("\t")
        fi.write("Z_imag_Meas [Ohm]")
        fi.write("\t")
        fi.write("\n")
            
        for i in range(len(Z_real_Meas)):
            fi.write(str(np.asscalar(Frequency_Array[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_real_Meas[i])))
            fi.write("\t")
            fi.write(str(np.asscalar(Z_imag_Meas[i])))
            fi.write("\t")
            fi.write("\n")
    
    root.destroy()


# In[61]:


root = Tk()
root.title("Polarographica")                         
root.geometry("800x500")

menu = Menu(root)                                               #Genereller Befehl um neuen Reiter zu erzeugen
root.config(menu=menu)


filemenu = Menu(menu)                                           #Reiter Filemenu und alles was es da drin gibt steht drunter
menu.add_cascade(label="File", menu=filemenu)           
filemenu.add_command(label="Open CV/LSV", command=LSV_DATA)
filemenu.add_command(label="Open Chronoamperometry", command=CA_DATA)
filemenu.add_command(label="Open Impedance/Nyquist", command=IMP_DATA)


SemiInfDiff   = Menu(menu)
menu.add_cascade(label="Classical Evaluations", menu=SemiInfDiff)
SemiInfDiff.add_command(label="Randles Sevcik Reversible", command=RandlesRevWindow)
SemiInfDiff.add_command(label="Randles Sevcik Irreversible", command=RS_IRR_Window)
SemiInfDiff.add_command(label="Koutecky-Levich Analysis", command=KouteckyLevichWindow)
SemiInfDiff.add_command(label="Tafel and Shape-Analysis", command=TafelWindow)
SemiInfDiff.add_command(label="Cottrell-Analysis", command=CottrellWindow)

Simulators   = Menu(menu)
menu.add_cascade(label="Stat.CV-Simulators", menu=Simulators)
Simulators.add_command(label="Planar Semi-Infinite Diffusion CV Simulator", command=Semi_Inf_Planar)
Simulators.add_command(label="Planar Finite Diffusion CV Simulator", command=Finit_Planar)
Simulators.add_command(label="External Cylindrical Semi-Infinite Diffusion CV Simulator", command=Semi_Inf_Zyl_Ext)
Simulators.add_command(label="External Cylindrical Finite Diffusion CV Simulator", command=Finit_Zyl_Ext)
Simulators.add_command(label="Internal Cylindrical Finite Diffusion CV Simulator", command=Finit_Zyl_Int)
Simulators.add_command(label="External Spherical Semi-Infinite Diffusion CV Simulator", command=Semi_Inf_Sphere_Ext)
Simulators.add_command(label="Internal Spherical Finite Diffusion CV Simulator", command=Finit_Sphere_Int)
Simulators.add_command(label="--------------------------------------------------------------------", command=Finit_Sphere_Int)
Simulators.add_command(label="Statistically Weighted Planar Finite Diffusion CV Simulator", command=Statistical_Finit_Planar)
Simulators.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion CV Simulator", command=Statistical_Finit_Zyl_Ext)
Simulators.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion CV Simulator", command=Statistical_Finit_Zyl_Int)
Simulators.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion CV Simulator", command=Statistical_Finit_Sphere_Int)

Fitters   = Menu(menu)
menu.add_cascade(label="Stat.CV-Fitters", menu=Fitters)
Fitters.add_command(label="Planar Semi-Infinite Diffusion CV Fitter", command=Semi_Inf_Planar_FITTER)
Fitters.add_command(label="Planar Finite Diffusion CV Fitter", command=Finit_Planar_FITTER)
Fitters.add_command(label="External Cylindrical Semi-Infinite Diffusion CV Fitter", command=Semi_Inf_Zyl_Ext_FITTER)
Fitters.add_command(label="External Cylindrical Finite Diffusion CV Fitter", command=Finit_Zyl_Ext_FITTER)
Fitters.add_command(label="Internal Cylindrical Finite Diffusion CV Fitter", command=Finit_Zyl_Int_FITTER)
Fitters.add_command(label="External Spherical Semi-Infinite Diffusion CV Fitter", command=Semi_Inf_Sphere_Ext_FITTER)
Fitters.add_command(label="Internal Spherical Finite Diffusion CV Fitter", command=Finit_Sphere_Int_FITTER)
Fitters.add_command(label="--------------------------------------------------------------------")
Fitters.add_command(label="Statistically Weighted Planar Finite Diffusion CV Fitter", command=Statistical_Finit_Planar_FITTER)
Fitters.add_command(label="Statistically Weighted External Cylindrical Finite Diffusion CV Fitter", command=Statistical_Finit_Zyl_Ext_FITTER)
Fitters.add_command(label="Statistically Weighted Internal Cylindrical Finite Diffusion CV Fitter", command=Statistical_Finit_Zyl_Int_FITTER)
Fitters.add_command(label="Statistically Weighted Internal Spherical Finite Diffusion CV Fitter", command=Statistical_Finit_Sphere_Int_FITTER)



Imp_simulators   = Menu(menu)
menu.add_cascade(label="Impedance Simulators", menu=Imp_simulators)
Imp_simulators.add_command(label="Planar Semi-Infinite Diffusion Impedance Simulator", command=Semi_Inf_Plan_Imp_Sim)
Imp_simulators.add_command(label="Planar Finite Transmissive Diffusion Impedance Simulator", command=Fin_Plan_Trans_Imp_Sim)
Imp_simulators.add_command(label="Planar Finite Reflective Diffusion Impedance Simulator", command=Fin_Plan_Ref_Imp_Sim)
Imp_simulators.add_command(label="Zylindrical Semi-Infinite Diffusion Impedance Simulator", command=Semi_Inf_Zyl_Imp_Sim)
Imp_simulators.add_command(label="Spherical Semi-Infinite Diffusion Impedance Simulator", command=Semi_Inf_Sphe_Imp_Sim)
Imp_simulators.add_command(label="Spherical Internal Finite Diffusion Impedance Simulator", command=Finit_Int_Sphe_Imp_Sim)


Imp_fitters   = Menu(menu)
menu.add_cascade(label="Impedance Fitters", menu=Imp_fitters)
Imp_fitters.add_command(label="Planar Semi-Infinite Diffusion Impedance Fitter", command=Semi_Inf_Plan_Imp_Fit)
Imp_fitters.add_command(label="Planar Finite Transmissive Diffusion Impedance Fitter", command=Fin_Plan_Trans_Imp_Fit)
Imp_fitters.add_command(label="Planar Finite Reflective Diffusion Impedance Fitter", command=Fin_Plan_Ref_Imp_Fit)
Imp_fitters.add_command(label="Zylindrical Semi-Infinite Diffusion Impedance Fitter", command=Semi_Inf_Zyl_Imp_Fit)
Imp_fitters.add_command(label="Spherical Semi-Infinite Diffusion Impedance Fitter", command=Semi_Inf_Sphe_Imp_Fit)
Imp_fitters.add_command(label="Spherical Internal Finite Diffusion Impedance Fitter", command=Finit_Int_Sphe_Imp_Fit)

DRT_Transformation  = Menu(menu)
menu.add_cascade(label="DRT-Transformation", menu=DRT_Transformation)
DRT_Transformation.add_command(label="DRT-Tools-NNLS-DRT", command=DRT_Tools_NNLS_DRT)
DRT_Transformation.add_command(label="JSTT-NNLS-DRT", command=JSTTNNLS_DRT)



mainloop()


# In[ ]:





# In[ ]:





# In[ ]:




