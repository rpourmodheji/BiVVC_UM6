import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib import pylab as plt
from matplotlib import rc
import matplotlib
import pandas as pd
import math


########################### For Plots ###################################
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
matplotlib.rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
matplotlib.rcParams['mathtext.it'] = 'Bitstream Vera Sans:italic'
matplotlib.rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'
matplotlib.rcParams['font.size'] = 14
rc('text', usetex=False)
Pa_to_mmHg = 0.00750062
mmHg_to_Pa = 1.0/Pa_to_mmHg
kPa_to_mmHg = Pa_to_mmHg*1000.0
mmHg_to_kPa = 1.0/kPa_to_mmHg

AO_npzfile = np.load("UM6_AO_DVP.npz")
PA_npzfile = np.load("UM6_PA_DVP.npz")

DAO = AO_npzfile["Diameter"]
VAO = AO_npzfile["Volume"]
PAO = AO_npzfile["Pressure"]
DPA = PA_npzfile["Diameter"]
VPA = PA_npzfile["Volume"]
PPA = PA_npzfile["Pressure"]


def PA_L_func(V, a, b):
    V = np.array(V)
    return a*V+b

def AO_L_func(V, a, b):
    V = np.array(V)
    return a*V+b

def PAAOV0_L_fun(V,C,V0):
    V = np.array(V)
    return (1.0/C)*(V-V0)





plt.figure(1)
plt.plot(VPA,PPA,'-o')
VPA_ED = np.min(VPA)
VPA_ES = np.max(VPA)
PPA_ED = np.min(PPA)
PPA_ES = np.max(PPA)
popt, _ = curve_fit(PA_L_func, [VPA_ED,VPA_ES], [PPA_ED,PPA_ES])
print(popt)
plt.plot([VPA_ED,VPA_ES],PA_L_func([VPA_ED,VPA_ES],popt[0],popt[1]),'-')
VPA0 = -popt[1]/popt[0]
print(VPA0)
print(PA_L_func(VPA0,popt[0],popt[1]))
plt.plot([VPA0,VPA_ES],PA_L_func([VPA0,VPA_ES],popt[0],popt[1]),'-')
C = 1.0/popt[0]
plt.plot([VPA0,VPA_ES],PAAOV0_L_fun([VPA0,VPA_ES],C,VPA0),'-')
print(f"C for PA is {C}")


plt.figure(2)
plt.plot(VAO,PAO,'-o')
VAO_ED = np.min(VAO)
VAO_ES = np.max(VAO)
PAO_ED = np.min(PAO)
PAO_ES = np.max(PAO)
popt, _ = curve_fit(AO_L_func, [VAO_ED,VAO_ES], [PAO_ED,PAO_ES])
print(popt)
plt.plot([VAO_ED,VAO_ES],AO_L_func([VAO_ED,VAO_ES],popt[0],popt[1]),'-')
VAO0 = -popt[1]/popt[0]
print(VAO0)
print(AO_L_func(VAO0,popt[0],popt[1]))
plt.plot([VAO0,VAO_ES],AO_L_func([VAO0,VAO_ES],popt[0],popt[1]),'-')
C = 1.0/popt[0]
plt.plot([VAO0,VAO_ES],PAAOV0_L_fun([VAO0,VAO_ES],C,VAO0),'-')
print(f"C for AO is {C}")

plt.show()
