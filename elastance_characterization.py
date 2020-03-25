import numpy as np
from scipy.optimize import curve_fit
# from dVdt_system import DVdt_system as system
# from elastance import Elastance
# from initiations import *
# from getvolumes import Getvolumes as getvolumes
# from writeloop import writeloop
# import yaml
# from stringsall import *
# import datetime
# import os as os
# import sys
# import numpy as np
# import yaml
# from scipy.interpolate import interp1d
# from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import pylab as plt
from matplotlib import rc
import matplotlib
# import scipy
import pandas as pd
# import pdb
# import copy
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


PVloopfile = "UM6_PVloops.xlsx"
cardiac_data = pd.read_excel(open(PVloopfile, 'rb'), sheet_name='Sheet1')
RVV = cardiac_data["RVV"]
RVP = cardiac_data["RVP"]*mmHg_to_Pa
LVV = cardiac_data["LVV"]
LVP = cardiac_data["LVP"]*mmHg_to_Pa


#
# print(cardiac_data)
plt.figure(1)
plt.plot(RVV,RVP,'-o')
plt.plot(LVV,LVP,'-o')

def ESPVR_func(V, E_es):
    V_0 = V[0]
    return E_es*(V-V_0)

def EDPVR_func(V,A_ed,B_ed):
    V_0 = V[0]
    return A_ed*( np.exp(B_ed*(V-V_0)) - 1  )


###### ESPVR RV #######
ES_V0_RV = 20.0
ES_VS_RV = 25.82
ES_P0_RV = 0.0
ES_PS_RV = 8651.0
popt, pcov = curve_fit(ESPVR_func, [ES_V0_RV,ES_VS_RV], [ES_P0_RV,ES_PS_RV])
E_es = popt[0]
print(f'E_es for RV is {E_es}')
S_x = np.linspace(ES_V0_RV,ES_V0_RV*1.5,20)
S_y = ESPVR_func(S_x,E_es)
plt.plot(S_x,S_y,'--')

###### EDPVR RV #######
ED_V0_RV = 20.0
ED_VS_RV = 60.2
ED_P0_RV = 0.0
ED_PS_RV = 1457.6
popt, pcov = curve_fit(EDPVR_func, [ED_V0_RV,ED_VS_RV], [ED_P0_RV,ED_PS_RV]  )
A_ed = popt[0]
B_ed = popt[1]
print(f'A_ed for RV is {A_ed}')
print(f'B_ed for RV is {B_ed}')
D_x = np.linspace(ED_V0_RV,ED_V0_RV*3.3,20)
D_y = EDPVR_func(D_x,A_ed,B_ed)
plt.plot(D_x,D_y,'--')


###### ESPVR LV #######
ES_V0_LV = 10.0
ES_VS_LV = 20.6
ES_P0_LV = 0.0
ES_PS_LV = 11669.8
popt, pcov = curve_fit(ESPVR_func, [ES_V0_LV,ES_VS_LV], [ES_P0_LV,ES_PS_LV])
E_es = popt[0]
print(f'E_es for LV is {E_es}')
S_x = np.linspace(ES_V0_LV,ES_V0_LV*2.7,20)
S_y = ESPVR_func(S_x,E_es)
plt.plot(S_x,S_y,'--')


###### EDPVR LV #######
ED_V0_LV = 10.0
ED_VS_LV = 62.2
ED_P0_LV = 0.0
ED_PS_LV = 2577.6
popt, pcov = curve_fit(EDPVR_func, [ED_V0_LV,ED_VS_LV], [ED_P0_LV,ED_PS_LV]  )
A_ed = popt[0]
B_ed = popt[1]
print(f'A_ed for LV is {A_ed}')
print(f'B_ed for LV is {B_ed}')
D_x = np.linspace(ED_V0_LV,ED_V0_LV*6.3,20)
D_y = EDPVR_func(D_x,A_ed,B_ed)
plt.plot(D_x,D_y,'--')



plt.show()
