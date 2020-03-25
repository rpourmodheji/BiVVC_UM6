import numpy as np
from dVdt_system import DVdt_system as system
from elastance import Elastance
from initiations import *
from getvolumes import Getvolumes as getvolumes
from writeloop import writeloop
import yaml
from stringsall import *
import datetime

import os as os
import sys
import numpy as np
import yaml
from scipy.interpolate import interp1d
from scipy import interpolate
import matplotlib.pyplot as plt
from matplotlib import pylab as plt
from matplotlib import rc
import matplotlib
import scipy
import pandas as pd
import pdb
import copy
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

RVVmin = np.min(RVV)
RVVmax = np.max(RVV)
RVPmin = np.min(RVP)
RVPmax = np.max(RVP)


LVVmin = np.min(LVV)
LVVmax = np.max(LVV)
LVPmin = np.min(LVP)
LVPmax = np.max(LVP)

lowerxlim = 0.7
upperxlim = 1.2
# print(cardiac_data)
# plt.figure(1)
# plt.plot( RVV , RVP ,'-o')
# plt.xlim((RVVmin*lowerxlim,RVVmax*upperxlim))
# plt.ylim((RVPmin*lowerxlim,RVPmax*upperxlim))
#
# plt.figure(2)
# plt.plot( LVV , LVP ,'-o')
# plt.xlim((LVVmin*lowerxlim,LVVmax*upperxlim))
# plt.ylim((LVPmin*lowerxlim,LVPmax*upperxlim))

# plt.show()
#
# pdb.set_trace()

AO_npzfile = np.load("UM6_AO_DVP.npz")
PA_npzfile = np.load("UM6_PA_DVP.npz")

DAO = AO_npzfile["Diameter"]
VAO = AO_npzfile["Volume"]
PAO = AO_npzfile["Pressure"]
DPA = PA_npzfile["Diameter"]
VPA = PA_npzfile["Volume"]
PPA = PA_npzfile["Pressure"]



print ("                                                                ")
print ("                      ** The Circulatory **                     ")
print ("                                                                ")

for i in range(num_cycle): #i starts form zero!
    tt = trial_time(i,cycle_time)
    for j in range(len(tt)):
        print(dashline)
        t_trial = tt[j]

        print("t_trial",t_trial)
        arguments = { "Time" : t_trial,
                      "Volumes" : V,
                      "Cycle Number" : (i+1),
                      "Cycle Time" : cycle_time
                      }
        # print("Volumes LV " , arguments["Volumes"])
        print(dashline)
        A = system(arguments)
        dVdt = A.dVdt()
        # print(" LV " , dVdt)

        B = getvolumes(arguments)
        V_total = B.gettotalvolume()
        print(f'The total volume is {V_total} ml')


        dVdt_sa = dVdt["dVdt_sa"]
        dVdt_sv = dVdt["dVdt_sv"]
        dVdt_RA = dVdt["dVdt_RA"]
        dVdt_RV = dVdt["dVdt_RV"]
        dVdt_PA = dVdt["dVdt_PA"]
        dVdt_pc = dVdt["dVdt_pc"]
        dVdt_pv = dVdt["dVdt_pv"]
        dVdt_LA = dVdt["dVdt_LA"]
        dVdt_LV = dVdt["dVdt_LV"]
        dVdt_AO = dVdt["dVdt_AO"]


        ij = i*num_inc_per_cycle + j

        t[ij] = t_trial

        V_sa_trial = V["sa"] + dVdt_sa * dt
        V_sv_trial = V["sv"] + dVdt_sv * dt
        V_RA_trial = V["RA"] + dVdt_RA * dt
        V_RV_trial = V["RV"] + dVdt_RV * dt
        V_PA_trial = V["PA"] + dVdt_PA * dt
        V_pc_trial = V["pc"] + dVdt_pc * dt
        V_pv_trial = V["pv"] + dVdt_pv * dt
        V_LA_trial = V["LA"] + dVdt_LA * dt
        V_LV_trial = V["LV"] + dVdt_LV * dt
        V_AO_trial = V["AO"] + dVdt_AO * dt




        V = {   "AO" : V_AO_trial,
                "sa" : V_sa_trial,
                "sv" : V_sv_trial,
                "RA" : V_RA_trial,
                "RV" : V_RV_trial,
                "PA" : V_PA_trial,
                "pc" : V_pc_trial,
                "pv" : V_pv_trial,
                "LA" : V_LA_trial,
                "LV" : V_LV_trial
                 }


        PVs = Elastance(arguments)

        ###### Right Heart
        V_PA[ij] = V_PA_trial #V["PA"] + dVdt_PA * dt
        V_pc[ij] = V_pc_trial #V["PA"] + dVdt_PA * dt
        V_pv[ij] = V_pv_trial #V["pv"] + dVdt_pv * dt
        V_RA[ij] = V_RA_trial #V["RA"] + dVdt_RA * dt
        V_RV[ij] = V_RV_trial #V["RV"] + dVdt_RV * dt
        P_PA[ij] = PVs.elastance_PA()
        P_pc[ij] = PVs.elastance_pc()
        P_pv[ij] = PVs.elastance_pv()
        P_RA[ij] = PVs.elastance_RA()
        P_RV[ij] = PVs.elastance_RV()

        ###### Left Heart
        V_AO[ij] = V_AO_trial #V["AO"] + dVdt_AO * dt
        V_sa[ij] = V_sa_trial #V["sa"] + dVdt_sa * dt
        V_sv[ij] = V_sv_trial #V["sv"] + dVdt_sv * dt
        V_LA[ij] = V_LA_trial #V["LA"] + dVdt_LA * dt
        V_LV[ij] = V_LV_trial #V["LV"] + dVdt_LV * dt

        P_AO[ij] = PVs.elastance_AO()
        P_sa[ij] = PVs.elastance_sa()
        P_sv[ij] = PVs.elastance_sv()
        P_LA[ij] = PVs.elastance_LA()
        P_LV[ij] = PVs.elastance_LV()



print ("                                                                ")
print ("                                                                ")







################ Waveform #####################
# volume Waveform
filtered_bool = t>(cycle_time*18)

plt.figure("Waveform Volume")
plt.plot(t[filtered_bool],V_AO[filtered_bool],label = "AO")
plt.plot(t[filtered_bool],V_sa[filtered_bool],label = "sa")
plt.plot(t[filtered_bool],V_sv[filtered_bool],label = "sv")
plt.plot(t[filtered_bool],V_RA[filtered_bool],label = "RA")
plt.plot(t[filtered_bool],V_RV[filtered_bool],label = "RV")
plt.plot(t[filtered_bool],V_PA[filtered_bool],label = "PA")
plt.plot(t[filtered_bool],V_pc[filtered_bool],label = "pc")
plt.plot(t[filtered_bool],V_pv[filtered_bool],label = "pv")
plt.plot(t[filtered_bool],V_LA[filtered_bool],label = "LA")
plt.plot(t[filtered_bool],V_LV[filtered_bool],label = "LV")
plt.legend()


plt.figure("Waveform Pressure Right")
plt.plot(t[filtered_bool],P_RA[filtered_bool],label = "RA")
plt.plot(t[filtered_bool],P_RV[filtered_bool],label = "RV")
plt.plot(t[filtered_bool],P_PA[filtered_bool],label = "PA")
plt.plot(t[filtered_bool],P_pc[filtered_bool],label = "pc")
plt.plot(t[filtered_bool],P_pv[filtered_bool],label = "pv")
plt.legend()

plt.figure("Waveform Pressure Left")
plt.plot(t[filtered_bool],P_LA[filtered_bool],label = "LA")
plt.plot(t[filtered_bool],P_LV[filtered_bool],label = "LV")
plt.plot(t[filtered_bool],P_AO[filtered_bool],label = "AO")
plt.plot(t[filtered_bool],P_sa[filtered_bool],label = "sa")
plt.plot(t[filtered_bool],P_sv[filtered_bool],label = "sv")

plt.legend()

# pressure Waveform
plt.figure("Waveform Pressure")
plt.plot(t[filtered_bool],P_AO[filtered_bool],label = "AO")
plt.plot(t[filtered_bool],P_sa[filtered_bool],label = "sa")
plt.plot(t[filtered_bool],P_sv[filtered_bool],label = "sv")
plt.plot(t[filtered_bool],P_RA[filtered_bool],label = "RA")
plt.plot(t[filtered_bool],P_RV[filtered_bool],label = "RV")
plt.plot(t[filtered_bool],P_PA[filtered_bool],label = "PA")
plt.plot(t[filtered_bool],P_pc[filtered_bool],label = "pc")
plt.plot(t[filtered_bool],P_pv[filtered_bool],label = "pv")
plt.plot(t[filtered_bool],P_LA[filtered_bool],label = "LA")
plt.plot(t[filtered_bool],P_LV[filtered_bool],label = "LV")
plt.legend()




# plt.plot(V_LV_total,P_LV_total)
# plt.plot(LVV,LVP,'o')
# plt.plot(V_art_total,P_art_total)

plt.figure("RV PV-Loop")
plt.plot(V_RV[filtered_bool],P_RV[filtered_bool])
plt.plot(RVV,RVP,'o')

# plt.figure("RA PV-Loop")
# plt.plot(V_RA[filtered_bool],P_RA[filtered_bool])
# plt.plot(RVV,RVP,'o')


plt.figure("LV PV-Loop")
plt.plot(V_LV[filtered_bool],P_LV[filtered_bool])
plt.plot(LVV,LVP,'o')


# plt.figure("LA PV-Loop")
# plt.plot(V_LA[filtered_bool],P_LA[filtered_bool])


plt.figure("AO PV")
plt.plot(V_AO[filtered_bool],P_AO[filtered_bool])
plt.plot(VAO,PAO,'*')

# plt.figure("sa sv PV")
# plt.plot(V_sa[filtered_bool],P_sa[filtered_bool],'-')
# plt.plot(V_sv[filtered_bool],P_sv[filtered_bool],'--')


plt.figure("PA PV")
plt.plot(V_PA[filtered_bool],P_PA[filtered_bool])
# plt.plot(V_pc[filtered_bool],P_pc[filtered_bool])
plt.plot(VPA,PPA,'*')

# plt.figure("pc PV")
# plt.plot(V_pc[filtered_bool],P_pc[filtered_bool])
# plt.plot(VPA,PPA,'*')
#
# plt.figure("pv PV")
# plt.plot(V_pv[filtered_bool],P_pv[filtered_bool])

plt.show()
