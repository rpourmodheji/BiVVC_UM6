import os as os
import sys
from read_elastance_parameters import *
import numpy as np
import math
import yaml

Pa_to_mmHg = 0.00750062
mmHg_to_Pa = 1.0/Pa_to_mmHg
kPa_to_mmHg = Pa_to_mmHg*1000.0
mmHg_to_kPa = 1.0/kPa_to_mmHg


class Elastance(object):
    """docstring for Elastance."""

    def __init__(self, arguments):
        self.volumes = arguments["Volumes"]
        self.t = arguments["Time"]
        self.cycle_num = arguments["Cycle Number"]
        self.cycle_time = arguments["Cycle Time"]

    def elastance_AO(self):
        V = self.volumes["AO"]
        P = (1.0/C_AO)*( V - V_AO0)           #Pa
        return P

    def elastance_sa(self):
        V = self.volumes["sa"]
        P = (1.0/C_sa)*( V - V_sa0)           #Pa
        return P

    def elastance_sv(self):
        V = self.volumes["sv"]
        P = (1.0/C_sv)*( V - V_sv0)             #Pa
        return P

    def elastance_RA(self):
        V = self.volumes["RA"]
        t = self.t
        cycle_num = self.cycle_num
        cycle_time = self.cycle_time
        te = t - (cycle_num-1.0)*cycle_time
        te = te + delay_RA
        if te > cycle_time:
            te = te - cycle_time
        if te == cycle_time:
            te = 0
        e = atrial_driving_function(te,t_max_RA,tau_RA)
        V_Sbar = V - V_RAS0
        V_Dbar = V - V_RAD0
        P_esRA = E_esRA*V_Sbar
        P_edRA = A_RA * ( math.exp(B_RA*V_Dbar) - 1.0   )
        P = e*P_esRA + (1.0-e)*P_edRA
        # P = (V-V_RA0)/0.001
        return P

    def elastance_RV(self):
        V = self.volumes["RV"]
        t = self.t
        cycle_num = self.cycle_num
        cycle_time = self.cycle_time
        te = t - (cycle_num-1.0)*cycle_time
        te = te + delay_RV
        if te > cycle_time:
            te = te - cycle_time
        if te == cycle_time:
            te = 0
        e = driving_function(te,t_max_RV,tau_RV)
        V_Sbar = V - V_RVS0
        V_Dbar = V - V_RVD0
        P_esRV = E_esRV*V_Sbar
        P_edRV = A_RV * ( math.exp(B_RV*V_Dbar) - 1.0   )
        P = e*P_esRV + (1.0-e)*P_edRV
        return P

    def elastance_PA(self):
        V = self.volumes["PA"]
        P = (1.0/C_PA)*( V - V_PA0)             #Pa
        return P

    def elastance_pc(self):
        V = self.volumes["pc"]
        P = (1.0/C_pc)*( V - V_pc0)             #Pa
        return P

    def elastance_pv(self):
        V = self.volumes["pv"]
        P = (1.0/C_pv)*( V - V_pv0)             #Pa
        return P

    def elastance_LA(self):
        V = self.volumes["LA"]
        t = self.t
        cycle_num = self.cycle_num
        cycle_time = self.cycle_time
        te = t - (cycle_num-1.0)*cycle_time
        te = te + delay_LA
        if te > cycle_time:
            te = te - cycle_time
        if te == cycle_time:
            te = 0
        e = atrial_driving_function(te,t_max_LA,tau_LA)
        V_Sbar = V - V_LAS0
        V_Dbar = V - V_LAD0
        P_esLA = E_esLA*V_Sbar
        P_edLA = A_LA * ( math.exp(B_LA*V_Dbar) - 1.0   )
        P = e*P_esLA + (1.0-e)*P_edLA
        # P = (V-V_LA0)/0.001
        return P

    def elastance_LV(self):
        V = self.volumes["LV"]
        t = self.t
        cycle_num = self.cycle_num
        cycle_time = self.cycle_time
        te = t - (cycle_num-1.0)*cycle_time
        te = te + delay_LV
        if te > cycle_time:
            te = te - cycle_time
        if te == cycle_time:
            te = 0
        e = driving_function(te,t_max_LV,tau_LV)
        V_Sbar = V - V_LVS0
        V_Dbar = V - V_LVD0
        P_esLV = E_esLV*V_Sbar
        P_edLV = A_LV * ( math.exp(B_LV*V_Dbar) - 1.0   )
        P = e*P_esLV + (1.0-e)*P_edLV
        return P

    def pressure_pack_Pa(self):
        P_AO = self.elastance_AO()
        P_sa = self.elastance_sa()
        P_sv = self.elastance_sv()
        P_RA = self.elastance_RA()
        P_RV = self.elastance_RV()
        P_PA = self.elastance_PA()
        P_pc = self.elastance_pc()
        P_pv = self.elastance_pv()
        P_LA = self.elastance_LA()
        P_LV = self.elastance_LV()
        P = {   "AO" : P_AO,
                "sa" : P_sa,
                "sv" : P_sv,
                "RA" : P_RA,
                "RV" : P_RV,
                "PA" : P_PA,
                "pc" : P_pc,
                "pv" : P_pv,
                "LA" : P_LA,
                "LV" : P_LV
                }
        # self.P = P
        return P

    def pressure_pack_mmHg(self):
        P_AO = self.elastance_AO()*Pa_to_mmHg
        P_sa = self.elastance_sa()*Pa_to_mmHg
        P_sv = self.elastance_sv()*Pa_to_mmHg
        P_RA = self.elastance_RA()*Pa_to_mmHg
        P_RV = self.elastance_RV()*Pa_to_mmHg
        P_PA = self.elastance_PA()*Pa_to_mmHg
        P_pc = self.elastance_pc()*Pa_to_mmHg
        P_pv = self.elastance_pv()*Pa_to_mmHg
        P_LA = self.elastance_LA()*Pa_to_mmHg
        P_LV = self.elastance_LV()*Pa_to_mmHg
        P_mmHg = {   "AO" : P_AO,
                "sa" : P_sa,
                "sv" : P_sv,
                "RA" : P_RA,
                "RV" : P_RV,
                "PA" : P_PA,
                "pc" : P_pc,
                "pv" : P_pv,
                "LA" : P_LA,
                "LV" : P_LV
                }
        return P_mmHg

def driving_function(t, t_max, tau):
    pi = math.pi
    if t <= 1.5*t_max:
        e = (math.sin(0.5*pi*t/t_max))**2.0   # 0.5 * ( math.sin(math.pi*t/t_max-math.pi/2.0)+ 1.0)
    else:
        e = 0.5 * math.exp( - (t-1.5*t_max)/tau  )
    return e

def atrial_driving_function(t, t_max, tau):
    pi = math.pi
    if t <= 1.5*t_max:
        e = (math.sin(0.5*pi*t/t_max))**2.0   # 0.5 * ( math.sin(math.pi*t/t_max-math.pi/2.0)+ 1.0)
    else:
        A = 1.0
        e = 0.5 * math.exp( - (t-1.5*t_max)/tau  )
    return e
