from elastance import Elastance
import yaml
from read_resistance_parameters import *

class Flows(object):
    """docstring for Flows."""

    def __init__(self, arguments):
        self.arguments = arguments
        self.volumes = arguments["Volumes"]
        self.t = arguments["Time"]
        self.cycle_num = arguments["Cycle Number"]
        self.cycle_time = arguments["Cycle Time"]
        pressure = Elastance(arguments)
        print(f'pressure is {pressure.pressure_pack_Pa()}')
        self.P_AO = pressure.elastance_AO()
        self.P_sa = pressure.elastance_sa()
        self.P_sv = pressure.elastance_sv()
        self.P_RA = pressure.elastance_RA()
        self.P_RV = pressure.elastance_RV()
        self.P_PA = pressure.elastance_PA()
        self.P_pc = pressure.elastance_pc()
        self.P_pv = pressure.elastance_pv()
        self.P_LA = pressure.elastance_LA()
        self.P_LV = pressure.elastance_LV()

    def flowrate_AO(self): # AO -^-^-^- sa
        P_AO = self.P_AO
        P_sa = self.P_sa
        if 1: #P_AO >= P_sa:
            q = (P_AO-P_sa)/R_AO
        else:
            q = 0.0
        return q

    def flowrate_sa(self): # sa -^-^-^- sv
        P_sa = self.P_sa
        P_sv = self.P_sv
        if 1: #P_sa >= P_sv:
            q = (P_sa-P_sv)/R_per
        else:
            q = 0.0
        return q

    def flowrate_sv(self): # sv -^-^-^- RA
        P_sv = self.P_sv
        P_RA = self.P_RA
        if 1: #P_sv >= P_RA:
            q = (P_sv-P_RA)/R_sv
        else:
            q = 0.0
        return q

    def flowrate_tr(self): # RA -^-^-^- RV
        P_RA = self.P_RA
        P_RV = self.P_RV
        if P_RA >= P_RV:
            q = (P_RA-P_RV)/R_tr
        else:
            q = 0.0 #(P_RA-P_RV)/R_tr/8.0
        return q

    def flowrate_pv(self): # RV -^-^-^- PA
        P_RV = self.P_RV
        P_PA = self.P_PA
        if P_RV >= P_PA:
            q = (P_RV-P_PA)/R_pv
        else:
            q = 0.0
        return q

    def flowrate_pc(self): # PA -^-^-^- pc
        P_PA = self.P_PA
        P_pc = self.P_pc
        if 1: # P_PA >= P_pv:
            q = (P_PA-P_pc)/R_pc
        else:
            q = 0.0
        return q

    def flowrate_cp(self): # pc -^-^-^- pv
        P_pc = self.P_pc
        P_pv = self.P_pv
        if 1: # P_PA >= P_pv:
            q = (P_pc-P_pv)/R_cp
        else:
            q = 0.0
        return q

    def flowrate_cv(self): # pv -^-^-^- LA
        P_pv = self.P_pv
        P_LA = self.P_LA
        if 1: #P_pv >= P_LA:
            q = (P_pv-P_LA)/R_cv
        else:
            q = 0.0
        return q

    def flowrate_mv(self): # LA -^-^-^- LV
        P_LA = self.P_LA
        P_LV = self.P_LV
        if P_LA >= P_LV:
            q = (P_LA-P_LV)/R_mv
        else:
            q = 0.0
        return q

    def flowrate_av(self): # LV -^-^-^- AO
        P_LV = self.P_LV
        P_AO = self.P_AO
        if P_LV >= P_AO:
            q = (P_LV-P_AO)/R_av
        else:
            q = 0.0
        return q


    def flowrate_pack(self):
        q_AO = self.flowrate_AO()
        q_sa = self.flowrate_sa()
        q_sv = self.flowrate_sv()
        q_tr = self.flowrate_tr()
        q_pv = self.flowrate_pv()
        q_pc = self.flowrate_pc()
        q_cp = self.flowrate_cp()
        q_cv = self.flowrate_cv()
        q_mv = self.flowrate_mv()
        q_av = self.flowrate_av()
        q = {   "q_AO"  : q_AO,
                "q_sa" : q_sa,
                "q_sv" : q_sv,
                "q_tr"  : q_tr,
                "q_pv"  : q_pv,
                "q_pc"  : q_pc,
                "q_cp" : q_cp,
                "q_cv" : q_cv,
                "q_mv"  : q_mv,
                "q_av"  : q_av
                }
        return q
