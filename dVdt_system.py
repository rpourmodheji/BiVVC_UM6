import numpy as np
from flows import Flows

class DVdt_system(object):
    """docstring for DVdt_system."""

    def __init__(self, arguments):
        self.arguments = arguments

    def dVdt(self):
        flowrate = Flows(self.arguments)
        q_AO = flowrate.flowrate_AO()
        q_sa = flowrate.flowrate_sa()
        q_sv = flowrate.flowrate_sv()
        q_tr = flowrate.flowrate_tr()
        q_pv = flowrate.flowrate_pv()
        q_pc = flowrate.flowrate_pc()
        q_cp = flowrate.flowrate_cp()
        q_cv = flowrate.flowrate_cv()
        q_mv = flowrate.flowrate_mv()
        q_av = flowrate.flowrate_av()


        dVdt_sa = q_AO - q_sa
        dVdt_sv = q_sa - q_sv
        dVdt_RA = q_sv - q_tr
        dVdt_RV = q_tr - q_pv
        dVdt_PA = q_pv - q_pc
        dVdt_pc = q_pc - q_cp
        dVdt_pv = q_cp - q_cv
        dVdt_LA = q_cv - q_mv
        dVdt_LV = q_mv - q_av
        dVdt_AO = q_av - q_AO



        dVdt_total = dVdt_sa + \
                     dVdt_sv +  \
                     dVdt_RA + \
                     dVdt_RV + \
                     dVdt_PA + \
                     dVdt_pc + \
                     dVdt_pv + \
                     dVdt_LA + \
                     dVdt_LV + \
                     dVdt_AO
        # print(f'=============== dVdt_total is {dVdt_total}')


        dVdt = { "dVdt_sa" : dVdt_sa,
                 "dVdt_sv" : dVdt_sv,
                 "dVdt_RA" : dVdt_RA,
                 "dVdt_RV" : dVdt_RV,
                 "dVdt_PA" : dVdt_PA,
                 "dVdt_pc" : dVdt_pc,
                 "dVdt_pv" : dVdt_pv,
                 "dVdt_LA" : dVdt_LA,
                 "dVdt_LV" : dVdt_LV,
                 "dVdt_AO" : dVdt_AO
                 }

        return dVdt
