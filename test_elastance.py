from elastance import Elastance

import os
import sys

import numpy as np
import math
import yaml


V = {   "AO" : 60.0,
        "sa" : 60.0,
        "sv" : 60.0,
        "RA" : 30.0,
        "RV" : 30.0,
        "PA" : 30.0,
        "pv" : 20.0,
        "LA" : 60.0,
        "LV" : 60.0
         }
t = 2.0
cycle_time = 10.0
i = 1.0


arguments = { "Time" : t,
              "Volumes" : V,
              "Cycle Number" : i,
              "Cycle Time" : cycle_time
              }

pressure = Elastance(arguments)
pressure_mmHg = pressure.pressure_pack_mmHg()
pressure_Pa = pressure.pressure_pack_Pa()

print(pressure_mmHg)
print(pressure_Pa)
