from dVdt_system import DVdt_system as dVdt

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
i = 1


arguments = { "Time" : t,
              "Volumes" : V,
              "Cycle Number" : i,
              "Cycle Time" : cycle_time
              }

dV = dVdt(arguments)
dVs = dV.dVdt()
print (dVs)
