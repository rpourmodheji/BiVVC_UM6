import numpy as np

#857
cycle_time = 857.0
num_cycle = 20                      #number of cycle is integer
total_time = num_cycle*cycle_time
dt = 1.0
t = np.arange(0.0,cycle_time,dt)
num_inc_per_cycle = len(t)

# 3081
total_volume = 3081.0

#Left
left = total_volume*0.65
# V_LA = 0.006*left
# V_LV = 0.01*left
# V_AO = 0.03*left
# V_sa = 0.45*left
# V_sv = 0.5*left
V_LA = 0.03*left
V_LV = 0.0*left
V_AO = 0.17*left
V_sa = 0.4*left
V_sv = 0.4*left

##Right
right = total_volume*0.35
V_RA = 0.4*right
V_RV = 0.0*right
V_PA = 0.1*right
V_pc = 0.2*right
V_pv = 0.3*right



V = {   "AO" : V_AO,
        "sa" : V_sa,
        "sv" : V_sv,
        "RA" : V_RA,
        "RV" : V_RV,
        "PA" : V_PA,
        "pc" : V_pc,
        "pv" : V_pv,
        "LA" : V_LA,
        "LV" : V_LV
         } #Initial Conditions They will be settled down after cycles

num_all = num_inc_per_cycle*num_cycle
# print(num_all)




t = np.zeros(  num_all  )


#### Right
V_PA = np.zeros(  num_all  )
V_pc = np.zeros(  num_all  )
V_pv = np.zeros(  num_all  )
V_RA = np.zeros(  num_all  )
V_RV = np.zeros(  num_all  )

P_PA = np.zeros(  num_all  )
P_pc = np.zeros(  num_all  )
P_pv = np.zeros(  num_all  )
P_RA = np.zeros(  num_all  )
P_RV = np.zeros(  num_all  )


##### Left
V_AO = np.zeros(  num_all  )
V_sa = np.zeros(  num_all  )
V_sv = np.zeros(  num_all  )
V_LA = np.zeros(  num_all  )
V_LV = np.zeros(  num_all  )

P_AO = np.zeros(  num_all  )
P_sa = np.zeros(  num_all  )
P_sv = np.zeros(  num_all  )
P_LA = np.zeros(  num_all  )
P_LV = np.zeros(  num_all  )




def trial_time(cycle_num,cycle_time): #This in fact should be cycle_num minus one
    # print("dt = " , dt)
    ti = cycle_num*cycle_time+dt
    te = ti + cycle_time
    t  = np.arange(ti,te,dt)
    return t
