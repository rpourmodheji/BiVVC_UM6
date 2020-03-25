dVdt_sa = dVdt["dV/dt Systemic Arteries"]
dVdt_sv = dVdt["dV/dt Systemic Veins"]
dVdt_RA = dVdt["dV/dt Right Atrium"]
dVdt_RV = dVdt["dV/dt Right Ventricle"]
dVdt_pa = dVdt["dV/dt Pulminory Arteries"]
dVdt_pv = dVdt["dV/dt Pulminory Veins"]
dVdt_LA = dVdt["dV/dt Left Atrium"]
dVdt_LV = dVdt["dV/dt Left Ventricle"]




t_data[i] = t

print(dVdt_sa)
print(dVdt_sv)
print(dVdt_RA)
print(dVdt_RV)
print(dVdt_pa)
print(dVdt_pv)
print(dVdt_LA)
print(dVdt_LV)

print(t[j])
print(V)

f.write(str(P_LV_data[i][j])+"\n")


file_LV.write("%f %f %f \n" % t_total[i] , V_LV_total[i] , P_LV_total[i] )

P_LV_data[i][j] = PVs.nonlinear_elastance_LV()

print(counter)
print(num_all)

print(pooldata_zeros)

print("The size = " ,  t_total.ndim)

t_total[ counter ] = t_trial
V_LV_total[ counter ] = V_LV
P_LV_total[ counter ] = P_LV


with open('PVt_LV.txt', 'w') as file_LV:
    for i in range(num_inc_per_cycle*num_cycle):
        file_LV.write("%f %f %f\n" % t_total[i], % V_LV_total[i], % P_LV_total[i]  )
file_LV.close()


V_LV_total = pooltotal_zeros
P_LV_total = pooltotal_zeros
V_sa_total = pooltotal_zeros
P_sa_total = pooltotal_zeros
V_RV_total = pooltotal_zeros
P_RV_total = pooltotal_zeros

flowrate = Flows(arguments)


q_av = flowrate.flowrate_av()
q_sa = flowrate.flowrate_sa()
q_sv = flowrate.flowrate_sv()
q_tv = flowrate.flowrate_tv()
q_pvv = flowrate.flowrate_pvv()
q_pa = flowrate.flowrate_pa()
q_pv = flowrate.flowrate_pv()
q_mv = flowrate.flowrate_mv()

dVdt_sa = q_av - q_sa
dVdt_sv = q_sa - q_sv
dVdt_RA = q_sv - q_tv
dVdt_RV = q_tv - q_pvv
dVdt_pa = q_pvv -  q_pa
dVdt_pv = q_pa - q_pv
dVdt_LA = q_pv - q_mv
dVdt_LV = q_mv - q_av


print(flowrate.flowrate_av())
print(flowrate.flowrate_sa())
print(flowrate.flowrate_sv())
print("*",flowrate.flowrate_tv())
print(flowrate.flowrate_pvv())
print(flowrate.flowrate_pa())
print(flowrate.flowrate_pv())
print(flowrate.flowrate_mv())


pooldata_zeros = np.zeros((num_cycle,num_inc_per_cycle))
t_data = pooldata_zeros
V_LV_data = pooldata_zeros
P_LV_data = pooldata_zeros
V_sa_data = pooldata_zeros
P_sa_data = pooldata_zeros
V_RV_data = pooldata_zeros
P_RV_data = pooldata_zeros


dash_line = '---------------------------------------------'
pin_line = '################################################'

print("The num_all = " ,  num_all)


import numpy as np

def set_times():
    cycle_time = 200.0
    num_cycle = 30                      #number of cycle is integer
    total_time = num_cycle*cycle_time
    dt = 1.0
    t = np.arange(0.0,cycle_time+dt,dt)
    num_inc_per_cycle = len(t)
    return cycle_time , num_cycle, total_time, dt, t, num_inc_per_cycle

print("V total = " , V_total)
