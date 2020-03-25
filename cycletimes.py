import numpy as np

def set_times():
    cycle_time = 200.0
    num_cycle = 30                      #number of cycle is integer
    total_time = num_cycle*cycle_time
    dt = 1.0
    t = np.arange(0.0,cycle_time+dt,dt)
    num_inc_per_cycle = len(t)
    return cycle_time , num_cycle, total_time, dt, t, num_inc_per_cycle
