import yaml

filename = "elastance_parameters.yaml"
parmfiles = open(filename,'r')
params = yaml.safe_load(parmfiles)

################## Aorta ####################
C_AO = params["C_AO"]
V_AO0 = params["V_AO0"]

################### Systemic Artery ##################
C_sa = params["C_sa"]
V_sa0 = params["V_sa0"]

################### Systemic Veins ##################
C_sv = params["C_sv"]
V_sv0 = params["V_sv0"]

################### Right Atrium ##################
E_esRA = params["E_esRA"]
t_max_RA = params["t_max_RA"]
tau_RA = params["tau_RA"]
A_RA = params["A_RA"]
B_RA = params["B_RA"]
V_RAS0 = params["V_RAS0"]
V_RAD0 = params["V_RAD0"]
delay_RA = params ["delay_RA"]


################### Right Ventricle ##################
E_esRV = params["E_esRV"]
t_max_RV = params["t_max_RV"]
tau_RV = params["tau_RV"]
A_RV = params["A_RV"]
B_RV = params["B_RV"]
V_RVS0 = params["V_RVS0"]
V_RVD0 = params["V_RVD0"]
delay_RV = params["delay_RV"]


################### Pulmonary Artery ##################
C_PA = params["C_PA"]
V_PA0 = params["V_PA0"]

################### Pulmonary Artery ##################
C_pc = params["C_pc"]
V_pc0 = params["V_pc0"]


################### Pulmonary Veins ##################
C_pv = params["C_pv"]
V_pv0 = params["V_pv0"]



################### LA ######################
E_esLA = params["E_esLA"]
t_max_LA = params["t_max_LA"]
tau_LA = params["tau_LA"]
A_LA = params["A_LA"]
B_LA = params["B_LA"]
V_LAS0 = params["V_LAS0"]
V_LAD0 = params["V_LAD0"]
delay_LA = params ["delay_LA"]



################### LV ######################
E_esLV = params["E_esLV"]
t_max_LV = params["t_max_LV"]
tau_LV = params["tau_LV"]
A_LV = params["A_LV"]
B_LV = params["B_LV"]
V_LVS0 = params["V_LVS0"]
V_LVD0 = params["V_LVD0"]
delay_LV = params ["delay_LV"]
