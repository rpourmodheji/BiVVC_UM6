import yaml

filename = "resistance_parameters.yaml"
parmfiles = open(filename,'r')
params = yaml.safe_load(parmfiles)

################### Artery ##################
R_AO = params["R_AO"]

################### Peripheral ##############
R_per = params["R_per"]

################### Veins ###################
R_sv = params["R_sv"]

################### Tricuspid ###################
R_tr = params["R_tr"]

############ Pulmonary Valve ##########
R_pv = params["R_pv"]

############ Pulmonary Valve ##########
R_pc = params["R_pc"]

############ Capillary Pulmonary ##########
R_cp = params["R_cp"]

############ Pulmonary Veins ##########
R_cv = params["R_cv"]

################### Mitral ##################
R_mv = params["R_mv"]

################### Aortic ##################
R_av = params["R_av"]
