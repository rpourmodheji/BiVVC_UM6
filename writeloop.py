import numpy as np
import datetime

def writeloop(t,V,P,compartment):
    now = datetime.datetime.now()
    filename = "PVt_"+compartment+str(now.year)+"_" \
                                 +str(now.year)+"_"\
                                 +str(now.day)+"_"\
                                 +str(now.hour)+"_"\
                                 +str(now.minute)+"_"\
                                 +str(now.second)+".txt"
    # filename = "PVt_"+compartment+str(now.year)+"_" \
    #                              +str(now.day)+"_"+".txt"
    filename = "PVt_"+compartment+".txt"
    text_file = open(filename, 'w')
    for i in range(len(t)):
        text_file.write("%f %f %f\n" % ( t[i], V[i] , P[i]) )
    text_file.close()
