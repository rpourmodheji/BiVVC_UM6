


LV_file = open('PVt_LV.txt', 'w')
for i in range(len(t_total)):
    LV_file.write("%f %f %f\n" % ( t_total[i], V_LV_total[i] , P_LV_total[i]) )
LV_file.close()
