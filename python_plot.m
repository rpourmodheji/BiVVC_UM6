clear
clc
close all
%%
load PVt_LV.txt
load PV_experiment.txt;

% PVt_LV = PVt_LV;
t = PVt_LV(:,1);
V_LV = PVt_LV(:,2);
P_LV = PVt_LV(:,3)*0.0075;

V_LV_experiment = PV_experiment(:,2);
P_LV_experiment = PV_experiment(:,3);

load PVt_art.txt
% PVt_art = PVt_art;
t = PVt_art(:,1);
V_art = PVt_art(:,2);
P_art = PVt_art(:,3)*0.0075;

figure
hold on
plot(V_LV,P_LV)
plot(V_LV_experiment,P_LV_experiment,'*')
xlim([20 70])



figure
hold on
plot(t,V_LV)

figure
hold on
plot(t,P_LV)
plot(t,P_art)



% load PVt_LV2019_2019_1_16_11_58.txt
% PVt_LV = PVt_LV2019_2019_1_16_11_58;
% load PVt_LV.txt
% t = PVt_LV(:,1);
% V = PVt_LV(:,2);
% P = PVt_LV(:,3);
% 
% figure
% hold on
% plot(t,V)
% plot(t,P*0.0075)
% 
% figure
% hold on
% plot(V,P*0.0075)
% xlim([30 130])
% ylim([0 130])
% plot(V,P)
