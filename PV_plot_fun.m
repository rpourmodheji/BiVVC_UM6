function [] =  PV_plot_fun(filename)
load filename
t = filename(:,1);
V = filename(:,2);
P = filename(:,3);
plot(t,P*0.0075)
end