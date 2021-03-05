%gen data
close all
n=100;
x=linspace(0,10,n);
y1=normrnd(4,1,n,1);
y2=normrnd(6,1,n,1);

m1=measDist(x,y1,[],'normal');
m2=measDist(x,y2,[],'normal');
obj=dvDists(m1,m2);
obj.plot_measurements();
obj.plot_conditionals();
obj.plot();
obj.plot_ROCs();
