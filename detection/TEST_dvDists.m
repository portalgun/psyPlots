d1=normrnd(0,1,200,1);
d2=normrnd(0,4,200,1);

x=linspace(min([d1;d2]),max([d1;d2]),100);
m1=measDist(x,d1,[],'normal');
m2=measDist(x,d2,[],'normal');
dv=dvDists(m1,m2);
dv.plot;
