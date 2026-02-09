load rawdata_for_miniex
r=2;
ts=1;
t1=3;
t2=1000;
xpos=[1 2];
lb=0;
ub=1;
norm=1;
opti="fmin";
os=0;
[cpn,J,u,x]=cpnTens.timeseries2CPN1approximation(d,t1,t2,xpos,ts,r,os,lb,ub,"fmin",norm)
%cpn.x0=x(1,:);
x0=x(1,:);
us=[(0:ts:length(u)-1)',u];

si = simMTI(mss(CPN1(cpn.U,cpn.phi)));
[xx,yy] = simulate(si,us(:,2:end)',us(end,1),x0,true);
si.tsim = us(:,1)';
%si=simulateMTI(cpn,us);
%simDiscMTI(si,us(end,1)); 
%figure
plot(si)
       %plotSim(si)