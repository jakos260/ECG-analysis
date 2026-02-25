
if 0
[VER,ITRI]=make_sphere(3);
[A,D]=graphdist(ITRI,VER,1);
else
    VER=GEOM.VER;
    ITRI=GEOM.ITRI;
    A=GEOM.neigh;
    D=GEOM.DIST2W;
end
dep =D(1,:)*30;
rep = 250 - 0.5*dep;
T = ones(length(VER),1) * 0:330;
SPECS =[0.001518 , -0.031537 ,0.072656];
tic
for k=1:4
for i=1:length(D)
    dep =D(i,:)*30;
    S0 =getSp(1,100,A,dep,rep,SPECS)';
end
end
toc
S0 =getSp(3,330,A,dep,rep,SPECS)';
S1 =getSp(2,330,A,dep,rep,[-SPECS(2)/2 SPECS(3)] )';
clf
plot(S0(1:40:end,:)','b')
hold on; plot(S1(1:40:end,:)','r')

 