% triplo.m
% basis of triplot; calling trigrid

er=0;
alpha=.5;
p=.0;
izin=1;
axis([-1 1 -1 1]);
axis('equal')
axis('square')
pp=[er alpha p izin];
subplot(1,2,1)
trigrid(VER,ITRI,pp)
subplot(1,2,1)
trigrid(VER,ITRI,pp)
subplot(1,2,2)
er=1;
pp=[er alpha p izin];
trigrid(VER,ITRI,pp)
