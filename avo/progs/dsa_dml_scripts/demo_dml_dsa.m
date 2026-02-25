% demo_dml_dsa
% A. van Oosterom;

clear

% demo of  dml and dsa  
% specify vertices of two triangles

% NBNB: all results expressed in units(4*pi) and with unit conductivity of
% an (infinite) medium

% parameters to play with: a; b; ITRI; P1; P2; nobs

a=1; b=1; 

VER=[-a/2  -b/2 0;
     -a/2   b/2 0;
      a/2   b/2 0;
      a/2  -b/2 0;];

figure(1)
clf
nver=size(VER,1);

for i=1:nver,
    plot3(VER(i,1),VER(i,2),VER(i,3),'r*')
    text(VER(i,1),VER(i,2),VER(i,3),num2str(i),'fontsize',12)
    hold on
end
grid on

plot3([0 1],[0 0], [0 0],'k' ,'linewidth',1.5)
plot3([0 0],[0 1], [0 0],'k' ,'linewidth',1.5)
plot3([0 0],[0 0], [0 1],'k' ,'linewidth',1.5)
xlabel('x-axis'); ylabel('y-axis'); zlabel('z-axis')


ITRI=[1 2 3;
      3 4 1];
ntri=size(ITRI,1);

% monolayer distributed over the two triangles, elsewhere zero
s=ones(4,1); % monolayer strength at vertices 1:4

% s=[0 0 1 0]'
% s=[0 1 1 0]'

for i=1:ntri,
    quatro=[ITRI(i,:) ITRI(i,1)];
    plot3(VER(quatro,1),VER(quatro,2), VER(quatro,3))
end

P1=2*[ 0.25  0.5  0.5];
P2=2*[-0.25 -0.5 -0.5];

plot3(P1(1),P1(2),P1(3),'*g')
text(P1(1),P1(2),P1(3),num2str(1),'fontsize',12)

plot3(P2(1),P2(2),P2(3),'*g') 
text(P2(1),P2(2),P2(3),num2str(2),'fontsize',12)

plot3([P1(1) P2(1)],[P1(2) P2(2)], [P1(3) P2(3)],':g') % line segment carrying the observation points 
% (((note that the line in this example is not normal to either of the
% source triangles

view(22,16)
axis equal
pause

nobs=401;
lambda=(0:nobs-1)/(nobs-1);
OBS=(1-lambda)'*P1+lambda'*P2;

phi=zeros(nobs,1);
psi=zeros(nobs,1);

small=1.e-4;

for i=1:nobs,
    
    obs=OBS(i,:);
    
    PSI=dml(VER,ITRI,obs); % size(ntri,3)
    PHI=dsa(VER,ITRI,obs,small)'; % size(ntri,3)
    
    for j=1:ntri,
        tris=ITRI(j,:);
        psi(i)=psi(i)+PSI(j,:)*s(tris);
        phi(i)=phi(i)+PHI(j,:)*s(tris);
    end
end


figure(2)
clf
plot(lambda,psi) % mono layer
xlabel('lambda(obs) along line connecting P1 to P2')
ylabel('potential')
title(['source strengths at vertices:  ' num2str(s')])
hold on
plot(lambda,phi,'r')  % double layer

% NOTE: obs crosses the source triangles triangles for lambda=0.5
% where phi is continuous, while its derivative is not

% locations of P1,P2, the vertices of the triangles,  as well as of the source strenghths at the vertices may be changed to
% study the nature of the field generated








