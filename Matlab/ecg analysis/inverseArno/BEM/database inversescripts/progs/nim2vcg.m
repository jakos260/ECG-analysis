% nim2vcg.m
% function VCG=nim2vcg(PHI)
% weights according to FRANK;; axes: frontal; left ; superior
% approximate Frank electrodes in old nim lead syst 

%              A  C    E   I   M   H    F
%franklabs1=[ 53 33   26   7  62  61   64];  %for nim64 signals 

%franklabs2=[ 53 33  147   7  76 135  204];  %for gutor300 nodes 

% 20050721


function VCG=nim2vcg(PHI)
franklabs=[53       33         26           64        61       7          62];
%           A        C          E            F         H       I           M  
%(alphabetical)

FRANK =[ -0.1330    0.2310    0.3740         0         0       0.2640    -0.7360;
          0.6100    0.1710         0         0         0      -0.7810         0;
          0         0              0   -0.6550         1.      0         -0.3450];

VCG=FRANK*PHI(franklabs,:);     


     
