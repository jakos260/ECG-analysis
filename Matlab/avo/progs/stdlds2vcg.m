% stdlds2vcg.m
%  VCG=stdlds2vcg(PHI)
% compute VCG from weighted standard lead potentials: 
% input PHI: 9 'unipolar' signals of the 12-leads 
% VR   VL  VF   V1  V2   V3    V4   V5  V6
% see leadconversions
% 20050721
function VCG=stdlds2vcg(PHI)


% testing:
%PHI=loadmat('../gu/phi_central_ventr_dip.mat'); % gu
%PHI=nim2std12(PHI);

  
%'Uijen et al. std2vcg; as published; leads x;y;z'  
%     VR     VL             V1     V2     V3     V4     V5     V6
% M=  [ -0.52   0.82         -0.01   0.04   0.05   0.05    0.07  0.37;
%      -1.53  -1.09          0.03  -0.02  -0.02   0.03   -0.07  0.08;
%       0.43  -0.01         -0.26  -0.28  -0.14   0.04   -0.15  0.34;];

% converted to VR VL VF data ; axes: anterior; left; superior;
% result independent of (reference electrode)

M=[   -0.4400         0   -0.0100    0.2600    0.2800    0.1400   -0.0400    0.1500   -0.3400;
      -0.8100    0.5300   -0.2900   -0.0100    0.0400    0.0500    0.0500    0.0700    0.3700;
       0.6667    0.2267   -0.8633   -0.0300    0.0200    0.0200   -0.0300    0.0700   -0.0800;];

VCG=M*PHI;

