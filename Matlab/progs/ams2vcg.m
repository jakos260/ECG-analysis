% ams2vcg.m
% function VCG=ams2vcg(PHI);
% extract XYZ-leads from 64 ams data

function VCG=ams2vcg(PHI);
 VCG(1,:)=PHI(61,:)*0.264 + PHI(19,:)*0.374 + PHI(31,:)*0.231;
 VCG(1,:)=VCG(1,:)        - PHI(45,:)*0.133 - PHI(58,:)*0.737;
 VCG(2,:)=PHI(45,:)*0.61  + PHI(31,:)*0.171 - PHI(61,:)*0.781;
 VCG(3,:)=PHI(10,:)       - PHI(53,:)*0.655 - PHI(58,:)*0.345;
