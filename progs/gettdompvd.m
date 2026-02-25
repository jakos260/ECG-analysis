% gettdom
% function [tdom,ttm]=gettdom(PHI,jpoint)
% 2011-06-25
% Estimate dominant T wave
% from a weighted mean of T waves contained in PHI
% jpoint is a time marker after which depolarization is
% assumed to be completed
% ttm is timing apex Tdom
% sigmoid curve used in estimate during phase coinciding with QRS


function [tdom,ttm]=gettdompvd(PHI,jpoint)


PSI = zeros(size(PHI));
PSI(:,jpoint:size(PHI,2))=PHI(:,jpoint:size(PHI,2));

% estimate of tdom from weighted leads
tdom = rms(PSI);%
% tdom = ones(1,nt)*PSI'*PSI;
[~, ttm]=max(tdom);
 
slope = 300*tdom(jpoint) / jpoint;
t = -(jpoint-2)*2/3:(jpoint-2)*1/3;
S= tdom(jpoint) ./(1+exp(-slope*t));
S = S-min(S);
S = S/max(S) * tdom(jpoint);
tdom(1:length(S)) = S;
plot(tdom)

ap = 1 - cumsum(tdom);
ap = ap-min(ap);
ap = ap/max(ap);
plot(ap)

 

    
   
