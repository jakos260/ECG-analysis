% baselcor.m
% potentials in interval of width: range around tb and te are forced to be zero

% SCRIPT version of the FUNCTION baselinecor
% specify tb, te, range
% use range>1 if noisy data are involved

range=1;

if range > 1,
   fbeg=mean(PHI(:,tb:tb+range-1)');
   fend=mean(PHI(:,te-range+1:te)');
else,
   fbeg=PHI(:,tb)';
   fend=PHI(:,te)';
end

   slope=(fend'-fbeg')/(te-tb);
   xt=1:nt;
   PHI=PHI-fbeg'*ones(1,nt)-slope*(xt-tb);
   
   rrms=rms(PHI);
   setplots

