% baselinecor.m
% 20070623; A. van Oosterom
% function PSI=baselinecor(PHI,tb,te,window)
% linear baselinecor
% potentials in interval of width: window, following tb and preceeding te, are forced to be zero
% for 50 Hz mains and 1 ms sampling: use window=20

function PSI=baselinecor(PHI,tb,te,window)

[nsig nt]=size(PHI);

if nargin<4, window=1; end

if nargin==1, tb=1; te=nt; end
if te>nt, te=nt; end

if window > 1,
   twe=min(tb+window-1,te);
   fbeg=mean(PHI(:,tb:twe)');
   
   twb=max(te-window+1,1);
   fend=mean(PHI(:,twb:te)');
else,
   fbeg=PHI(:,tb)';
   fend=PHI(:,te)';
end

   slope=(fend'-fbeg')/(te-tb);
   xt=1:nt;
   PSI=PHI-fbeg'*ones(1,nt) - slope*(xt-tb);
   
      





