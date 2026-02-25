% rms.m
% function produces rowvector: r=rms(PHI,mode); 
% mode==1: rms of columns of PHI  (default)
% mode==2: rms after shift to zeromean of columns of PHI (rmszm)
%          in this case: rms(t)=std(t)  i.e. the standard deviation of the lead signals as a funtion of t
% If the reference potential is NOT implied in the data, add dummy row with
%    zero entries to PHI; 
% If WCT is implied, ensure that the limb-leads, if included, are unweighted
% Use common reference leads only. i. e. exclude leads I, II, III

% 2012-09-29

function r=rms(PHI,mode);
n=size(PHI,1);

if nargin<=1, mode=1; end
if mode >2 | mode<1, mode==1; end
    
if mode==1
   r=sqrt(sum(PHI.*PHI,1)/n);
else
    Z=zeromean(PHI);
    r=sqrt(sum(Z.*Z,1)/n);
end




    


