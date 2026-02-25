% nim2std12.m
% function PSI=nim2std12(PHI,mode);
% extract  std-12 lead signals from nim data
% order output: I II III VR VL VF V1.. V6
% if mode='a' : augmented

% 2015_06_10


function PSI=nim_rev2std12(PHI,mode)


%        VR VL VF  V1 V2   V3 V3   V4 V5 V6 
% labs=[  1  2 56  19 26   33 34   41 48 54]; NB: signal labels; 


PSI(  1,:)=PHI( 65,:)-PHI(66,:);  % lead I
PSI(  2,:)=PHI(67,:)-PHI(66,:);   % lead II
PSI(  3,:)=PHI(67,:)-PHI(65,:);   % lead III

PSI(4:6,:)=PHI([65 66 67],:);

PSI(7:8,:)=PHI([19 26],:);
PSI(9,:)=(PHI(33,:)+PHI(34,:))/2;  % V31.
PSI(10:12,:)=PHI([41 48 54],:);

if nargin>1,
   if mode=='a',
      'refered to WCT'
      wct=mean(PHI([1 2 ilf],:));
      PSI(4:12,:)=PSI(4:12,:)-ones(9,1)*wct;
      nt=size(PSI,2); 
      if norm((PHI(5,:)-PHI(4,:)-PSI(1,:))/nt<1.e-3);
         PSI(4:6,:)=1.5*PSI(4:6,:);
         'extremity leads augmented'
      end
       
   end
end







