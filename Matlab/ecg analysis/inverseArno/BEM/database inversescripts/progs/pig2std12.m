% pig2std12.m
% function PSI=nim2std12(PHI,mode);
% extract  std-12 lead signals from 64 nim data
% order output: I II III VR VL VF V1.. V6
% if mode='a' : augmented
% 2012_12_06


function PSI=pig2std12(PHI,mode)

%     VR VL VF  V1  V2   V3 V3   V4 V5 V6 
 l=[  1  2 56   19  26   33 34   41 48 54]; % nim2std12
 l=[  1  2 56   19  26   33 34   41 48 54];


PSI(  l(1),:)=PHI(l(2),:)-PHI(l(1),:); % lead I
PSI(  l(2),:)=PHI(l(3),:)-PHI(l(1),:); % lead II
PSI(  l(3),:)=PHI(l(3),:)-PHI(l(2),:); % lead III

PSI(l(4:6),:)=PHI(l(1:3),:);

PSI(l(7:8),:)=PHI(l(4:5),:);
PSI(9,:)=(PHI(33,:)+PHI(34,:))/2;  % V31.
PSI(10:12,:)=PHI([41 48 54],:);

if nargin>1,
   if mode=='a',
      'refered to WCT'
      wct=mean(PHI(l(1:3),:));
      PSI(l(4:12),:)=PSI(l(4:12),:)-ones(9,1)*wct;
      nt=size(PSI,2); 
      if norm((PHI(l(5),:)-PHI(l(4),:)-PSI(1,:))/nt<1.e-3);
         PSI(4:6,:)=1.5*PSI(4:6,:);
         'extremity leads augmented'
      end
       
   end
end







