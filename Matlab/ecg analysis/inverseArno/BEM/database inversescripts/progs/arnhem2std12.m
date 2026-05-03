% arnhem2std12.m
% function PSI=arnhem2std12(PHI,mode);
% extract 12 std lead signals from 64 nim data
% order in output: I II III VR VL VF V1.. V6
% if mode='a' : augmented
% correction 2005-07-20; sign of lead I changed; typo:
%                                                PSI(9,:)=(PHI(33,:)+PHI(35,:) corrected

function PSI=nim2std12(PHI,mode);

%        VR VL VF  V1 V2   V3 V3   V4 V5 V6 
labs=[  98 97 99  20 28   36 37   45 53 61];

fac=1; 
if nargin>1,
   if mode=='a',
       fac=1.5;
       'extremity leads augmented'
   end
end
PHI=zeromean(PHI);
PSI(  1,:)=PHI(labs(1),:)-PHI(labs(2),:);
PSI(  2,:)=PHI(labs(3),:)-PHI(labs(1),:);
PSI(  3,:)=PHI(labs(3),:)-PHI(labs(1),:);
PSI(4:6,:)=PHI(labs(1:3),:)*fac;

if norm(PSI(4:6,:),'fro')<= sqrt(eps),
    'output has WCT reference',
else,
    'output NOT refered to WCT'
end

PSI(7:8,:)=PHI(labs(4:5),:);
PSI(9,:)=mean(PHI(labs(6:7),:));
PSI(10:12,:)=PHI(labs(8:10),:);
