% ams2std12.m
% function PSI=ams2std12(PHI,mode);
% extract 12 std lead signals from 64 ams data
% if mode='a' : produces augmented extremity leads, else: VR, VL and VF
% output I, II, III VR, VL VF, V1....  V6
% assumes WCT reference

% 20090421 A. vaqn Oosterom

function PSI=ams2st12(PHI,mode);
labs=[12 18 25 31 40 45 63 64];
%     V1 V2 V3 V4 V5 V6 VL VR      !!!!!!!! NB: VL VR


PHI=PHI(labs,:);
PSI(1,:)=PHI(7,:)-PHI(8,:);         % I
PSI(6,:)=-(PHI(7,:)+PHI(8,:));      % produces VF % assumng WCT
PSI(2,:)=PSI(6,:)-PHI(8,:);         % II
PSI(3,:)=PSI(6,:)-PHI(7,:);         % III 
PSI(7:12,:)=PHI(1:6,:);             % V1 ....V6
PSI(4:5,:)=PHI([8 7],:);

if nargin>1,
   if mode=='a',
      PSI(4:6,:)=PHI(4:6,:)*3/2;    % apply augmentation
   end
end 
