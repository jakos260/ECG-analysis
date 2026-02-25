% schiller2_16.m
% expands Schiller  I, II, V1, V2, V3, V4, V5, V6, A1,A2, A3 A4;
% to I II III VR VL VF V1 ....V6 A1 ... A4
% assumes, i.e., demands,  WCT during recording
% A. van Oosterom; 2012_12_07


function PSI=schiller2_16(PHI)

% Schiller ful records 
% input: I, II, V1, V2, V3, V4, V5, V6, A1,A2, A3 A4;
% converted to 
% leads: I, II, III, VR, VL, VF, V1, V2, V3, V4, V5, V6, A1,A2, A3 A4;
nt=size(PHI,2);

PSI=zeros(16,nt);

PSI(1:3,:) =[PHI(1:2,:);
            PHI(2,:)-PHI(1,:)];
PSI(4:16,:)=[-(PHI(1,:) +   PHI(2,:))/3;     % VR
            (2*PHI(1,:) -   PHI(2,:))/3;     % VL
            ( -PHI(1,:) + 2*PHI(2,:))/3;     % VF 
             PHI(3:12,:);];

        
            



