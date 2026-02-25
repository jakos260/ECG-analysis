% nim2std12.m
% function PSI=nim2std12(PHI,mode);
% extract  std_12 lead signals from nim data
% input 
% order output: I II III VR VL VF V1.. V6
% if mode='a' : augmented
% 2017_03_27


function PSI=nim2std12(PHI,mode)

if size(PHI,1)==65,
    i_lf=65;
    %'in nim2std12 i_lf==65'
else, 
    %'in nim2std12 i_lf==56'
    i_lf=56;  
end

PSI(  1,:)=PHI(   2,:)-PHI(1,:); % lead I
PSI(  2,:)=PHI(i_lf,:)-PHI(1,:); % lead II
PSI(  3,:)=PHI(i_lf,:)-PHI(2,:); % lead III

PSI(4:6,:)=PHI([1 2 i_lf],:); % extremity leads
PSI(7:8,:)=PHI([19 26],:);

V3=PHI(33,:); % V3=mean(PHI([33 34],:)); % alternative

PSI(9,:)=V3; 

PSI(10:12,:)=PHI([41 48 54],:);

if nargin>1,
   if mode=='a',
      'refer to WCT'
      wct=mean([PHI(1:2,:);PHI(i_lf,:)] );
      PSI(4:12,:)=PSI(4:12,:)-ones(9,1)*wct;
      nt=size(PSI,2); 
      if norm((PSI(5,:)-PHI(4,:)-PSI(1,:))/nt<1.e-3);
         PSI(4:6,:)=1.5*PSI(4:6,:);
         'extremity leads augmented'
      end
       
   end
end







