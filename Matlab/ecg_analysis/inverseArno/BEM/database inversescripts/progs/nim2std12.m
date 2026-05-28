% nim2std12.m
% function PSI=nim2std12(PHI,mode);
% extract  std-12 lead signals from nim data
% order output: I II III VR VL VF V1.. V6
% if mode='a' : augmented
% 2012_12_06


function PSI=nim2std12(PHI,mode)

if size(PHI,1)==65,
    ilf=65;
    %'in nim2std12 i_lf==65'
else, 
    %'in nim2std12 i_lf==56'
    ilf=56;  
end

%        VR VL VF  V1 V2   V3 V3   V4 V5 V6 
% labs=[  1  2 56  19 26   33 34   41 48 54];


PSI(  1,:)=PHI( 2,:)-PHI(1,:); % lead I
PSI(  2,:)=PHI(ilf,:)-PHI(1,:); % lead II
PSI(  3,:)=PHI(ilf,:)-PHI(2,:); % lead III

PSI(4:6,:)=PHI([1 2 ilf],:);

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







