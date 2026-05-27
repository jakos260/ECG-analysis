% naf2std12.m

function PHI=naf2std12(PHI);
% create pure, i.e. non-augmented,  extremity leads
PHI(4:6,:)=PHI(4:6,:)*2/3;

% derive extremity lead VF from VR and VL
PHI(6,:)=-(PHI(4,:)+PHI(5,:));

% derive "bipolar" leads I II III

PHI(1,:)=PHI(5,:)-PHI(4,:);
PHI(2,:)=PHI(6,:)-PHI(4,:);
PHI(3,:)=PHI(6,:)-PHI(5,:);

