% line_up_contour_segments.m
% 
% 20130723

% form a contiunous string from the individual (scattered) segments of the
% contour lines X,Y,Z
% this version assumes that contours (stripes) form closed loops
% containing more than 2 line segments

crit=1.e-5;

STRINGS=[];

for icco=1:nstripes,  
    contlev=stripes(icco);
    XC=X(:,C(1,:)==contlev);
    YC=Y(:,C(1,:)==contlev);
    ZC=Z(:,C(1,:)==contlev);
    nc=size(XC,2);
    
    LIST=[[XC(1,:)' ;XC(2,:)'] [YC(1,:)' ;YC(2,:)'] [ZC(1,:)' ;ZC(2,:)']];
    col1=(1:2*nc)';
    LIST=[col1 LIST];
    
    bead=1;
    string=bead;
    for j=1:nc,
        k=find(abs(LIST(:,2)-LIST(bead,2))<=crit &  abs(LIST(:,3)-LIST(bead,3))<=crit  ...
            & abs(LIST(:,4)-LIST(bead,4))<=crit);                  
        k(k==bead)=[];
        if k>nc, 
            bead=k-nc;
        else, 
            bead=k+nc; 
        end 
        string=[string;bead];  
    end
     
    STRING=LIST(string,2:4);
    ns=size(STRING,1);
    
    while norm3d(STRING(1,:)- STRING(ns,:))>crit;
          STRING(ns,:)=[];
          ns=ns-1;
    end
    STRINGS=[STRINGS;STRING ones(ns,1)*icco];
end       
   

        

 
