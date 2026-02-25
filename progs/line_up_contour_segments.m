% line_up_contour_segments.m
%
% 20151118

% form a contiunous string from the individual (scattered) segments of the
% contour lines X,Y,Z
% this version requires that contours (stripes) form closed loops
% and contains more than 2 line segments

crit=1.e-5;
STRINGS=[];


for icco=1:nstripes,
    
    cont_val=stripes(icco);
    
    XC=X(:,C(1,:)==cont_val);
    YC=Y(:,C(1,:)==cont_val);
    ZC=Z(:,C(1,:)==cont_val);
    
    nc=size(XC,2);
    col1=(1:2*nc)';
    
    LIST=[col1 [XC(1,:)' ;XC(2,:)'] [YC(1,:)' ;YC(2,:)'] [ZC(1,:)' ;ZC(2,:)']];
   
    bead=1;
    
    string=bead;
    strings=[];
    nsubs=1;
    
    remo=ones(2*nc,1);
    
    while sum(remo)>0,
        
        % find contours at selected cont_val
        
        while ~isempty(bead),
            k=find(abs(LIST(:,2)-LIST(bead,2))<=crit &  abs(LIST(:,3)-LIST(bead,3))<=crit  ...
                & abs(LIST(:,4)-LIST(bead,4))<=crit);
                k(k==bead)=[];
                
            if k>nc,
                bead=k-nc;
            else,
                bead=k+nc;
            end
            remo([k bead])=zeros(2,1); 
            string=[string;bead];
            if size(string,1)>2 & string(1)==string(end), break,end     
        end
       
        ns=size(string,1);
        STRING=[LIST(string,1:4) ones(ns,1)*[nsubs icco]]';
        STRINGS=[STRINGS; STRING'];
        
        if sum(remo)==0, break,end
        k=find(remo>0);
        bead=k(1);
        string=bead;
        nsubs=nsubs+1;        
    end
  
end





