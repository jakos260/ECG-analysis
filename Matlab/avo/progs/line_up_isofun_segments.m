% line_up_isofun_segments.m
% and resample contours

% script called by remesh_atri
% 20110829

% forms contiunous strings from the individual (scattered) segments of the
% contour lines X,Y,Z at level: stripes(icco)

% assumes all segments form closed loops

crit=1.e-5;

nnewn=0;
ALLNEWC=[];

nclosed=0;

for icco=1:nstripes,
    
    contlev=stripes(icco);
    
    if contlev==0,
       % force NEWC to coincide with the border segments
       for i=border_segments,
           border=LABS(ismember(LABS(:,2),i),1);
           nb=size(border,1);
           NEWC=VER(border,:);
           nclosed=nclosed+1;
           ALLNEWC=[ALLNEWC; NEWC ones(nb,1)*[icco nclosed]];
           nnewn=nnewn+nb-1;
           figure(1)
           plot3(NEWC(:,1),NEWC(:,2),NEWC(:,3),'r','linewidth',1.5)
           ['closed contour: ' num2str(nclosed) '  ilev: ' num2str(icco) ' funval: ' num2str(contlev)]
           figure(2)
           plot3(NEWC(:,1),NEWC(:,2),NEWC(:,3),'r+-','linewidth',1.5)
           %pause
        end
    else 
        XC=X(:,C(1,:)==contlev);
        YC=Y(:,C(1,:)==contlev);
        ZC=Z(:,C(1,:)==contlev);
        LIST=[[XC(1,:)' ;XC(2,:)'] [YC(1,:)' ;YC(2,:)'] [ZC(1,:)' ;ZC(2,:)']];
        nc=size(XC,2);
        col1=(1:2*nc)';
        LIST=[col1 LIST];
        bead=1;
        string=bead;
        used=1;
        STRINGS=[];
        for j=1:nc,
           k=find(abs(LIST(:,2)-LIST(bead,2))<=crit &  abs(LIST(:,3)-LIST(bead,3))<=crit & abs(LIST(:,4)-LIST(bead,4))<=crit);
           used=[used; k];
           k(k==bead)=[];
           if k>nc, 
               bead=k-nc;
           else, 
               bead=k+nc; 
           end 
           string=[string;bead];
        
           if string(1)==string(end),
               % string is closed
               ns=size(string,1)
               if ns>nskip,
                  nclosed=nclosed+1;
                  STRING=LIST(string,2:4);
                  length=sum(norm3d(STRING(2:ns,:)-STRING(1:ns-1,:))) 
                  nrec=round(length/6);  % delta=6 mm
                  nrec=max(6,nrec);
               
                  figure(2)
                  plot3(STRING(:,1),STRING(:,2),STRING(:,3),'b+','linewidth',1)
                  axis equal
                  axis square
                  hold on
                  NEWC=resample_contour(STRING,nrec,2);
                  nrecon=size(NEWC,1)
                  ALLNEWC=[ALLNEWC; NEWC ones(nrecon,1)*[icco nclosed]];
                  nnewn=nnewn+nrecon
             
                  figure(1)
                  plot3(NEWC(:,1),NEWC(:,2),NEWC(:,3),'r','linewidth',1)
                 ['closed contour: ' num2str(nclosed) '   icont: ' num2str(icco) '  funval: ' num2str(contlev) ]
                  %pause
               end
        
               if size(used,1)<2*nc, 
                  % search new starting bead
                  starters=col1;
                  starters(used)=[];
                  bead=starters(1);
                  string=bead;
               end  
           end
         end
     end
  end

 
                          
%                nh=10;
%                nh=min(10,round(ns/2)-1)              
%                if icco==12, nh=1, end
%                if icco==11, nh=3, end
%                NEWC=fourier_contours(STRING,nh,2,nrec,1);    
    

