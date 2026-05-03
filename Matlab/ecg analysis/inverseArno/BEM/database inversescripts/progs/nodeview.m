% nodeview.m
% function view=nodeview(VER,ITRI,j1,j2)
% view =1 if line connecting node j1 to node j2 of
% the triangulated surface specified by (VER and ITRI) lies entirely in 
% its interior; else view=0;

% function called by: tri2gra.m
% calls: solida.m

% A van Oosterom; 20010802


function view=nodeview(VER,ITRI,j1,j2)
view=1;
dim=size(ITRI);
ntri=dim(1);
node1=VER(j1,:);
node2=VER(j2,:);
a3=node1-node2;
ns=1; alphas=[]; alph=[];
alpha(1)=0;
% identify intersections in between node1 and node2

for i=1:ntri,
  itri=ITRI(i,:);
  if all(ITRI(i,:)~=j1) & all(ITRI(i,:)~=j2),
     a1=VER(ITRI(i,2),:)-VER(ITRI(i,1),:);
     a2=VER(ITRI(i,3),:)-VER(ITRI(i,1),:);
     deter=det([a1; a2; a3]);
     if abs(deter) > 1.e-9,
        % no parallel detected
        % solution via Cramer's rule
        sav=node1-VER(ITRI(i,1),:);
        lambda=det([sav;a2;a3])/deter;
        if lambda >= 0 & lambda <= 1,
           mu=det([a1; sav; a3])/deter;
           if mu >= 0 & mu <= 1 & lambda + mu <= 1,
              alpha=det([a1;a2;sav])/deter; 
              if alpha > 0 & alpha < 1, 
                 ns=ns+1; alphas(ns)=alpha;            
              end 
           end
        end
     end  
  end
end

% sort for increasing values of alpha
ns=ns+1;
alphas(ns)=1;
alphas=sort(alphas);
na=1;
alph(na)=alphas(1);

for n=2:ns,
    if alphas(n)~=alphas(n-1), na=na+1;, alph(na)=alphas(n); end
end

% alph
% test if positions halfway the intersections are interior or exterior
% if any of these is exterior, nodes j1 and j2 cannot be connected

dis=VER(j2,:)-VER(j1,:);
for k=2:na,
  position=VER(j1,:)+(alph(k)+alph(k-1))/2*dis;
  sa=solida(VER,ITRI,position);
  hoek=sum(sa);
  %[k,hoek]
  if abs(hoek) < .1, view=0; break, end
end

   

