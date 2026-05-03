% gets.m
% variant of gets
% function S=gets(t,dep,rep,p,mode)
% t: rowvector; dep and rep: columnvectors
% A. van Oosterom; 2008-12-04

function S=getspvdnew(T,dep,rep,p,VER,ITRI,mode)
nt=size(T,2);
nn=size(T,1);

if size(dep,1)<size(dep,2), dep=dep'; end
if size(rep,1)<size(rep,2), rep=rep'; end


TDEP=T-dep*ones(1,nt);


% p(1): slope upstroke)*4;
% p(2): init value Y_dom, setting initial downslope of TMP
% p(3): parameter setting leading curvature of Tdom;
% p(4): parameter setting trailing curvature of Tdom;
% p(5): extra shift easing the timing of the apex to coincide toward that of
%       apex Tdom

if mode~=1,
    TREP=T-rep*ones(1,nt)-p(5);
    p3=ones(nn,1)*p(3);
    p4=ones(nn,1)*p(4);
    
    % compute (-1*) derivative of the downward slope of the TMP;
    Y=(p(2)+1./(1+exp(diag(p3)*TREP)))./(1+exp(diag(p4)*TREP));
    
    % compute TMP; unit upstroke
    Y=(1-diag(1./sum(Y'))*cumsum(Y')');
    % apply the up-slope: logistic shape;
    
    
    
    %     S=Y./(1+exp(-p(1)*TDEP));
    S=Y;
    TTRIDEP= ones(size(ITRI,1),1) * T(1,:) - min( dep(ITRI),[],2) * ones(1,size(T,2));
    Stridep = 1./(1+exp(-p(1)*TTRIDEP));
    if 1
        [area, tottriarea]= triareas(VER,ITRI);
        verareas=zeros(size(area,1),size(S,2));
        for i = 1:size(ITRI,1)
            i1 = ITRI(i,find(dep(ITRI(i,:)) == min(dep(ITRI(i,:)))));
            i3 = ITRI(i,find(dep(ITRI(i,:)) == max(dep(ITRI(i,:)))));
            if length(i1) > 1
                i2 = i1(2);
                i1 = i1(1);
            elseif length(i3) > 1
                i2 = i3(1);
                i3 = i3(2);
            else
                i2 = ITRI(i,:);
                i2(i2==i1 | i2==i3)=[];
            end
            for t=0:max(dep)+1
                if t < dep(i1)
                    verareas([i1 i2 i3],t+1) = 0;
                else
                    if t >= dep(i3)
                        verareas(i1,t+1:end) = verareas(i1,t+1:end) + tottriarea(i)* Stridep(i,t+1:end);
                        verareas(i2,t+1:end) = verareas(i2,t+1:end) + tottriarea(i)* Stridep(i,t+1:end);
                        verareas(i3,t+1:end) = verareas(i3,t+1:end) + tottriarea(i)* Stridep(i,t+1:end);
                        break;
                    end
                    if t <= dep(i2)
                        verareas([i2 i3],t+1) = 0;
                        rm = (VER(i2,:) - VER(i1,:)) * min(1.0,max(0,(t-dep(i1))) / (dep(i2)-dep(i1)));
                        rp = (VER(i3,:) - VER(i1,:)) * min(1.0,max(0,(t-dep(i1))) / (dep(i3)-dep(i1)));
                        triarea = norm(cross(rp,rm))/2.0;
                        triarea = triarea * Stridep(i,t);
                        verareas(i1,t+1) = verareas(i1,t+1) + triarea* Stridep(i,t);
                    else
                        verareas(i3,t+1) = 0;
                        rm = (VER(i2,:) - VER(i3,:)) * min(1.0,max(0,(dep(i3)-t)) / (dep(i3)-dep(i2)));
                        rp = (VER(i1,:) - VER(i3,:)) * min(1.0,max(0,(dep(i3)-t)) / (dep(i3)-dep(i1)));
                        triarea = tottriarea(i) - norm(cross(rp,rm))/2.0;
                        triarea = triarea * Stridep(i,t);
                        verareas(i1,t+1) = verareas(i1,t+1) + triarea* Stridep(i,t);
                        verareas(i2,t+1) = verareas(i2,t+1) + triarea* Stridep(i,t);
                    end
                end
            end
        end   
        S = S .* verareas ./ (3.0*area * ones(1,size(verareas,2)));
    end
    % re-establish unit upstroke
    S=diag(1./max(S,[],2)) * S;
else
    S=1./(1+exp(-p(1)*TDEP));
end

