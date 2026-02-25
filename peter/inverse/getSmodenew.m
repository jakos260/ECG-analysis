function S=getSmode(T,dep,rep,pS,VER,ITRI,mode)
% variant of gets dedicated to ventricles
% function S=getSvnotch(T,dep,rep,pS,notchpot,scaleAmpl,mode)
% t: rowvector; dep, rep, notchpot, and scaleAmpl: columnvectors
% mode: 0 is dep only else AP complete
% PM van Dam; 2008-10-11
pup(1)= -4;	% determines the slope of the upstroke 
%plateau

nt=size(T,2);
[ni nj]=size(dep);if ni<nj, dep=dep';end
[ni nj]=size(rep);if ni<nj, rep=rep';end


%% depolarization
TDEP=T-dep*ones(1,nt);
S=1./(1+exp(pup(1)*TDEP));	

[area, tottriarea]= triareas(VER,ITRI);
area = area * 3;
triarea = ones(size(tottriarea,1),round(max(dep)) + 1);
verareas=zeros(size(area,1),size(S,2));

for i = 1:length(ITRI)
    i1 = ITRI(i,find(dep(ITRI(i,:)) == min(dep(ITRI(i,:)))));    
    i3 = ITRI(i,find(dep(ITRI(i,:)) == max(dep(ITRI(i,:)))));
    i2 = ITRI(i,:);
    i2(i2==i1 | i2==i3)=[];           a=[];
    for t=0:max(dep)+1
        if t < dep(i1)
            verareas([i1 i2 i3],t+1) = 0;
        else
            if t >= dep(i3)
                verareas(i1,t+1:end) = verareas(i1,t+1:end) + tottriarea(i);
                verareas(i2,t+1:end) = verareas(i2,t+1:end) + tottriarea(i);
                verareas(i3,t+1:end) = verareas(i3,t+1:end) + tottriarea(i);
                break;
            end
    
            if t <= dep(i2)
                verareas([i2 i3],t+1) = 0;
                rm = (VER(i2,:) - VER(i1,:)) * min(1.0,max(0,(t-dep(i1))) / (dep(i2)-dep(i1)));
                rp = (VER(i3,:) - VER(i1,:)) * min(1.0,max(0,(t-dep(i1))) / (dep(i3)-dep(i1)));
                triarea = norm(cross(rp,rm))/2.0;
                verareas(i1,t+1) = verareas(i1,t+1) + triarea;
                a=[a triarea];
            else
                verareas(i3,t+1) = 0;
                rm = (VER(i2,:) - VER(i3,:)) * min(1.0,max(0,(dep(i3)-t)) / (dep(i3)-dep(i2)));
                rp = (VER(i1,:) - VER(i3,:)) * min(1.0,max(0,(dep(i3)-t)) / (dep(i3)-dep(i1)));
                triarea = tottriarea(i) - norm(cross(rp,rm))/2.0;
                a=[a triarea];
                verareas(i1,t+1) = verareas(i1,t+1) + triarea;
                verareas(i2,t+1) = verareas(i2,t+1) + triarea;
            end      
        end
    end
end
S = S .* verareas ./ (area * ones(1,size(verareas,2)));

% S(:,1:size(verareas,2)) = S(:,1:size(verareas,2)) .* verareas;


if mode==1	
	S=diag(1./max(S,[],2))*S;
	return;
end
%% repolarization
% down
if mode>=4
    pdp(2)=pS(1);    % INPUT determines the slope leading up to the apex
    pdp(3)=pS(2);    % INPUT determines the (negative) slope following the apex
    TREP=T-rep*ones(1,nt); 

    Y= (1./(1+exp(pdp(2)*TREP))).*1./(1+exp(pdp(3)*TREP));
%     Y= (1./(1+exp(p3*TREP)))./(1+exp(p4*TREP)); 

    S=S.*Y;

    S=diag(1./max(S,[],2))*S;
end
return
%% plateau
if mode == 4
    if ~isempty(scaleAmpl)
        S=S.*(scaleAmpl * ones(1,size(S,2)));
    end
elseif mode==2
	plp=scaleAmpl;
	TPLATEAU=T-(dep+6)*ones(1,nt);	
	P=(plp*ones(1,nt))./(1+exp(-1*TPLATEAU)); 
	temp=max(S.*(1-P),[],2);
	for i=1:size(P,1)
		P(i,:)=1-(P(i,:)/temp(i));
	end
	S=S.*P;
	
	S=diag(1./max(S,[],2))*S;
	return
elseif mode==3
	N= (1./(1+exp(-0.05*(TDEP-25))))+(1.0./(1+exp(0.025*(TDEP-60))));
	N= N-1;
	N=1+scaleAmpl*diag(1./max(N,[],2))*N;
	S=S.*N;
	S=diag(1./max(S,[],2))*S;		
% 	return
elseif mode==5
	TPLATEAU=T-(dep+5)*ones(1,nt);	
	P=1-(scaleAmpl*ones(1,nt))./(1+exp(-1*TPLATEAU)); 
	S=S.*P;
	S=diag(1./max(S,[],2))*S;
	
	scaleAmpl=max(scaleAmpl*100/(rep-dep),0);
	N= (1./(1+exp(-0.15*(TDEP-25))))+(1.0./(1+exp(0.03*(TDEP-50))));
	N= N-1;
	N=1+scaleAmpl*diag(1./max(N,[],2))*N;
 	S=S.*N;
	S=diag(1./max(S,[],2))*S;	
end

	

	

% if size(T,1)==1, plot(Y,'g','linewidth',2);hold on; plot(P,'r','linewidth',2);end

% if exist('dodeprep','var') && dodeprep==1 && ~isempty(scaleAmpl)
% 	for i=1:size(S,1)
% 		S(i,:)=S(i,:)*scaleAmpl(i);
% 	end
% end

% S=S-0.85;