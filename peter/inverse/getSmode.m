function S=getSmode(T,dep,rep,SPECS,mode,varargin)
% variant of gets dedicated to ventricles
% function S=getSvnotch(T,dep,rep,pS,notchpot,scaleAmpl,mode)
% t: rowvector; dep, rep, notchpot, and scaleAmpl: columnvectors
% mode: 0 is dep only else AP complete
% PM van Dam; 2008-10-11

% pS(3);    % INPUT determines the slope leading up to the apex
% pS(4);    % INPUT determines the (negative) slope following the apex


nt=size(T,2);
if size(dep,1) < size(dep,2) 
    dep=dep';
end
if size(rep,1) < size(rep,2) 
    rep=rep';
end


%% depolarization
% if size(SPECS.depSlope,1) == 1
%     slope = SPECS.depSlope * ones(size(dep));
% else
%     slope = SPECS.depSlope;
% end
% TDEP = T - dep*ones(1,nt);    
% GEOM=varargin{1};
% dept = (dep(GEOM.ITRI(:,1),:) + dep(GEOM.ITRI(:,2),:) + dep(GEOM.ITRI(:,3),:))/3;
% for i=1:length(dep)
%     [a,~]=find(GEOM.ITRI==i);
%     b = dept(a) - dep(i);
%     dT = diff(range(b)); 
%     slope(i) = min(slope(i), 6.0/dT);
% end
% S = 1./ (1 + exp(bsxfun(@times,TDEP,-slope(:))));

% if length(varargin) == 1
%     
%     neigh=varargin{1};
%     S = getSpc(T,dep,rep,SPECS,mode,neigh);
% %     S = getSmodeHakkel(T,dep,rep,SPECS,mode,varargin{1});
%     return;
% end

TDEP = T - dep*ones(1,nt);
if size(SPECS.depSlope,1) == 1
    S=1./(1+exp(-SPECS.depSlope.*TDEP));
else
    S = 1./ (1 + exp(bsxfun(@times,TDEP,-SPECS.depSlope(1))));
end




% TDEP = T - dep*ones(1,nt);
% if size(SPECS.depSlope,1) == 1
%     S=1./(1+exp(-SPECS.depSlope.*TDEP));
% else
%     S = 1./ (1 + exp(bsxfun(@times,TDEP,-SPECS.depSlope(:))));
% end



%% repolarization
% down
if mode>=4
    TREP=T - rep*ones(1,nt); 
       

    if SPECS.useCumsum
        Y = (SPECS.initialSlope + 1./(1 + exp( SPECS.plateauslope.* TREP )))./(1 + exp(SPECS.repslope .*TREP)); 
        Y = 1 - bsxfun(@times,cumsum(Y,2),1./sum(Y,2));
    else
%         Y = -(SPECS.initialSlope + 1./(1 + exp( SPECS.plateauslope.* TREP )))./(1 + exp(SPECS.repslope .*TREP)); 
        Y = (1./(1 + exp( SPECS.plateauslope.* TREP )))./ (1 + exp(SPECS.repslope .*TREP)); %ECGSIM
        Y = bsxfun(@times,Y,1./Y(:,1));
    end
    
    % apply the up-slope: logistic shape;   
    S = S .* Y;
    if isfield(SPECS,'scaleAmpl') && ~isempty(SPECS.scaleAmpl)
        S=bsxfun(@rdivide,S,max(S,[],2)./SPECS.scaleAmpl); % oostep1: 20121015 much faster
    else
        S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return
% 
% if 0
%     VER= GEOM.VER;
%     ITRI= GEOM.ITRI;
%     [area, tottriarea]= triareas(VER,ITRI);
%     verareas=zeros(size(area,1),size(S,2));
%     for i = 1:size(ITRI,1)
%         i1 = ITRI(i,find(dep(ITRI(i,:)) == min(dep(ITRI(i,:)))));
%         i3 = ITRI(i,find(dep(ITRI(i,:)) == max(dep(ITRI(i,:)))));
%         if length(i1) > 1
%             i2 = i1(2);
%             i1 = i1(1);
%         elseif length(i3) > 1
%             i2 = i3(1);
%             i3 = i3(2);
%         else
%             i2 = ITRI(i,:);
%             i2(i2==i1 | i2==i3)=[];
%         end
%         for t=0:max(dep)+1
%             if t < dep(i1)
%                 verareas([i1 i2 i3],t+1) = 0;
%             else
%                 if t >= dep(i3)
%                     verareas(i1,t+1:end) = verareas(i1,t+1:end) + tottriarea(i);
%                     verareas(i2,t+1:end) = verareas(i2,t+1:end) + tottriarea(i);
%                     verareas(i3,t+1:end) = verareas(i3,t+1:end) + tottriarea(i);
%                     break;
%                 end
%                 if t <= dep(i2)
%                     verareas([i2 i3],t+1) = 0;
%                     rm = (VER(i2,:) - VER(i1,:)) * min(1.0,max(0,(t-dep(i1))) / (dep(i2)-dep(i1)));
%                     rp = (VER(i3,:) - VER(i1,:)) * min(1.0,max(0,(t-dep(i1))) / (dep(i3)-dep(i1)));
%                     triarea = norm(cross(rp,rm))/2.0;
%                     triarea = triarea;
%                     verareas(i1,t+1) = verareas(i1,t+1) + triarea;
%                 else
%                     verareas(i3,t+1) = 0;
%                     rm = (VER(i2,:) - VER(i3,:)) * min(1.0,max(0,(dep(i3)-t)) / (dep(i3)-dep(i2)));
%                     rp = (VER(i1,:) - VER(i3,:)) * min(1.0,max(0,(dep(i3)-t)) / (dep(i3)-dep(i1)));
%                     triarea = tottriarea(i) - norm(cross(rp,rm))/2.0;
%                     triarea = triarea;
%                     verareas(i1,t+1) = verareas(i1,t+1) + triarea;
%                     verareas(i2,t+1) = verareas(i2,t+1) + triarea;
%                 end
%             end
%         end
%     end   
% %     S = S .* verareas ./ (3.0*area * ones(1,size(verareas,2)));
%     verareas=bsxfun(@rdivide,verareas,3.0*area); % oostep1: 20121015 much faster
%     S = S .* verareas;
% end
% if mode==1
%     % 	S=diag(1./max(S,[],2))*S;
%     S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster
% end
% 
% %% plateau
% if mode == 4
%     if ~isempty(scaleAmpl)
%         S=S.*(scaleAmpl * ones(1,size(S,2)));
%     end
% elseif mode==2
%     plp=scaleAmpl;
%     TPLATEAU=T-(dep+6)*ones(1,nt);
%     P=(plp*ones(1,nt))./(1+exp(-1*TPLATEAU));
%     temp=max(S.*(1-P),[],2);
%     for i=1:size(P,1)
%         P(i,:)=1-(P(i,:)/temp(i));
%     end
%     S=S.*P;
%     
%     % 	S=diag(1./max(S,[],2))*S;
%     S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster
%     
%     
%     return
% elseif mode==3
%     N= (1./(1+exp(-0.05*(TDEP-25))))+(1.0./(1+exp(0.025*(TDEP-60))));
%     N= N-1;
%     N=1+scaleAmpl*diag(1./max(N,[],2))*N;
%     S=S.*N;
%     % 	S=diag(1./max(S,[],2))*S;
%     S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster
%     
%     % 	return
% elseif mode==5
%     TPLATEAU=T-(dep+5)*ones(1,nt);
%     P=1-(scaleAmpl*ones(1,nt))./(1+exp(-1*TPLATEAU));
%     S=S.*P;
%     % 	S=diag(1./max(S,[],2))*S;
%     S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster
%     
%     
%     scaleAmpl=max(scaleAmpl*100/(rep-dep),0);
%     N= (1./(1+exp(-0.15*(TDEP-25))))+(1.0./(1+exp(0.03*(TDEP-50))));
%     N= N-1;
%     N=1+scaleAmpl*diag(1./max(N,[],2))*N;
%     S=S.*N;
%     % 	S=diag(1./max(S,[],2))*S;
%     S=bsxfun(@rdivide,S,max(S,[],2)); % oostep1: 20121015 much faster
%     
% end
% 
% 
% 
% 
% 
% % if size(T,1)==1, plot(Y,'g','linewidth',2);hold on; plot(P,'r','linewidth',2);end
% 
% % if exist('dodeprep','var') && dodeprep==1 && ~isempty(scaleAmpl)
% % 	for i=1:size(S,1)
% % 		S(i,:)=S(i,:)*scaleAmpl(i);
% % 	end
% % end
% 
% % S=S-0.85;