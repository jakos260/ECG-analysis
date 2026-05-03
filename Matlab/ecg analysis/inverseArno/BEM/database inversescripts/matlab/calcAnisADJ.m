function ADJ2=calcAnisADJ(varargin)

%   ('ADJ2=calcAnisADJ(GEOM,anisotropyRatio) or') 
%   ('ADJ2=calcAnisADJ(VER,ITRI,ADJ,RVER,LVER,anisotropyRatio)')    
%   ('ADJ2=calcAnisADJ(VER,ITRI,RVER,LVER,anisotropyRatio)') 
%   ('ADJ2=calcAnisADJ(VER,ITRI,ADJ,vertype,anisotropyRatio)')            
%   ('ADJ2=calcAnisADJ(VER,ITRI,vertype,anisotropyRatio)')  
%   ('vertype==0 epiver; vertype==1 is right endover; vertype==2 left endover');



if length(varargin)<2 
    disp('ADJ2=calcAnisADJ(GEOM,anisotropyRatio) or') 
    disp('ADJ2=calcAnisADJ(VER,ITRI,ADJ,RVER,LVER,anisotropyRatio)')    
    disp('ADJ2=calcAnisADJ(VER,ITRI,RVER,LVER,anisotropyRatio)') 
    disp('ADJ2=calcAnisADJ(VER,ITRI,ADJ,vertype,anisotropyRatio)')            
    disp('ADJ2=calcAnisADJ(VER,ITRI,vertype,anisotropyRatio)')  
    disp('vertype==0 epiver; vertype==1 is right endover; vertype==2 left endover');
	error('Not enough parameters!!!  ');
elseif length(varargin)>=3 
    type=[];
	VER=varargin{1};
	ITRI=varargin{2};
    tmp=varargin{3};
    pp=3;
    if size(tmp,1)==size(tmp,2) && size(tmp,1)==length(VER) % next parameter is ADJ
        ADJ=varargin{pp};
        pp=pp+1;
    end
    tmp=varargin{pp};
    if min(size(tmp))==3  % next parameters are the left and right vertices
        RVER=varargin{pp};
        pp=pp+1;
        if length(varargin)<pp+1
            error('not enough parameters')
        end
        if min(size(varargin{pp}))==3
            LVER=varargin{pp};pp=pp+1;
        else
             error('unknown parameter')
        end
        anisotropyRatio=varargin{pp};
    elseif min(size(tmp))==1 %read vertype
        type=varargin{pp};
        pp=pp+1;
        anisotropyRatio=varargin{pp};
        RendoVER=zeros(size(type));
        LendoVER=zeros(size(type));
        RendoVER(type==1)=1;
        LendoVER(type==2)=1;
    else
        error('unknown parameter');
    end
else %ADJ2=calcAnisADJ(GEOM,anisotropyRatio)
   type=[];
   GEOM=varargin{1};
   VER=GEOM.VER;
   RVER=GEOM.RVER;   
   LVER=GEOM.LVER;
   ITRI=GEOM.ITRI;
   ADJ=GEOM.ADJ;
   anisotropyRatio=varargin{2};
end


if anisotropyRatio==1
	ADJ2=ADJ;
	return;
end

%% determine the endocardial and epicardial parts
if isempty(type)
    RendoVER=zeros(1,length(VER));
    for i=1:length(VER)
        if ~isempty(find(RVER(:,1)==VER(i,1) & ...
                        RVER(:,2)==VER(i,2) & ...
                        RVER(:,3)==VER(i,3),1))
            RendoVER(i)=1;
        end
    end
    LendoVER=zeros(1,length(VER));
    for i=1:length(VER)
        if ~isempty(find(LVER(:,1)==VER(i,1) & ...
                         LVER(:,2)==VER(i,2) & ...
                        LVER(:,3)==VER(i,3),1))
            LendoVER(i)=1;
        end
    end
    endoVER=LendoVER+RendoVER;
% epiVER=ones(size(endoVER));epiVER(endoVER==1)=0;
    type=endoVER+LendoVER;
end 

% type(0)= epiver; type(1) is right endover; type(2)= left endover

%%
RendoITRI=zeros(length(ITRI),1);
for i=1:length(ITRI)
    if sum(RendoVER(ITRI(i,:)))==3
        RendoITRI(i)=1;
    end
end
LendoITRI=zeros(length(ITRI),1);
for i=1:length(ITRI)
    if sum(LendoVER(ITRI(i,:)))==3
        LendoITRI(i)=1;
    end
end
epiITRI=ones(length(ITRI),1);
epiITRI(LendoITRI==1)=0;
epiITRI(RendoITRI==1)=0;

lver=length(VER);

%% calculate wall thicknesses

% [wd,wopVer]=wallthickness(VER,ITRI);
if sum(RendoITRI) > 0
    [dR,pntR]=pnttridist(VER,ITRI(RendoITRI==1,:),VER);
end
if sum(LendoITRI) > 0
    [dL,pntL]=pnttridist(VER,ITRI(LendoITRI==1,:),VER);
end

[dE,pntE]=pnttridist(VER,ITRI(epiITRI==1,:),VER);


%% 
adj0=graphdist(ITRI,VER,4);
if exist('ADJ')
	adj0w=ADJ;
else
	adj0w=graphdist(ITRI,VER,3);
end
%% remove connections between points on the same type of card, i.e.
%  endocard to endocard, epicard to epicard

for i=1:lver
	adj0w(i,type==type(i))=0;
	adj0w(type==type(i),i)=0;
end

adj0w(adj0>0)=adj0(adj0>0);
ADJ2=adj0w;

%% Adapt the connections of the adajcency matrix through the wall

anis=calcVlongwall(anisotropyRatio);
for i=1:length(ADJ2)
	for j=i:length(ADJ2)
		if ADJ2(i,j)>0 && type(i)~=type(j)
            switch type(j) %type(j) indicates the ..card i is connected to
				case 0
					mdi = dE(i);
					pnti=pntE(i,:);
				case 1
					mdi = dR(i);
					pnti=pntR(i,:);
				case 2
					mdi = dL(i);
					pnti=pntL(i,:);
				otherwise
					disp('error')
            end
            switch type(i) %type(i) indicates the ..card j is connected to
                case 0
					mdj = dE(j);
					pntj=pntE(j,:);
				case 1
					mdj = dR(j);
					pntj=pntR(j,:);
				case 2
					mdj = dL(j);
					pntj=pntL(j,:);
				otherwise
					disp('error')
			end
			
			if 0
				% TODO This part assumes an fiber rotation of 120 degrees and
				% takes into account the reduction inporpgation velocity
				% caused by this fiber rotation
				d1=anis*norm3d(VER(j,:)-pnti);
				d2=anis*norm3d(VER(i,:)-pntj);
			else
				% an error is caused because the distance along the 
				% wall surface isapproximated by a straight line. 
				% This accounted for by a factor 1.1, assuming an average
				% error of 10%
				d1 = 1.1*norm3d(VER(j,:)-pnti);
				d2 = 1.1*norm3d(VER(i,:)-pntj);                
            end
            h = mean([mdi mdj]);
            md = sqrt(adj0w(i,j)^2 -  h^2); 
            
            if exist('GEOM','var') && isfield(GEOM,'Rfreewallver') && GEOM.Rfreewallver(i) && GEOM.Rfreewallver(j)
%                 d=max(adj0w(i,j),sqrt(min([d1,d2])^2+(anisotropyRatio * mean([min(mdi,4.0),min(mdj,4.0)]))^2));
                d=sqrt(min([d1,d2])^2+(anisotropyRatio * mean([min(mdi,4.0),min(mdj,4.0)]))^2);
            else
%                 d=max(adj0w(i,j),sqrt(min([d1,d2])^2+(anisotropyRatio*mean([mdi,mdj]))^2));

%                 d=sqrt(min([d1,d2])^2+(anisotropyRatio *h)^2);
                d=sqrt(md^2+(anisotropyRatio *h)^2);
            end
            ADJ2(i,j) = d; 
            ADJ2(j,i) = d;                
		end
	end
end

%%
function anis=calcVlongwall(anisotropyRatio)
alpha=0:pi/3000:pi*120/180;  % assume rotation of 120 degrees
anis=1/mean(sqrt(cos(alpha).^2+(sin(alpha)./anisotropyRatio).^2));
