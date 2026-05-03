function [VER, ITRI, funout] = makecatheter(POS,id,funin,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Name:    makecatheter
% 
% Desc:    Morph a cylinder to match catehter electrode positions
% 
% Inputs:  POS: position of electrodes, funin: values (e.g. activation time) at electrode
% positions, set unknown values to NaN. If funin is empty of all NaN, all cathere positions will be drawn, else only from first ~NaN to last ~NaN
% , remainng NaN values in funin will be interpolated. r: radius of catheter. 
% 
% Outputs: VER, ITRI geometry of morphed 'catheter', funout: function
% values at vertices in VER. Empty if no funin values ~nan were given.
% 
% Created: 06-Nov-2013 12:33:40
% 
% Author:  peteroosterhoff
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f=3; % interpolation factor

if ~exist('funin','var') || isempty(funin)
    funin=nan(1,size(POS,1));
end
firstl=1;
lastl=size(POS,1);

if ~isempty(funin) && any(~isnan(funin))
    % interpolate, but do not extrapolate
    firstel=find(~isnan(funin),1,'first');
    lastel=find(~isnan(funin),1,'last');
    POS=POS(firstel:lastel,:);
    funin=funin(firstel:lastel);
    % create intermediate points, TBD: only if dist>xmm/every xmm
%     funin=[funin;nan(1,length(funin))];
%     funin=funin(1:end-1);
%     intPOS=[(POS(2:end,:)+POS(1:end-1,:))/2; 0 0 0];
%     POS=reshape([POS(:) intPOS(:)]',2*size(POS,1),[]);
%     POS=POS(1:end-1,:);
    
    dpos=cumsum([ 0 rownorm(POS(2:end,:)-POS(1:end-1,:))]); % distance along catheter
    ii=1:1/f:length(dpos);
    dposi=interp1(1:length(dpos),dpos,ii,'linear');
%     dposi=lowpassma(reshape([dpos;dpos2],1,[]),2); % add intermediate points
%     dposi(1)=[];
%     POSi=interp1(dpos,POS,dposi,'pchip');
    
    
    funin=interp1(dpos(~isnan(funin)),funin(~isnan(funin)),dposi,'pchip'); % linear interpolation along catheter
end

nphi=10;


[VER,ITRI]=make_snake(POS,
% [x, y, z]= gencyl( POS' , r*ones(1,size(POS,1)) , f , nphi );
% VER=[reshape(x',[],1),reshape(y',[],1),reshape(z',[],1)];
% 
% VER(nphi+1:nphi+1:end,:)=[]; % remove double vertices 


if isempty(funin)|| all(isnan(funin)) 
    funout=[]; %
else
    
    for i=1:(size(VER,1)-2)/nphi
        funout(2+(i-1)*nphi:1+i*nphi)=funin(i);
    end
    funout(1)=funin(1);
    funout(end+1)=funin(end);
end

    
    



