

function varargout = wallthickness(varargin)

% types identifies atria (1-3), ventricles(4-6) or other.
% only for the atria and ventricles the is determined
% for teh atria the wall thickness is set to 2 mm except for wall
% thicknesses > 8 mm (septum) Type is only used for one output parameter,
% i.e the wall thickness

maxWallThick = 30;
VER = varargin{1};
ITRI = varargin{2};
if nargin == 3
    maxWallThick = varargin{3};
end


if nargout == 3
    wd=zeros(length(VER),1);
    wopVer=zeros(length(VER),3);
    woptri=zeros(length(VER),1);
    for i=1:length(wd)
        normal = VertexNormal(VER,ITRI,i);
        v1 = VER(i,:);
        v2 = v1 - normal;
        TR = linetris(VER,ITRI,v1,v2);
        TR(TR(:,5)<0.01,:)=[];
        TR(TR(:,5) > min(TR(:,5)),:)=[];
        if ~isempty(TR)
            wd(i) = TR(1,5);
            wopVer(i,:) = v1 - normal * TR(5);
            woptri(i,:) = TR(1);
        end
    end
    varargout{1}= wd;
    varargout{2}= wopVer;
    varargout{3}= woptri;   
elseif nargout==2
    wd=zeros(length(VER),1);
    wopVer=zeros(length(VER),3);
    for i=1:length(wd)
        normal = VertexNormal(VER,ITRI,i);
        v1 = VER(i,:);
        v2 = v1 - normal;
        TR = linetris(VER,ITRI,v1,v2);
        TR(TR(:,5)<0.01,:)=[];
        TR(TR(:,5) > min(TR(:,5)),:)=[];
        if ~isempty(TR)
            wd(i)=TR(1,5);
            wopVer(i,:) = v1 - normal * TR(5);
        end
    end
    varargout{1}= wd;
    varargout{2}= wopVer;
elseif nargout==1
    wd=zeros(length(VER),1);
    for i=1:length(wd)
        normal = VertexNormal(VER,ITRI,i);
        v1 = VER(i,:);
        v2 = v1 + normal;
        TR = linetris(VER,ITRI,v1,v2);
        TR(abs(TR(:,5))<1e-6,:)=[];
        TR((TR(:,5))>1e-6,:)=[];
        TR(abs(TR(:,5))>maxWallThick,:)=[];
        if ~isempty(TR) 
            % the oposite wall should be a negative distance value
            % (oposite normal direction) 
            if all(TR(:,5)> 0 )
                a =TR(TR(:,5) == min(TR(TR(:,5)> 0 ,5)),5);                    
                wd(i)= -a(1);
            else
                a = TR(TR(:,5) == max(TR(TR(:,5)<0,5)),5);
                wd(i) = -a(1);
            end
        else
            wd(i) = maxWallThick;
        end
    end
    varargout{1}= wd;
else
    error('incorrect number of output arguments');
end    