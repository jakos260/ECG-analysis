function varargout=getAmplInfo(varargin) 


if length(varargin) < 1
	error('This routine needs the depolarization moments (1) and the adajcency matrix (2)');
else
	phih=varargin{1};
	if length(varargin) > 1
		givenJpoint=ceil(varargin{2});
	end
end
%  1	maximum P wave
%  2	vertex with maximum P wave


Pend=140;
PtoS=120;
info=zeros(1,8);
dJ=60;
%filter the egm with 100Hz
phih=lowpassma(phih,10);
absPhi=rms(phih);absPhi=baselinecor(absPhi);absPhi=abs(absPhi);
dPhi=diffrows(absPhi);
P=find(absPhi(1:Pend)==max(absPhi(1:Pend)));
R=find(absPhi==max(absPhi(1:Pend+PtoS)));
tPtoS=PtoS;
while R>=Pend+tPtoS || length(R)>1
	tPtoS=tPtoS-5;
	R=find(absPhi==max(absPhi(1:Pend+tPtoS)));
end
[mi,imi]=min(absPhi(R:R+150));
if imi>150 % find local minimum
	Jpoint=givenJpoint;
% 	[mi,imi]=min(baselinecor(absPhi(R:R+200)));
else
	Jpoint=imi+R;			
end

% future: add basleine correction
Send=Jpoint+dJ+30;
if Send>=length(absPhi)
	top=1
end

Tmommax=Send+find(dPhi(Send:end)==max(dPhi(Send:end)));
Tmommin=Send+find(dPhi(Send:end)==min(dPhi(Send:end)))-1;
if Tmommax >= length(phih)
	Tmom=Tmommin;
elseif Tmommin >= length(phih)
	Tmom=Tmommax;
else
	Tmom=max(Tmommax,Tmommin);
end
Tmom=min(Tmom,length(absPhi));

info(1)=absPhi(P);
info(2)=absPhi(R);
info(3)=absPhi(Tmom);
info(4)=absPhi(Send);
info(5)=P;
info(6)=R;
info(7)=Send;
info(8)=Tmom;

	
% clf;plot([1 :length(absPhi)],absPhi,'b:',[1 :length(absPhi)],phih,'b',imi+R-1,absPhi(imi+R-1),'ro',Send,absPhi(Send),'gs',Tmom,absPhi(Tmom),'kd','linewidth',2)	
% clf;plot([1 :length(absPhi)],absPhi,'b',P,absPhi(P),'yo',R,absPhi(R),'bo',Jpoint,absPhi(Jpoint),'ro',Send,absPhi(Send),'gs',Tmom,absPhi(Tmom),'kd','linewidth',2)	
	if Send > 400
		figure(400)
		clf;plot(1 :length(absPhi),phih(i,:),'b',1:length(absPhi),absPhi,'b:',P,absPhi(P),'yo',R,absPhi(R),'bo',Jpoint,absPhi(Jpoint),'ro',Send,absPhi(Send),'gs',Tmom,absPhi(Tmom),'kd','linewidth',2)	
	end

varargout{1}=info;
