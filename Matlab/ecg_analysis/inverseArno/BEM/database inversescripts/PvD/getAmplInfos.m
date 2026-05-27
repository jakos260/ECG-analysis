function varargout=getAmplInfos(varargin) 


if length(varargin) < 1
	error('This routine needs the depolarization moments (1) and the adajcency matrix (2)');
else
	phih=varargin{1};
	dogiven=0;
	if length(varargin) > 1
		givenJpoint=varargin{2};
		dogiven=1;
	end
end
%  1	maximum P wave
%  2	vertex with maximum P wave


Pend=140;
PtoS=160;
info=zeros(size(phih,1),4);
dJ=60;
%filter the egm with 200Hz
PHIH=lowpassma(phih,5);
phih=lowpassma(phih,15);

for i=1:size(phih,1)
	phi=(phih(i,:));
	if sum(phih(i,:))==0
		continue;
	end
	absPhi=baselinecor(phi);
	absPhi=abs(absPhi);
	dPhi=diffrows(absPhi);
	
	P=find(absPhi(1:Pend)==max(absPhi(1:Pend)));
	R=find(absPhi==max(absPhi(1:Pend+PtoS)));	
	tPtoS=PtoS;
	while R>=Pend+tPtoS +40 || length(R)>1
		tPtoS=tPtoS-5;
		R=find(absPhi==max(absPhi(1:Pend+tPtoS)));
	end
	if length(R) >1, R=R(1);elseif isempty(R), R=40;end
	% find J point in this signal
	if dogiven % find local minimum
		Jpoint=floor(givenJpoint);
	else
% 		[mi,imi]=min(absPhi(R:R+150));
		[mi,imi]=min(dPhi(R:R+150));
		if imi>150
			[mi,imi]=max(dPhi(R:R+150));
		end
		Jpoint=imi+R;			
	end
	Send=Jpoint+dJ;
	Tmom1=Send+find(dPhi(Send:end)==max(dPhi(Send:end)));	if Tmom1>=length(absPhi),Tmom1=Send;end
	Tmom2=Send+find(dPhi(Send:end)==min(dPhi(Send:end)));	if Tmom2>=length(absPhi),Tmom2=Send;end
	if absPhi(Tmom1) < absPhi(Tmom2)/10
		Tmom=Tmom2;
	elseif absPhi(Tmom2) < absPhi(Tmom1)/10
		Tmom=Tmom1;
	else
		Tmom=max(Tmom1,Tmom2);
	end
	
	info(i,1)=PHIH(i,P);
	info(i,2)=PHIH(i,R);
	info(i,3)=PHIH(i,Tmom);
	info(i,4)=PHIH(i,Send);
	info(i,5)=P;
	info(i,6)=R;
	info(i,8)=Send;
	info(i,7)=Tmom;
	if Send > 400
		figure(400)
		clf;plot(1 :length(absPhi),phih(i,:),'b',1:length(absPhi),absPhi,'b:',P,absPhi(P),'yo',R,absPhi(R),'bo',Jpoint,absPhi(Jpoint),'ro',Send,absPhi(Send),'gs',Tmom,absPhi(Tmom),'kd','linewidth',2)	
	end
	
% clf;plot([1 :length(absPhi)],absPhi,'b:',[1 :length(absPhi)],phih,'b',imi+R-1,absPhi(imi+R-1),'ro',Send,absPhi(Send),'gs',Tmom,absPhi(Tmom),'kd','linewidth',2)	
% clf;plot([1 :length(absPhi)],absPhi,'b:',P,absPhi(P),'yo',R,absPhi(R),'bo',Jpoint,absPhi(Jpoint),'ro',Send,absPhi(Send),'gs',Tmom,absPhi(Tmom),'kd','linewidth',2)	
end

varargout{1}=info;
