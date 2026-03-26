function varargout=nrClusters(GEOM,foci)
clusterDist=30;
varargout{1}=0;
if nargout>1
    varargout{2}=-1; 
end

if ~isempty(foci) & foci>0
	D=GEOM.DISTsurf(foci,foci);D(D>clusterDist)=0;D(D>0)=1;
	D=D+eye(size(D));
	CLUST=[];
	id=zeros(size(foci));
	idN=1;
	while sum(sum(D))>0
		sD=sum(D);
		i=find(sD==max(sD));i=i(1);
		ii=[i find(D(i,:)>0)];
		id(ii)=idN;idN=idN+1;
		clust.foci=foci(ii);
		clust.center=find(sum(GEOM.DISTsurf(foci(ii),:))==min(sum(GEOM.DISTsurf(foci(ii),:))));
		clust.center=clust.center(1);
		CLUST=[CLUST clust];
		rm=find(D(i,:)>0);
		D(i,:)=0;D(:,i)=0;
		D(rm,:)=0;	D(:,rm)=0;
	end
	varargout{1}=length(CLUST);
	if nargout>1
		varargout{2}=id; 
	end
end