function [elcs,ver,itri]=readMarkelcs(fn,VER)
elcs=[];
fid=fopen(fn);
if (fid==-1)
%   error(['Cannot open: ' fn]);
	return
end
ver=[];
itri=[];
id=[];
for i=1:6,fgetl(fid);end
while ~feof(fid)
	id=[id; fscanf(fid,'%d',1)];
	ver=[ver; fscanf(fid,'%f',3)'];
	itri=[itri; fscanf(fid,'%f',1)'];	
	fgetl(fid);
end
fclose(fid)
elcs=zeros(size(id));
for i=1:length(ver)	
	d=norm3d([VER(:,1)-ver(i,1) VER(:,2)-ver(i,2) VER(:,3)-ver(i,3)]);
	disp(num2str(min(d)));
	elcs(i)=find(d==min(d));
end

