function status=CheckInverseFiles(subject,type,geomdir,amapath)



% date of ventricle or atrium
type=regexprep(type,'^_',''); % remove leading _
D=dir(fullfile(geomdir,subject,[subject '_' type '.tri']));
sourcegeomdate=D.datenum;

% dates of distance and adjacency matrices
extlist={'.adj2d','.dst2d','.adj3d','.dst3d','.adjanis','.dstanis','.adjneigh','.ordneigh'};
for i=1:length(extlist)
    dstadjpath=fullfile(geomdir,subject,[subject '_' type,extlist{i}]);
    dstadjinfo=dir(dstadjpath);
    if isempty(dstadjinfo) || dstadjinfo.datenum-sourcegeomdate<-60/(24*3600) % max 60 time diff
        error('Distance/adjacency matrix %s outdated or missing. Run export and create in GeomPeacs',dstadjinfo.name);
    end
end

% date of amatrix
D=dir(amapath);
if isempty(D)
    error('A-matrix is missing');
end
amadate=D(1).datenum;
fid=fopen(amapath);
amatxt=fread(fid,'uint8=>char')';
fclose(fid);

% compare amadate to all used tri files. Type files ignored
D=dir(fullfile(geomdir,subject,[subject '_*.tri']));
if ~isstruct(D)
    D={D};
end
for i=1:length(D)
    triinfo=dir(fullfile(geomdir,subject,D(i).name));
    if triinfo.datenum>amadate && ~isempty(strfind(amatxt,triinfo.name))
        error('A matrix older than %s. Run doama/double.', triinfo.name);
    end
end
status=1;



