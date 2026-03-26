function nrrdwrite(X,meta,filename)

% this function writes a nrrd file  with the NRRD0001 format.
% usage: nrrdwrite(X,meta,filename) where: 
%   X - data vector to save out.  should be the shaped geometrically, ie.
%   3D or 4D.  
%   meta - header information as a structure with strings.
%   filename - file name for writing the nrrd



obs.size=size(X);
obs.dim=length(obs.size);


if ~isfield(meta, 'nrrd_type') || length(meta.nrrd_type)~=8 ||...
        ~strcmp(meta.nrrd_type(1:4),'NRRD')
    meta.nrrd_type='NRRD0001';
end

if ~isfield(meta, 'encoding')
    meta.encoding='raw';
end

if ~isfield(meta, 'endian')
    [~,~,endian] = computer();
    if endian=='L'
        meta.endian='little';
    elseif endian=='B'
        meta.endian='big';
    else
        error('cannot get machine format')
    end
end

if ~isfield(meta, 'type')
    tmps=whos('X');
    meta.type=tmps.class;
end


if ~isfield(meta, 'dimension')
    meta.dimension = sprintf('%u',obs.dim);
elseif isfield(meta, 'dimension') && min(obs.dim==str2double(meta.dimension))==0
    warning('Fixing meta dimension.')
    meta.dimension = sprintf('%u',obs.dim);
end


if ~isfield(meta, 'sizes')
    meta.sizes= sprintf('%u ',obs.size);
else
    if ischar(meta.sizes)
        declared_size=str2num(meta.sizes);
    else
        declared_size=meta.sizes;
    end
    
    if min(declared_size==obs.size)==0 && min(obs.size)==1 && obs.dim==2
        X=reshape(X,declared_size);
        obs.size=size(X);
        obs.dim=length(obs.size);
    elseif min(declared_size==obs.size)==0 
        declared_size=obs.size;
        meta.sizes=sprintf('%u ',obs.size);
        warning('Meta size did not match the size of matrix.  Using size of matrix')
    end
end
    

header=sprintf([ meta.nrrd_type '\n'...
    '# Complete NRRD file format specification at:\n'...
    '# http://teem.sourceforge.net/nrrd/format.html;\n']);


meta=rmfield(meta,'nrrd_type');

fields=fieldnames(meta);

for k=1:length(fields)
    field=fields{k};
    field_content=meta.(fields{k});
    
    field(field=='_')=' ';
    
    if ~ischar(field_content)
        if min(size(field_content))>1 && length(size(field_content))>2
            error('any meta data field must be 1D');
        end
        
        if size(field_content,1)~=1
            field_content=field_content';
        end
        
        field_content=num2str(field_content);
        
    end
    
    header=[header sprintf('%s: %s\n',field,field_content)];
end

header=sprintf('%s\n',header);


    
if strcmpi(meta.endian,'big')
    endian_str='b';
elseif strcmpi(meta.endian,'little')
    endian_str='l';
else
    error('cannot recognize machine format.  try "little" or "big"')
end

    


fid=fopen(filename,'wt',endian_str);
fprintf(fid,'%s',header);
fseek(fid,0,'eof');
fwrite(fid,X,meta.type);
fclose(fid);


end
        
        
    
    


