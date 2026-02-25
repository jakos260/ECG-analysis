function V = lsh_read(fname);

% open file
fid = fopen(fname,'r');

%%%% read header

s = fgetl(fid);
[type,s] = strtok(s(2:end));
[fmt,s] = strtok(s);

fgetl(fid);

s = fgetl(fid);
[geom,s] = strtok(s(2:end));
if strcmp(geom,'tissue2d')
	t_size = sscanf(s,'%f %f');
elseif strcmp(geom,'tissue1d')
	t_size = sscanf(s,'%f');
else
	%disp('lsh_read: unknown geometry');
end

fgetl(fid);

s = fgetl(fid);
n = sscanf(s(2:end),'%d');


if strcmp(type,'vmmap')

	if fmt == 'txt'
		V = load(fname);
	elseif fmt == 'bin'
		for i=1:n, fgetl(fid); end
		V = fread(fid,inf,'float');
	elseif fmt == 'mtx'
		s = fgetl(fid);
		Ne = sscanf(s(2:end),'%d');
		for i=1:n-1
			fgetl(fid);
		end
		V = fread(fid,inf,'float');
		Vn = (Ne+1) * floor(length(V)/(Ne+1));
		V = reshape(V(1:Vn),Ne+1,floor(length(V)/(Ne+1)))';
	end
	
	if strcmp(geom,'tissue2d')
		V = reshape(V,t_size(2),t_size(1));
	elseif strcmp(geom,'tissue1d')
	else
		%disp('lsh_read: unknown geometry');
	end

elseif strcmp(type,'cmap')

	fgetl(fid);
	s = fgetl(fid);
	x = sscanf(s(2:end),'%f %f %f');
	
	V = fread(fid,inf,'char');
	V = x(2) + V/255 * (x(3)-x(2));
	
	if strcmp(geom,'tissue2d')
		V = reshape(V,t_size(2),t_size(1));
	elseif strcmp(geom,'tissue1d')
	else
		%disp('lsh_read: unknown geometry');
	end

elseif strcmp(type,'spoon')
	
	fgetl(fid);
	s = fgetl(fid);
	Ne = sscanf(s(2:end),'%d');
	for i=1:n-2
		fgetl(fid);
	end
	
	V = fread(fid,inf,'float');
	Vn = (Ne+1) * floor(length(V)/(Ne+1));
	V = reshape(V(1:Vn),Ne+1,floor(length(V)/(Ne+1)))';
	
else
	%disp('lsh_read: unknown file format');
	V = load(fname);
end

fclose(fid);

