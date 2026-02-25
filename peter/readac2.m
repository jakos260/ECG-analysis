function M = readac2(filename)

f=fopen(filename,'r','b');
if (f==-1)
  fprintf('\nCannot open %s\n\n', name);
  M=0;
  extraresult='';
  return;
end

numMuxChannels = fread(f,1,'uint16');
numheadbytes   = fread(f,1,'uint16');
A = fread(f,602,'uint8');
nrleads = fread(f,1,'uint16');
nrFrames = fread(f,1,'uint');
A = fread(f,22,'uint8');
sampRate = fread(f,1,'int16');
gain = fread(f,1,'int16');
A = fread(f,numheadbytes-638,'uint8');


M=fread(f,'int16');
M = reshape(M,numMuxChannels,length(M)/numMuxChannels);
fclose(f);