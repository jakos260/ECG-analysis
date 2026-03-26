function [M,R] = loadJHUECG(filename)


fd = fopen(filename);
if fd == -1
    error('could not open file');
end

% 19 Sample Interval: 2 ms 

% 21 Median signal: 600 X 12 

for i=1:19
    str = fgetl(fd);
end

sampleT = str2double(str(17:end-3));
str = fgetl(fd);
str = fgetl(fd);
nSamp = str2double(str(15:end-5));
str = fgetl(fd);
M=[];
for i=1:nSamp
    str= fgetl(fd);
    M=[M; str2num(str)];
end
str = fgetl(fd);
str = fgetl(fd);
str = fgetl(fd);
nSamp = str2double(str(15:end-5));
str = fgetl(fd);
R=[];
for i=1:nSamp
    str= fgetl(fd);
    R=[R; str2num(str)];
end
M=M';
R=R';



