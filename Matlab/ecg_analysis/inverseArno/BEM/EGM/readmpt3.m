% function Point = readmpt3(filename)
%
% Trace signals in order of appearance: I,II,III,aVR,aVL,aVF,
% V1,V2,V3,V4,V5,V6,M1,M2,M3,M4,R1,R2,M1-M2,M3-M4,R1-R2,?,?,?,?,?,?
%
% M1 is tip of catheter, bipolar leads are filtered differently from
% unipolar, hence e.g. signal M1-M2 is not the same as difference of
% M1 and M2.
% ldata is index of rdata that corresponds to sample in Trace.
% basically [1:10:2500];
% rdata 1:3 probably 3D pos of reference
% rdata 4:6 probably 3D orientation of reference (in r,theta,phi ?)
% rdata 7:9 probably 3D pos of catheter tip
% rdata 10:12 probably 3D orientation of tip (in r,theta,phi ?)
% f1 is probably UnitPerBit.
%
% l1 unknown (pointer?)
%
%
% Version 1: 01-apr-2015

function Point = readmpt3(filename)

fid         = fopen(filename);
dummy       = fread(fid,14,'ulong');
Point.Tag1  = fread(fid,1,'ulong');
N           = fread(fid,1,'ulong');
M           = fread(fid,1,'ulong'); % number of channels

for lead = 1:M,
    C                   = fread(fid,1,'ulong'); % C == lead
    Point.Tag2(lead)    = fread(fid,1,'ulong'); % should be '0x0000803f'
    Point.Trace(lead,:) = fread(fid,N,'short');
    
    disp(sprintf('%x',ftell(fid)));
    
    Point.Unkown1(lead) = fread(fid,1,'long');
    Point.Name(lead,:)  = fread(fid,66,'char');
    
    disp(sprintf('%x',ftell(fid)));
end

Point.l1    = fread(fid,1,'long');
ndata       = fread(fid,1,'ulong');

return

if ndata ~= 251, keyboard; end

for k = 1:ndata,
    Point.ldata(1,k) = fread(fid,1,'long');
    Point.ldata(2,k) = fread(fid,1,'long');
    Point.rdata(1,k) = fread(fid,1,'float');
    Point.rdata(2,k) = fread(fid,1,'float');
    Point.rdata(3,k) = fread(fid,1,'float');
    Point.rdata(4,k) = fread(fid,1,'float');
    Point.rdata(5,k) = fread(fid,1,'float');
    Point.rdata(6,k) = fread(fid,1,'float');
end

ndata = fread(fid,1,'ulong');

if ndata ~= 251, keyboard; end

for k = 1:ndata,
    Point.ldata(3,k)    = fread(fid,1,'long');
    Point.ldata(4,k)    = fread(fid,1,'long');
    Point.rdata(7,k)    = fread(fid,1,'float');
    Point.rdata(8,k)    = fread(fid,1,'float');
    Point.rdata(9,k)    = fread(fid,1,'float');
    Point.rdata(10,k)   = fread(fid,1,'float');
    Point.rdata(11,k)   = fread(fid,1,'float');
    Point.rdata(12,k)   = fread(fid,1,'float');
end

offset      = ftell(fid);
Point.Tail  = fread(fid);

fclose(fid);