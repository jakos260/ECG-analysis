function DATA = readGeomPeacsAma(DATA,aAmaName,vAmaName)

% modeldir = [dirname  '/model/'];
% ecgdir = [dirname  '/ecgs/'];
% abeatsdir = [dirname  '/atrial_beats/beat1/'];
% vbeatsdir = [dirname  '/ventricular_beats/beat1/'];

wct = DATA.GEOM.wct;

if ~isempty(vAmaName)
    
    AV = 40 * loadmat(vAmaName);
    
    n = length(DATA.GEOM.thorax.VER);
    DATA.VENTR.THORAX = AV(1:n,:);
    
    DATA.VENTR.LCAV= (AV(n+1: n + length(DATA.GEOM.lcav.VER),:));
    n = n + length(DATA.GEOM.lcav.VER);
    
    DATA.VENTR.RCAV= (AV(n+1: n + length(DATA.GEOM.rcav.VER),:));
    n = n + length(DATA.GEOM.rcav.VER);
    
    DATA.VENTR.RLUNG= (AV(n+1: n + length(DATA.GEOM.rlung.VER),:));
    n = n + length(DATA.GEOM.rlung.VER);
    
    DATA.VENTR.LLUNG= (AV(n+1: n + length(DATA.GEOM.llung.VER),:));
    n = n + length(DATA.GEOM.llung.VER);
    if isfield(DATA.GEOM,'ribcage') && size(AV,1) >  n + length(DATA.GEOM.ribcage.VER)
        DATA.VENTR.RIBCAGE= (AV(n+1: n + length(DATA.GEOM.ribcage.VER),:));
        n = n + length(DATA.GEOM.ribcage.VER);
    end
    if isfield(DATA.GEOM,'fatpad_1') && size(AV,1) >  n + length(DATA.GEOM.fatpad_1.VER)
        DATA.VENTR.FATPAD_1= (AV(n+1: n + length(DATA.GEOM.fatpad_1.VER),:));
        n = n + length(DATA.GEOM.fatpad_1.VER);
    end

    if isfield(DATA.GEOM,'fatpad_2') && size(AV,1) >  n + length(DATA.GEOM.fatpad_2.VER)
        DATA.VENTR.FATPAD_2= (AV(n+1: n + length(DATA.GEOM.fatpad_2.VER),:));
        n = n + length(DATA.GEOM.fatpad_2.VER);
    end

    
    DATA.VENTR.VENTRICLES = zeromean(AV(size(AV,1)-length(DATA.GEOM.ventr.VER)+1:end,:));
    Awct = calcAwct(DATA.VENTR.THORAX,wct);
    DATA.VENTR.THORAX = doWCT(DATA.VENTR.THORAX, Awct);
    % DATA.VENTR.ATRIA = doWCT(DATA.VENTR.ATRIA,   Awct);
    DATA.VENTR.VENTRICLES = doWCT(DATA.VENTR.VENTRICLES,   Awct);
    DATA.VENTR.RLUNG = doWCT(DATA.VENTR.RLUNG,   Awct);
    DATA.VENTR.LLUNG = doWCT(DATA.VENTR.LLUNG,   Awct);
    DATA.VENTR.RCAV = doWCT(DATA.VENTR.RCAV,   Awct);
    DATA.VENTR.LCAV = doWCT(DATA.VENTR.LCAV,   Awct);
    if isfield(DATA.VENTR,'RIBCAGE')
        DATA.VENTR.RIBCAGE = doWCT(DATA.VENTR.RIBCAGE,   Awct);
    end
    if isfield(DATA.VENTR,'FATPAD_1')
        DATA.VENTR.FATPAD_1 = doWCT(DATA.VENTR.FATPAD_1,   Awct);
    end
    if isfield(DATA.VENTR,'FATPAD_2')
        DATA.VENTR.FATPAD_2 = doWCT(DATA.VENTR.FATPAD_2,   Awct);
    end
end


if ~isempty(aAmaName)
    
    AA = 40 * loadmat(aAmaName);
    
    n = length(DATA.GEOM.thorax.VER);
    DATA.ATRIA.THORAX = AA(1:n,:);
    DATA.ATRIA.LCAV= (AA(n+1: n + length(DATA.GEOM.lcav.VER),:));
    n = n + length(DATA.GEOM.lcav.VER);
    
    DATA.ATRIA.RCAV= (AA(n+1: n + length(DATA.GEOM.rcav.VER),:));
    n = n + length(DATA.GEOM.rcav.VER);
    
    DATA.ATRIA.RLUNG= (AA(n+1: n + length(DATA.GEOM.rlung.VER),:));
    n = n + length(DATA.GEOM.rlung.VER);
    
    DATA.ATRIA.LLUNG= (AA(n+1: n + length(DATA.GEOM.llung.VER),:));
    n = n + length(DATA.GEOM.llung.VER);
    
    if isfield(DATA.GEOM,'ribcage') && size(AA,1) >  n + length(DATA.GEOM.ribcage.VER)
        DATA.ATRIA.RIBCAGE= (AA(n+1: n + length(DATA.GEOM.ribcage.VER),:));
        n = n + length(DATA.GEOM.ribcage.VER);
    end
    if isfield(DATA.GEOM,'fatpad_1') && size(AA,1) >  n + length(DATA.GEOM.fatpad_1.VER)
        DATA.ATRIA.FATPAD_1= (AA(n+1: n + length(DATA.GEOM.fatpad_1.VER),:));
        n = n + length(DATA.GEOM.fatpad_1.VER);
    end

    if isfield(DATA.GEOM,'fatpad_2') && size(AA,1) >  n + length(DATA.GEOM.fatpad_2.VER)
        DATA.ATRIA.FATPAD_2= (AA(n+1: n + length(DATA.GEOM.fatpad_2.VER),:));
        n = n + length(DATA.GEOM.fatpad_2.VER);
    end
    
    DATA.ATRIA.ATRIA = (AA(size(AA,1)-length(DATA.GEOM.atria.VER)+1:end,:));
    
    
    Awct = calcAwct(DATA.ATRIA.THORAX,wct);
    DATA.ATRIA.THORAX = doWCT(DATA.ATRIA.THORAX, Awct);
    DATA.ATRIA.ATRIA = doWCT(DATA.ATRIA.ATRIA,   Awct);
    % DATA.ATRIA.VENTRICLES = doWCT(DATA.ATRIA.VENTRICLES,   Awct);
    DATA.ATRIA.RLUNG = doWCT(DATA.ATRIA.RLUNG,   Awct);
    DATA.ATRIA.LLUNG = doWCT(DATA.ATRIA.LLUNG,   Awct);
    DATA.ATRIA.RCAV = doWCT(DATA.ATRIA.RCAV,   Awct);
    DATA.ATRIA.LCAV = doWCT(DATA.ATRIA.LCAV,   Awct);
    if isfield(DATA.ATRIA,'RIBCAGE')
        DATA.ATRIA.RIBCAGE = doWCT(DATA.ATRIA.RIBCAGE,   Awct);
    end
    if isfield(DATA,'ATRIA.FATPAD_1')
        DATA.ATRIA.FATPAD_1 = doWCT(DATA.ATRIA.FATPAD_1,   Awct);
    end
    if isfield(DATA,'ATRIA.FATPAD_2')
        DATA.ATRIA.FATPAD_2 = doWCT(DATA.ATRIA.FATPAD_2,   Awct);
    end
    
end












%%=======================================================================

function A = doWCT(Ain,Awct)

if ~isempty(Ain)
    A = Ain - ones(size(Ain,1),1) * Awct;
else
    A= Ain;
end

function Awct = calcAwct(AthorsoIn,wct)

Awct = mean(AthorsoIn(wct,:));

%%=======================================================================
function [M, extraresult]=loadmat(name)

% LOADMAT	Load a MFBF matrix file.
%
%		Usage: m         = loadmat('file');
%		   or  [m,extra] = loadmat('file');
%
%		LOADMAT('file') returns the matrix stored in 'file' and
%		the extra information stored at the bottom of that file.
%		LOADMAT works for binary as well as asci matrix files.
%
%		See also SAVEMAT.
%
%		Thom Oostendorp, MF&BF University of Nijmegen, the Netherlands
% 20060816; echo (S) switched off

f=fopen(name);
if (f==-1)
    fprintf('\nCannot open %s\n\n', name);
    M=0;
    extraresult='';
    return;
end

[N,nr]=fscanf(f,'%d',2);
if (nr~=2)
    fclose(f);
    f=fopen(name);
    [magic ,nr]=fread(f,8,'char');
    if (char(magic')==';;mbfmat')
        fread(f,1,'char');
        hs=fread(f,1,'long');
        fread(f,1,'char');
        fread(f,1,'char');
        fread(f,1,'char');
        N=fread(f,2,'long');
        M=fread(f,[N(2),N(1)],'double');
    else
        fclose(f);
        f=fopen(name);
        N=fread(f,2,'long');
        M=fread(f,[N(2),N(1)],'float');
    end
else
    M=fscanf(f,'%f',[N(2) N(1)]);
end
[extra,nextra]=fread(f,1000,'char');
fclose(f);
extra = char(extra);

if ~all(isspace(extra))
    %   S=sprintf('%s contains the following extra information:\n', name);
    %   disp(S);
    %   disp(extra');
else
    extra=[];
end
M=M';
extraresult=extra;

%% zeromean.m
% shift potentials to zero mean
% function PSI=zeromean(PHI)
function PSI=zeromean(PHI)
[nlds ntims]=size(PHI);
% shift to zero mean
PSI=PHI-ones(nlds,1)*mean(PHI);
