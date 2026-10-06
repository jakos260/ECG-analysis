function DATA = readGeomPeacsModel_prof(dirname,subject,varargin)

% modeldir = [dirname  '/model/'];
% ecgdir = [dirname  '/ecgs/'];
% abeatsdir = [dirname  '/atrial_beats/beat1/'];
% vbeatsdir = [dirname  '/ventricular_beats/beat1/'];

modeldir = fullfile(dirname ,subject);
modeldir = [modeldir '_'];
if exist([modeldir  'atria.adj2d'])    
    [DATA.ATRIA.geom.VER,DATA.ATRIA.geom.ITRI] = loadtri([modeldir  'atria.tri']);
    [DATA.ATRIA.geom.NORMV,DATA.ATRIA.geom.NORMT] = trinormals(DATA.ATRIA.geom.VER,DATA.ATRIA.geom.ITRI);

    DATA.ATRIA.ADJsurf = loadmat([modeldir  'atria.adj2d']);
    DATA.ATRIA.ADJ3D = loadmat([modeldir  'atria.adj3d']);
    DATA.ATRIA.DISTsurf = loadmat([modeldir  'atria.dst2d']);
    DATA.ATRIA.DIST3D = loadmat([modeldir  'ATRIA.dst3d']);
    DATA.ATRIA.ADJANIS = loadmat([modeldir  'ATRIA.adjanis']);
    DATA.ATRIA.DISTANIS = loadmat([modeldir  'ATRIA.dstanis']);
    [DATA.GEOM.atria.VER, DATA.GEOM.atria.ITRI] = loadtri([modeldir  'atria.tri']);
    DATA.GEOM.atria.geomtyp = loadmat([modeldir  'atria.typ']);
end
[DATA.VENTR.geom.VER, DATA.VENTR.geom.ITRI] = loadtri([modeldir  'ventricles.tri']);
[DATA.VENTR.geom.NORMV,DATA.VENTR.geom.NORMT] = trinormals(DATA.VENTR.geom.VER,DATA.VENTR.geom.ITRI);
DATA.VENTR.ADJsurf = loadmat([modeldir  'ventricles.adj2d']);
DATA.VENTR.ADJ3D = loadmat([modeldir  'ventricles.adj3d']);
DATA.VENTR.DISTsurf = loadmat([modeldir  'ventricles.dst2d']);
DATA.VENTR.DIST3D = loadmat([modeldir  'ventricles.dst3d']);
DATA.VENTR.ADJANIS = loadmat([modeldir  'ventricles.adjanis']);
DATA.VENTR.DISTANIS = loadmat([modeldir  'ventricles.dstanis']);
DATA.VENTR.HEARTDIST = loadmat([modeldir  'ventricles.heartdist']);

[DATA.GEOM.ventr.VER,DATA.GEOM.ventr.ITRI] = loadtri([modeldir  'ventricles.tri']);
[DATA.GEOM.ventr.NORMV,DATA.GEOM.ventr.NORMT] = trinormals(DATA.GEOM.ventr.VER,DATA.GEOM.ventr.ITRI);
DATA.GEOM.ventr.geomtyp = loadmat([modeldir  'ventricles.typ']);
DATA.GEOM.ventr.scar = loadmat([modeldir  'ventricles.scar']);
DATA.GEOM.ventr.walls= loadmat([modeldir  'ventricles.walls']);
DATA.GEOM.ventr.segments= loadmat([modeldir  'ventricles.segments']);



[DATA.GEOM.lcav.VER, DATA.GEOM.lcav.ITRI]      = loadtri([modeldir  'lcav.tri']);
[DATA.GEOM.rcav.VER, DATA.GEOM.rcav.ITRI]      = loadtri([modeldir  'rcav.tri']);
if exist([modeldir  'llung.tri'])
    [DATA.GEOM.llung.VER,DATA.GEOM.llung.ITRI]     = loadtri([modeldir  'llung.tri']);
    
end
if exist( [modeldir  'rlung.tri'])
    [DATA.GEOM.rlung.VER,DATA.GEOM.rlung.ITRI]     = loadtri([modeldir  'rlung.tri']);
end
[DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI]   = loadtri([modeldir  'thorax.tri']);
[DATA.GEOM.thorax.NORMV,DATA.GEOM.thorax.NORMT] = trinormals(DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI);

[DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI]   = loadtri([modeldir  'thorax.tri']);
if exist([modeldir  'liver.tri'])
    [DATA.GEOM.liver.VER,DATA.GEOM.liver.ITRI]   = loadtri([modeldir  'liver.tri']);
end

if exist([modeldir  'ribcage.tri'])
    [DATA.GEOM.ribcage.VER,DATA.GEOM.ribcage.ITRI]   = loadtri([modeldir  'ribcage.tri']);
end
if exist([modeldir  'fatpad_1.tri'])
    [DATA.GEOM.fatpad_1.VER,DATA.GEOM.fatpad_1.ITRI]   = loadtri([modeldir  'fatpad_1.tri']);
end
if exist([modeldir  'fatpad_2.tri'])
    [DATA.GEOM.fatpad_2.VER,DATA.GEOM.fatpad_2.ITRI]   = loadtri([modeldir  'fatpad_2.tri']);
end



if nargin ==3
    if ischar(varargin{end})
        dd = dir([varargin{end} '/*.lead']);
        for i=1:length(dd)
            if ~isempty(strfind(dd(i).name,'standard12lead'))
                break;
            end
        end
        ver = loadmat(fullfile(varargin{end}, dd(i).name));
        ver = ver(:,2:4);
        meanTh = mean(DATA.GEOM.thorax.VER);
        wct=[];
        for i=1:3
            
            TRIS= linetris(DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI,ver(i,:),[meanTh(1) meanTh(2) ver(i,3)]);
            if isempty(TRIS)
                TRIS= linetris(DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI,ver(i,:),meanTh);
                TRIS(TRIS(:,2)<0,:) = [];
            end
            TRIS(abs(TRIS(:,end))> min(abs(TRIS(:,end))),:)=[];
            %         itri = DATA.GEOM.thorax.ITRI(TRIS(1),:);
            wct = [wct DATA.GEOM.thorax.ITRI(TRIS(1),1)];
        end
    else
        wct = varargin{end};
    end
else
    dd = dir([dirname '/*.lead']);
    if ~isempty(dd)
        for i=1:length(dd)
            if ~isempty(strfind(dd(i).name,'standard12lead'))
                break;
            end
        end
        ver = loadmat(fullfile(dirname, dd(i).name));
        ver = ver(:,2:4);
        meanTh = mean(DATA.GEOM.thorax.VER);
        wct=[];
        for i=1:3
            TRIS= linetris(DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI,ver(i,:),[meanTh(1) meanTh(2) ver(i,3)]);
            TRIS(abs(TRIS(:,end))> min(abs(TRIS(:,end))),:)=[];
            %         itri = DATA.GEOM.thorax.ITRI(TRIS(1),:);
            wct = [wct DATA.GEOM.thorax.ITRI(TRIS(1),1)];
        end
    end
end

if exist('wct')
    DATA = splitAMA(modeldir, DATA,wct,1);   
    DATA = splitAMA(modeldir,DATA,wct,0);    
end

%%
if nargin == 3 && ischar(varargin{end})
    dd = dir([varargin{end} '/*.lead']);
    for i=1:length(dd)
        if ~isempty(strfind(dd(i).name,'standard12lead'))
            break;
        end
    end
    ver = loadmat(fullfile(varargin{end}, dd(i).name));
    ver = ver(:,2:4);
    name = 'standard12lead';
    iVer =zeros(length(ver),1);
    for i= 1:length(ver)
        d = norm3d(bsxfun(@minus, DATA.GEOM.thorax.VER, ver(i,:)));
        iVer(i) = find(d==min(d));
    end
    
    
    
    %     [AMA_A,AMA_V] = getAMA(DATA,ver);
    %     eval(['DATA.VENTR.' name '=AMA_V;']);
    %     eval(['DATA.ATRIA.' name '=AMA_A;']);
    
    eval(['DATA.VENTR.LEADPOS.' name '=ver;']);
    eval(['DATA.VENTR.ILEADPOS.' name '=iVer;']);

    eval(['DATA.ATRIA.LEADPOS.' name '=ver;']);
    eval(['DATA.ATRIA.ILEADPOS.' name '=iVer;']);
else
    dd = dir([dirname '/*.lead']);
    
    for i=1:length(dd)
        ver = loadmat(fullfile(dirname, dd(i).name));
        ver = ver(:,2:4);
        
        name = dd(i).name(1:end-5);
        name = strrep(name,' ','');
        name = strrep(name,subject,'');
        name = ['A' name];
        name = strrep(name,'+','_');

        [AMA_A,AMA_V] = getAMA(DATA,ver);
        eval(['DATA.VENTR.' name '=AMA_V;']);
        eval(['DATA.VENTR.LEADPOS.' name '=ver;']);
        if isfield(DATA,'ATRIA')
            eval(['DATA.ATRIA.' name '=AMA_A;']);               
            eval(['DATA.ATRIA.LEADPOS.' name '=ver;']);
        end
    end
end
%%
function DATA = splitAMA(modeldir,DATA,wct,usetriag)
extentionA = 'aedl';
extentionV = 'vedl';
if ( usetriag )
    extentionA = ['t' extentionA];
    extentionV = ['t' extentionV];
end
        
if exist([modeldir 'all.' extentionA],'file')
    AA = 40 * loadmat([modeldir 'all.' extentionA]);
elseif exist([modeldir 'thorax.' extentionA],'file')
    AA = 40 * loadmat([modeldir 'thorax.' extentionA]);
else
    AA=[];
end
if exist([modeldir 'all.' extentionV],'file')
    AV = 40*loadmat([modeldir 'all.' extentionV]);
elseif exist([modeldir 'thorax.' extentionV],'file')
    AV = 40*loadmat([modeldir 'thorax.' extentionV]);
else
    AV=[];
end

if ~isempty(AV)
    n = length(DATA.GEOM.thorax.VER);    
    AV = doWCT(AV,calcAwct(AV(1:n,:),wct));
    if usetriag
        DATA.VENTR.tTHORAX = AV(1:n,:);
    else
        DATA.VENTR.THORAX = AV(1:n,:);
    end
    if size(AV,1)>n
        if usetriag
            DATA.VENTR.tRCAV= (AV(n+1: n + length(DATA.GEOM.rcav.VER),:));
        else
            DATA.VENTR.RCAV= (AV(n+1: n + length(DATA.GEOM.rcav.VER),:));
        end
        n = n + length(DATA.GEOM.rcav.VER);
    end
    if size(AV,1)>n
        if usetriag
            DATA.VENTR.tLCAV= (AV(n+1: n + length(DATA.GEOM.lcav.VER),:));
        else
            DATA.VENTR.LCAV= (AV(n+1: n + length(DATA.GEOM.lcav.VER),:));
        end
        n = n + length(DATA.GEOM.lcav.VER);
    end


    if size(AV,1)>n && isfield(DATA.GEOM,'rlung')
        if usetriag
            DATA.VENTR.tRLUNG= (AV(n+1: n + length(DATA.GEOM.rlung.VER),:));
        else
            DATA.VENTR.RLUNG= (AV(n+1: n + length(DATA.GEOM.rlung.VER),:));
        end
        n = n + length(DATA.GEOM.rlung.VER);
    end
    if size(AV,1)>n && isfield(DATA.GEOM,'llung')
        if usetriag
            DATA.VENTR.tLLUNG= (AV(n+1: n + length(DATA.GEOM.llung.VER),:));
        else
            DATA.VENTR.LLUNG= (AV(n+1: n + length(DATA.GEOM.llung.VER),:));
        end
        n = n + length(DATA.GEOM.llung.VER);

        if isfield(DATA.GEOM,'liver') 
            if usetriag
                DATA.VENTR.tLIVER= (AV(n+1: n + length(DATA.GEOM.liver.VER),:));
            else
                DATA.VENTR.LIVER= (AV(n+1: n + length(DATA.GEOM.liver.VER),:));
            end
            n = n + length(DATA.GEOM.liver.VER);
        end

        if isfield(DATA.GEOM,'ribcage')
            if usetriag
                DATA.VENTR.tRIBCAGE= (AV(n+1: n + length(DATA.GEOM.ribcage.VER),:));
            else
                DATA.VENTR.RIBCAGE= (AV(n+1: n + length(DATA.GEOM.ribcage.VER),:));
            end
            n = n + length(DATA.GEOM.ribcage.VER);
        end
        if isfield(DATA.GEOM,'fatpad_1') 
            if usetriag
                DATA.VENTR.tFATPAD_1= (AV(n+1: n + length(DATA.GEOM.fatpad_1.VER),:));
            else
                DATA.VENTR.FATPAD_1= (AV(n+1: n + length(DATA.GEOM.fatpad_1.VER),:));
            end
            n = n + length(DATA.GEOM.fatpad_1.VER);
        end

        if isfield(DATA.GEOM,'fatpad_2') 
            if usetriag
                DATA.VENTR.tFATPAD_2= (AV(n+1: n + length(DATA.GEOM.fatpad_2.VER),:));
            else
                DATA.VENTR.FATPAD_2 = (AV(n+1: n + length(DATA.GEOM.fatpad_2.VER),:));
            end
            n = n + length(DATA.GEOM.fatpad_2.VER);
        end
    end
    if  size(AV,1)-n == length(DATA.GEOM.ventr.VER)   %size(AV,1) > n  
        if usetriag
            DATA.VENTR.tVENTRICLES = (AV(end-length(DATA.GEOM.ventr.VER)+1:end,:));
        else
            DATA.VENTR.VENTRICLES = (AV(end-length(DATA.GEOM.ventr.VER)+1:end,:));
        end
    end
end
n = length(DATA.GEOM.thorax.VER);
if size(AA,1) > n
    n = length(DATA.GEOM.thorax.VER);
    AA = doWCT(AA,calcAwct(AA(1:n,:),wct));
    if usetriag
        DATA.ATRIA.tTHORAX = AA(1:n,:);
    else
        DATA.ATRIA.THORAX = AA(1:n,:);
    end

    if usetriag
        DATA.ATRIA.tLCAV= (AA(n+1: n + length(DATA.GEOM.lcav.VER),:));
    else
        DATA.ATRIA.LCAV= (AA(n+1: n + length(DATA.GEOM.lcav.VER),:));
    end

    n = n + length(DATA.GEOM.lcav.VER);
    if usetriag
        DATA.ATRIA.tRCAV= (AA(n+1: n + length(DATA.GEOM.rcav.VER),:));
    else
        DATA.ATRIA.RCAV= (AA(n+1: n + length(DATA.GEOM.rcav.VER),:));
    end
    n = n + length(DATA.GEOM.rcav.VER);

    if isfield(DATA.GEOM,'rlung')
        if usetriag
            DATA.ATRIA.tRLUNG= (AA(n+1: n + length(DATA.GEOM.rlung.VER),:));
        else
            DATA.ATRIA.RLUNG= (AA(n+1: n + length(DATA.GEOM.rlung.VER),:));
        end
        n = n + length(DATA.GEOM.rlung.VER);
    end
    if isfield(DATA.GEOM,'llung')                    
        if usetriag
            DATA.ATRIA.tLLUNG= (AA(n+1: n + length(DATA.GEOM.llung.VER),:));
        else
            DATA.ATRIA.LLUNG= (AA(n+1: n + length(DATA.GEOM.llung.VER),:));
        end
        n = n + length(DATA.GEOM.llung.VER);
    end
    if isfield(DATA.GEOM,'liver')
        if usetriag
            DATA.VENTR.tLIVER= (AV(n+1: n + length(DATA.GEOM.liver.VER),:));
        else
            DATA.VENTR.LIVER= (AV(n+1: n + length(DATA.GEOM.liver.VER),:));
        end
        n = n + length(DATA.GEOM.liver.VER);
    end

    if isfield(DATA.GEOM,'ribcage')
        if usetriag
            DATA.ATRIA.tRIBCAGE= (AA(n+1: n + length(DATA.GEOM.ribcage.VER),:));
        else
            DATA.ATRIA.RIBCAGE= (AA(n+1: n + length(DATA.GEOM.ribcage.VER),:));
        end
        n = n + length(DATA.GEOM.ribcage.VER);
    end
    if isfield(DATA.GEOM,'fatpad_1')
        if usetriag
            DATA.ATRIA.tFATPAD_1= (AA(n+1: n + length(DATA.GEOM.fatpad_1.VER),:));
        else
            DATA.ATRIA.FATPAD_1= (AA(n+1: n + length(DATA.GEOM.fatpad_1.VER),:));
        end
        n = n + length(DATA.GEOM.fatpad_1.VER);
    end
    if isfield(DATA.GEOM,'fatpad_2')
        if usetriag
            DATA.ATRIA.tFATPAD_2= (AA(n+1: n + length(DATA.GEOM.fatpad_2.VER),:));
        else
            DATA.ATRIA.FATPAD_2= (AA(n+1: n + length(DATA.GEOM.fatpad_2.VER),:));
        end
        n = n + length(DATA.GEOM.fatpad_2.VER);
    end
    if usetriag
        DATA.ATRIA.tATRIA = (AA(end-length(DATA.GEOM.atria.VER)+1: end,:));
    else
        DATA.ATRIA.ATRIA = (AA(end-length(DATA.GEOM.atria.VER)+1: end,:));
    end
end


%%
DATA.GEOM.ventr.endoVER=zeros(1,length(DATA.GEOM.ventr.VER));
for i=1:length(DATA.GEOM.ventr.VER)
    if any(DATA.GEOM.lcav.VER(:,1)==DATA.GEOM.ventr.VER(i,1) & ...
            DATA.GEOM.lcav.VER(:,2)==DATA.GEOM.ventr.VER(i,2) & ...
            DATA.GEOM.lcav.VER(:,3)==DATA.GEOM.ventr.VER(i,3))
        DATA.GEOM.ventr.endoVER(i) = 2;
    end
end

for i=1:length(DATA.GEOM.ventr.VER)
    if any(DATA.GEOM.rcav.VER(:,1)==DATA.GEOM.ventr.VER(i,1) & ...
            DATA.GEOM.rcav.VER(:,2)==DATA.GEOM.ventr.VER(i,2) & ...
            DATA.GEOM.rcav.VER(:,3)==DATA.GEOM.ventr.VER(i,3))
        DATA.GEOM.ventr.endoVER(i) = 1;
    end
end

DATA.GEOM.ventr.type = DATA.GEOM.ventr.endoVER;

ADJ = DATA.VENTR.ADJ3D;
if size(ADJ,1) == length(DATA.GEOM.ventr.VER)
    ADJ(DATA.VENTR.ADJsurf>0)=0;
    for i=1:length(DATA.GEOM.ventr.VER)
        lvD = min(ADJ(i,DATA.GEOM.ventr.endoVER==2 & ADJ(i,:)>0 & ADJ(i,:)< 35));
        rvD = min(ADJ(i,DATA.GEOM.ventr.endoVER==1 & ADJ(i,:)>0 & ADJ(i,:)< 35));
        evD = min(ADJ(i,DATA.GEOM.ventr.endoVER==0 & ADJ(i,:)>0 & ADJ(i,:)< 35));
        if isempty(evD) && DATA.GEOM.ventr.endoVER(i) ~= 0
            if DATA.GEOM.ventr.endoVER(i) ==1
                DATA.GEOM.ventr.type(i) = 3;% rvseptum endo
            else
                DATA.GEOM.ventr.type(i) = 4;% lvseptum endo
            end
        elseif isempty(lvD) && DATA.GEOM.ventr.endoVER(i) ~= 2
            if DATA.GEOM.ventr.endoVER(i) == 0
                DATA.GEOM.ventr.type(i) = 5;% rvfree epi
            else
                DATA.GEOM.ventr.type(i) = 6;% rvfree endo
            end
        elseif isempty(rvD) && DATA.GEOM.ventr.endoVER(i) ~= 1
            if DATA.GEOM.ventr.endoVER(i) == 0
                DATA.GEOM.ventr.type(i) = 7;% lvfree epi
            else
                DATA.GEOM.ventr.type(i) = 8;% lvfree endo
            end
        else
            if lvD < rvD * 3 % this is lv
                if DATA.GEOM.ventr.endoVER(i) == 0 %epi
                    DATA.GEOM.ventr.type(i) = 7;% lvfree epi
                elseif DATA.GEOM.ventr.endoVER(i) == 1 %rv endo
                    DATA.GEOM.ventr.type(i) = 3;% rvseptum endo
                else
                    DATA.GEOM.ventr.type(i) = 4;% lvseptum endo
                end
            else
                if DATA.GEOM.ventr.endoVER(i) == 0 %epi
                    DATA.GEOM.ventr.type(i) = 5;% rvfree epi
                elseif DATA.GEOM.ventr.endoVER(i) == 1 %rv endo
                    if isempty(lvD)
                        DATA.GEOM.ventr.type(i) = 6;% rvfree endo
                    else
                        DATA.GEOM.ventr.type(i) = 3;% rvfree endo
                    end
                else
                    if isempty(rvD)
                        DATA.GEOM.ventr.type(i) = 8;% lvfree endo
                    else
                        DATA.GEOM.ventr.type(i) = 4;% lvseptum
                    end
                end
            end
        end
    end
    
end

if isfield(DATA.GEOM,'atria') && isfield(DATA.GEOM.atria, 'VER')
    DATA.GEOM.atria.endoVER=zeros(1,length(DATA.GEOM.atria.VER));
    for i=1:length(DATA.GEOM.atria.VER)
        if any(DATA.GEOM.lcav.VER(:,1)==DATA.GEOM.atria.VER(i,1) & ...
                DATA.GEOM.lcav.VER(:,2)==DATA.GEOM.atria.VER(i,2) & ...
                DATA.GEOM.lcav.VER(:,3)==DATA.GEOM.atria.VER(i,3))
            DATA.GEOM.atria.endoVER(i) = 2;
        end
    end
    for i=1:length(DATA.GEOM.atria.VER)
        if any(DATA.GEOM.rcav.VER(:,1)==DATA.GEOM.atria.VER(i,1) & ...
                DATA.GEOM.rcav.VER(:,2)==DATA.GEOM.atria.VER(i,2) & ...
                DATA.GEOM.rcav.VER(:,3)==DATA.GEOM.atria.VER(i,3))
            DATA.GEOM.atria.endoVER(i) = 1;
        end
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
%function [M, extraresult]=loadmat(name);
function [M, extraresult]=loadmat(name);

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
%   fprintf('\nCannot open %s\n\n', name);
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
    pos = ftell(f);
    [M,count]=fread(f,[N(2),N(1)],'double');    
    if count ~= N(2) * N(1)
        status = fseek(f,pos,'bof');
        M=fread(f,[N(2),N(1)],'float');    
    end
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

%%========================================================================
function [pnt, dhk] = loadtri(fn)

fid = fopen(fn, 'rt');
if fid~=-1
    
    % read the vertex points
    Npnt = fscanf(fid, '%d', 1);
    pnt  = fscanf(fid, '%f', [4, Npnt]);
    pnt  = pnt(2:4,:)';
    
    % if present, read the triangles
    if (~(feof(fid)))
        [Ndhk, count] = fscanf(fid, '%d', 1);
        if (count ~= 0)
            dhk = fscanf(fid, '%d', [4, Ndhk]);
            dhk = dhk(2:4,:)';
        end
    else
        dhk = [];
    end
    fclose(fid);
    
else
    error(['unable to open file: ' fn]);
end

%%
function [AMA_A,AMA_V] = getAMA(DATA,ver)

meanTh = mean(DATA.GEOM.thorax.VER);
if isfield(DATA,'ATRIA') && isfield(DATA.ATRIA,'THORAX')
    AMA_A = zeros(size(ver,1),size(DATA.ATRIA.THORAX,2));
else
    AMA_A = [];
end
AMA_V = zeros(size(ver,1),size(DATA.VENTR.THORAX,2));
for i=1:size(ver,1)
    
    TRIS= linetris(DATA.GEOM.thorax.VER,DATA.GEOM.thorax.ITRI,ver(i,:),[meanTh(1) meanTh(2) ver(i,3)]);
    TRIS(abs(TRIS(:,end))> min(abs(TRIS(:,end))),:)=[];
    itri = DATA.GEOM.thorax.ITRI(TRIS(1),:);
    if isfield(DATA,'ATRIA') && isfield(DATA.ATRIA,'THORAX')
        AMA_A(i,:) = DATA.ATRIA.THORAX(itri(1),:) + ...
            (DATA.ATRIA.THORAX(itri(2),:) - DATA.ATRIA.THORAX(itri(1),:)) * TRIS(3) + ...
            (DATA.ATRIA.THORAX(itri(3),:) - DATA.ATRIA.THORAX(itri(1),:)) * TRIS(4);
    end
    AMA_V(i,:) = DATA.VENTR.THORAX(itri(1),:) + ...
        (DATA.VENTR.THORAX(itri(2),:) - DATA.VENTR.THORAX(itri(1),:)) * TRIS(3) + ...
        (DATA.VENTR.THORAX(itri(3),:) - DATA.VENTR.THORAX(itri(1),:)) * TRIS(4);
    
end
