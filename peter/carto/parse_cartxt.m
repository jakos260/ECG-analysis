function DATA = parse_cartxt(pathname, filename)

mydir=pwd;
meshname = [filename(1:end-8) '.mesh'];
[DATA.VER,DATA.ITRI]=meshparse(pathname,meshname);
cd(mydir);
if isempty(DATA.VER)
    DATA=[];
    return;
end
fid = fopen(fullfile(pathname,filename));
if fid>0 
    str = fgetl(fid);
    %                             1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
    CAR = textscan(fid, '%c %d %d %d %f %f %f %f %f %f %f %f %d %d %d %d %d %c %d %d %d');
    fclose(fid);
else
    warning(['could not open ' filename] );
    DATA =[];
    return 
end
DATA.name = filename(1:end-8);
DATA.ipoint = cell2mat(CAR(3));
DATA.points = cell2mat(CAR(5:7));
DATA.data = [cell2mat(CAR(8:12)) double(cell2mat(CAR(13:17)))];
DATA.annotated = cell2mat(CAR(18));

keep=[];
use = zeros(length(DATA.VER(:,1)),1);
for i=1:length(DATA.points(:,1))
    d = norm3d([DATA.points(i,1) - DATA.VER(:,1) DATA.points(i,2) - DATA.VER(:,2) DATA.points(i,3) - DATA.VER(:,3)]);
    
    di= find(d==min(d(use==0)));di=di(1);
    use(di)=1;
    keep = [keep; [di i]];
end
if ~isempty(keep)
    DATA.mapping = keep;
    DATA.T=intripol(DATA.VER,DATA.ITRI,DATA.mapping(:,1));

    for i=1:length(DATA.ipoint)
        [ECG,DATA.channels] = parseCartoECG(fullfile(pathname,[ filename(1:end-8) '_P' num2str(DATA.ipoint(i)) '_ECG_Export.txt']));
        DATA.ECG_all{i} = ECG;
        DATA.ECGM1(i,:) = ECG(1,:);
        DATA.ECGM2(i,:) = ECG(2,:);    
        DATA.ECGM3(i,:) = ECG(2,:);    
        DATA.ECGM4(i,:) = ECG(2,:);    
    end
end
