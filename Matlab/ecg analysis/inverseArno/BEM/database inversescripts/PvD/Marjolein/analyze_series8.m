%% MRI
% A MRI of the heart is analyzed to extract the ratio between gray and white brain
% tissue in the head of the subject.
% The following steps are taken:
% * Access DICOM data to create 3D data set
% * Visualize data to explore in 3D
% * Segment White and Gray Brain Tissue using morphology
% * Measure Results

%% Data Access

% filename convention used in image series
% [filenames,dirname] = uigetfile('*.*','MultiSelect','on');

% first filename in series
fname = [dirname filenames{1}];

% examine file header
info = dicominfo(fname)

% extract size info from metadata
voxel_size = double([info.PixelSpacing; info.SliceThickness]')

% extract normals info from metadata
normX=info.ImageOrientationPatient(1:3);
normY=info.ImageOrientationPatient(4:6);
normZ=cross(normX,normY);

for x=0:info.Rows-1
	for y=0:info.Columns-1
		x=double(x);y=double(y);
		X0(x+1,y+1)=double(voxel_size(1)*normX(1)*x+voxel_size(2)*normY(1)*y);
		Y0(x+1,y+1)=double(voxel_size(1)*normX(2)*x+voxel_size(2)*normY(2)*y);
		Z0(x+1,y+1)=double(voxel_size(1)*normX(3)*x+voxel_size(2)*normY(3)*y);		
	end
end

% [X,Y,Z]=meshgrid(

%%
%read slice images; populate XYZ matrix
hWaitBar = waitbar(0,'Reading DICOM files');
D = uint16(zeros(info.Rows,info.Columns,length(filenames)));
X=zeros(info.Rows,info.Columns,length(filenames));Y=X;Z=X;
for i=length(filenames):-1:1
	fname = [dirname filenames{i}];
	Dinfo = dicominfo(fname);
	DP(i,:)=Dinfo.ImagePositionPatient;
	D(:,:,i) = uint16(dicomread(fname));
	m(i) = max(max(D(:,:,i)));
	X(:,:,i)=DP(i,1)+X0;
	Y(:,:,i)=DP(i,2)+Y0;
	Z(:,:,i)=DP(i,3)+Z0;
	waitbar((length(filenames)-i)/length(filenames))
end
delete(hWaitBar)
k=1;
dist=norm3d(DP(length(filenames),:)-DP(1,:));
for i=length(filenames):-1:2
	XI(:,:,i)=double(i)*normZ(1)+X0;
	YI(:,:,i)=double(i)*normZ(2)+Y0;
	ZI(:,:,i)=double(i)*normZ(3)+Z0;
end

DI=newtfit(X,Y,Z,D,XI,YI,ZI);
whos D

%% Visualization

% explore image data using Image Viewer GUI tool
i = 15;  %middle slice
im = D(:,:,i);
max_level = double(max(im(:))); 
imtool(im,[0 max_level])

%%
% explore 3D volume (new Slice-O-Matic viewer) 
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=764

sliceomatic(double(D),0:255,0:255,[0:size(D,3)]*voxel_size(3))
daspect(1./voxel_size)
hSlico2 = gcf;

% %%
% %reorient data for easier interpretation (stand patient up)
% 
% D = permute(D,[3 2 1]);
% voxel_size = voxel_size([1 3 2]);
% for i=1:3
%   D1 = flipdim(D,i);
% end
% whos D
% 
% %explore rotated 3D volume (new Slice-O-Matic viewer) 
% sliceomatic(double(D1))
% daspect(1./voxel_size)
% hSlico2 = gcf;
% %set(hSlico2,'position',[455 63 560 420])

%%
%intensity distribution also useful (more custom graphics)

max_level = double(max(im(:))); 
my_map = jet(max_level);
fig2 = figure; 

%intensity distribution - top 2/3 
subplot(3,1,1:2)
hist(double(im(:)),max_level)
axis([0 max_level 0 900])
title('Distribution')

%color scale - bottom 1/3 
subplot(3,1,3)
imagesc(1:max_level)
colormap(my_map)
xlim([0 max_level])
set(gca,'ytick',[])
ylabel('Color Map')
xlabel('Intensity')
set(fig2,'position',[22 60 560 300],'render','zbuffer')


%% Segmentation

%ignore low levels (backround air, CSF & other soft? tissues)
%using custom GUI tool to select best threshold level 
im = imrotate(squeeze(D(:,:,15)),90);

%custom GUI tool 
thresh_tool(im)

%duplicate original data set for later reference 
D1 = D;

%%
%apply some thresholding rules to ignore certain parts of data
D(D<=34) = 0;       %ignore low levels (CSF & air)
%D(D>=85) = 0;      %ignore high levels (skull & other hard? tissues)
%D(:,:,1:60) = 0;    %ignore spatially low positions (below brain mass)
figure;subplot(2,2,1);imagesc(squeeze(D(:,:,15)));title('slide 30 after thresholding D > 34')

%erode away thick layer (dissolve thin surrounding tissues)
blk = ones([9 9 9]);
D = imerode(D,blk);
subplot(2,2,2);imagesc(squeeze(D(:,:,15)));title('slide 30 after eroding thin tissues')

%isolate brain mass (bwlabeln)
lev = graythresh(double(im)/max_level) * max_level;
bw = (D>=lev+40);
L = bwlabeln(bw); 

%connected region properties - how many, how big?
stats = regionprops(L)

%remove smaller scraps
BW2 = bwareaopen(bw,200,6);
D(~BW2) = 0;
subplot(2,2,3);imagesc(squeeze(D(:,:,15)));title('slide 30 after removing small objects')

%grow back main region (brian mass) 
D = uint16(imdilate(BW2,blk)).*D1;
subplot(2,2,4);imagesc(squeeze(D(:,:,15)));title('slide 30 segmented contents of the heart')

%%
figure
gBrain = patch(isosurface(smooth3(double(D),'box',5),0.1));
set(gBrain,'FaceColor','red','EdgeColor','none');
daspect(1./voxel_size)
view(3); axis tight
camlight 
lighting gouraud

