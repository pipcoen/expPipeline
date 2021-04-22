
%% 2021-02-20 seeing if I can do gradients of cortical layers

% to do:
% was probably a mistake to focus on isocortex
%
%

%%

addpath(genpath('../allenCCF'));
addpath(genpath('../npy-matlab'));

% Structure Tree: a table that explains the labels
st = loadStructureTree('../Allen/structure_tree_safe_2017.csv'); 

% Template Volume: gray-scale "background signal intensity" [AP,DV,ML]
tv = readNPY('../Allen/template_volume_10um.npy'); 

% Annotation Volume: the number at each pixel labels the area [AP,DV,ML]
av = readNPY('../Allen/annotation_volume_10um_by_index.npy'); 

% Allen CCF-bregma transform (estimated from eyeballing Paxinos->CCF)
% [AP,DV,ML]
% this is a function: AllenCCFBregma
bregma = [540,0,570];

%% SKIP THIS SECTION these are the functions that would be called 

% allen_ccf_npx(tv,av,st);

file_save_location = 'C:\Histology\Mouse1'; % where will the probe locations be saved
probe_name = 'test'; % name probe to avoid overwriting

f = allenAtlasBrowser(tv, av, st, file_save_location, probe_name);

plotBrainGrid();

%% example of navigation through the structure tree

idx = 85;
st.structure_id_path(idx)

idx1 = find(st.id == st.parent_structure_id(idx ));
idx2 = find(st.id == st.parent_structure_id(idx1));
idx3 = find(st.id == st.parent_structure_id(idx2));
idx4 = find(st.id == st.parent_structure_id(idx3));
idx5 = find(st.id == st.parent_structure_id(idx4));
idx6 = find(st.id == st.parent_structure_id(idx5));

st.name(idx )
st.name(idx1)
st.name(idx2)
st.name(idx3)
st.name(idx4) % isocortex
st.name(idx5)
st.name(idx6)

%% make a volume with layers (from 1 to 5)

figure; clf
img = imagesc( squeeze( av( bregma(1),:,:) ) );
axis image

cmap = allen_ccf_colormap('2017');
colormap(cmap)

st(strcmp(st.name,'Isocortex'),:) % this is at depth = 5

id_isocortex = st.id(strcmp(st.name,'Isocortex')); % the id of isocortex is 315
ii_isocortex = contains(st.structure_id_path, ['/' num2str(id_isocortex) '/']);

LayerNames = {'layer 1','layer 2','layer 4','layer 5','layer 6'};

LayerVolume = zeros(size(av),'uint8');
for iLayer = 1:5
    ii_layer = contains(st.name,LayerNames{iLayer},'IgnoreCase',true)&ii_isocortex;
    LayerRows = find(ii_layer); % rows of st that have this layer
    for iLayerRow = 1:length(LayerRows)
        iRow = LayerRows(iLayerRow);
        disp(st.name(iRow));
        LayerVolume(av==iRow)=iLayer;
    end
end

figure; clf
img = imagesc( squeeze( LayerVolume( bregma(1),:,:) ) );
axis image

save('LayerVolume','LayerVolume');

%% load LayerVolume

load LayerVolume;
[nx,ny,nz] = size(LayerVolume);

figure; clf
imagesc( squeeze( LayerVolume( bregma(1),:,:) ) );
axis image

figure; clf
imagesc( squeeze( LayerVolume( bregma(1)+200,:,:) ) );
axis image



%% downsample it and smooth it

% the topic in Matlab is 3-D volumetric image processing

% decimate it brutally

DecLayerVolume = double(LayerVolume(1:10:nx,1:10:ny,1:10:nz));


% blur it
H = fspecial3('ellipsoid',[5 5 5]);
volSmooth = imfilter(DecLayerVolume,H,'replicate');


figure; clf
imagesc( squeeze( DecLayerVolume( (bregma(1)+200)/10,:,:) ) );
axis image

figure; clf
imagesc( squeeze( volSmooth( (bregma(1)+200)/10,:,:) ) );
axis image

%% find the gradients

[ DnGradX, DnGradY, DnGradZ] = gradient(volSmooth,100);

iX = (bregma(1)+200)/10;

figure; clf; ax = [];
ax(1) = subplot(3,1,1); imagesc( squeeze( DnGradX(iX,:,:) ) ); axis image
ax(2) = subplot(3,1,2); imagesc( squeeze( DnGradY(iX,:,:) ) ); axis image
ax(3) = subplot(3,1,3); imagesc( squeeze( DnGradZ(iX,:,:) ) ); axis image


figure; clf
imagesc( squeeze( volSmooth(iX,:,:) ) ); hold on
yy = 1:5:ny/10;
zz = 1:5:nz/10;

quiver( squeeze(DnGradY(iX,yy,zz)), squeeze(DnGradZ(iX,yy,zz)) )
axis image

figure; clf
slice(double(DecLayerVolume),bregma(1)/10,ny/20*(1:3),[]);
shading flat

% this upsamples it! yuck!
DnLayerVolume = interp3(double(LayerVolume),4); % reduce by a factor of 2^4 = 16

% for some reason this takes a ton of memory:
figure
slice(LayerVolume,bregma(1),bregma(2),bregma(3));
shading flat

man interp3

% this is crap. rather, see the manual page for interp3
% DnLayerVolume = downsample(LayerVolume,10); % one point every 100 um
% 
% figure; clf
% img = imagesc( squeeze( DnLayerVolume( bregma(1)/10,:,:) ) );
% axis image
% 
% 
% DnGradX







