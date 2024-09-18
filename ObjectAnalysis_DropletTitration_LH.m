clear all

% Specify the directory that contains the extracted files from the last
% step, where you extracted from the raw files obtained from the microscope
sourceDirectory = ...
    './ExtractedStacks_Titration/**/';

% Channels for segmentation
Condensate_SegChannel = 2; % Channel used to detect droplets
Surface_SegChannel = 1; % Channel used to detect surface structures

% Save images of the clusters
ImgSquareExtension = 0; % pixels for cut-out image extension, set 0 for no images
% Which image channels to store in example images
storeImgChannels = [];
numStoreChannels = numel(storeImgChannels);

% Target channels for intensity quantification, applied for all objects
quantChannels = [1,2];
quantBlurSigma = [0,0];

% Parameters for the blurring during image preprocessing
Condensate_segBlurSigma_object = 0.5; % in microns
Condensate_segBlurSigma_BG_removal = 1.5; % in microns
% Number of standard deviations above background for the robust threshold
Condensate_seg_numStdDev = 12.0;

% Blurring and segmentation parameters for the Surface
Surface_segBlurSigma_object = 0.5; % in microns
Surface_segBlurSigma_BG_removal = 1.5; % in microns
Surface_seg_numStdDev = 12.0; % number of standard deviations in robust threshold

% Minimum volumes for condensate and surface objects
Condensate_minVol = 0.02; % cubic microns
Surface_minVol = 0.02; % cubic microns

% end of analysis parameter section, do not change anything else in
% this section, all necessary parameters are listed above


%% --- analysis procedure begins here

% Recursive directory search to find all .mat files we previously saved
% during extraction
listing = rdir([sourceDirectory,'*Image*.mat']);
numFiles = numel(listing);

% Condition index retrieval
condInds = [];
condNames = {};
for ff = 1:numFiles
	thisFilePath = listing(ff).name;
	thisCondInd = load(thisFilePath,'condInd');
	thisCondInd = thisCondInd.condInd;
	condInds = [condInds,thisCondInd];
	thisCondName = load(thisFilePath,'condName');
	thisCondName = thisCondName.condName;
	condNames = [condNames,thisCondName];
end


%% --- analyze image stacks one by one

% Flag to log files that were successfully analyzed
validFileFlag = false(1,numFiles);

% Variables to store properties
numSurf_vec = zeros(1,numFiles);
numCond_vec = zeros(1,numFiles);
cond_intCell = cell(1,numFiles);
drop_intCell = cell(1,numFiles);
pixelSize_vec = zeros(1,numFiles);
zStepSize_vec = zeros(1,numFiles);

% Variable to store the pixel sizes
Condensate_xyVoxelSizeCell = cell(1,numFiles);
Condensate_zVoxelSizeCell = cell(1,numFiles);
Surface_xyVoxelSizeCell = cell(1,numFiles);
Surface_zVoxelSizeCell = cell(1,numFiles);

% Variables to store properties of objects inside nuclei
Condensate_volCell = cell(1,numFiles);
Condensate_solCell = cell(1,numFiles);
Condensate_eloCell = cell(1,numFiles);
Condensate_intCell = cell(1,numFiles);
Condensate_centCell = cell(1,numFiles);
Condensate_imgCell = cell(1,numFiles);
Condensate_Surface_distCell = cell(1,numFiles);

Surface_volCell = cell(1,numFiles);
Surface_solCell = cell(1,numFiles);
Surface_eloCell = cell(1,numFiles);
Surface_intCell = cell(1,numFiles);
Surface_centCell = cell(1,numFiles);
Surface_imgCell = cell(1,numFiles);
Surface_Condensate_distCell = cell(1,numFiles);

numQuantChannels = numel(quantChannels);

for ff = 1:numFiles
	
	fprintf('Processing file %d of %d\n',ff,numFiles)
	
	thisCondInd = condInds(ff);	
	thisFilePath = listing(ff).name;
	
	loadStruct = load(thisFilePath,...
		'imgStack','imgSize','pixelSize','zStepSize');
	imgStack = loadStruct.imgStack;
	imgSize = loadStruct.imgSize;
	pixelSize = loadStruct.pixelSize;
	zStepSize = loadStruct.zStepSize;

    pixelSize_vec(ff) = pixelSize;
    zStepSize_vec(ff) = zStepSize;

	% Surface segmentation
	segImg_surf = imgStack{Surface_SegChannel};
    if Surface_segBlurSigma_object>0
        segImg_surf = ...
            + imgaussfilt(segImg_surf,Surface_segBlurSigma_object./pixelSize) ...
            - imgaussfilt(segImg_surf,Surface_segBlurSigma_BG_removal./pixelSize);
    else
        segImg_surf = ...
            + segImg_surf ...
            - imgaussfilt(segImg_surf,Surface_segBlurSigma_BG_removal./pixelSize);
    end

    % Condensate segmentation
    segImg_cond = imgStack{Condensate_SegChannel};
    if Condensate_segBlurSigma_object>0
        segImg_cond = ...
            + imgaussfilt(segImg_cond,Surface_segBlurSigma_object./pixelSize) ...
            - imgaussfilt(segImg_cond,Surface_segBlurSigma_BG_removal./pixelSize);
    else
        segImg_cond = ...
            + segImg_cond ...
            - imgaussfilt(segImg_cond,Condensate_segBlurSigma_BG_removal./pixelSize);
    end

    seg_intensities_surf = segImg_surf(:);
    seg_mean_surf = mean(seg_intensities_surf);
    seg_std_surf = std(seg_intensities_surf);
    SurfaceSegMask = segImg_surf>(seg_mean_surf+Surface_seg_numStdDev.*seg_std_surf);

    seg_intensities_cond = segImg_cond(:);
    seg_mean_cond = mean(seg_intensities_cond);
    seg_std_cond = std(seg_intensities_cond);
    CondensateSegMask = segImg_cond>(seg_mean_cond+Condensate_seg_numStdDev.*seg_std_cond);    
    
    % display raw image, preprocessed image, and segmented image for both
    % channels

	subplot(4,3,1)
	imagesc(squeeze(imgStack{Surface_SegChannel}(:,:,ceil(imgSize(3)./2))))
	axis tight equal
    title('Raw image')
    ylabel('Surface, central section')
	
	subplot(4,3,2)
	imagesc(squeeze(segImg_surf(:,:,ceil(imgSize(3)./2))))
	axis tight equal
    title('Processed image')

	subplot(4,3,3)
	imagesc(squeeze(SurfaceSegMask(:,:,ceil(imgSize(3)./2))))
	axis tight equal
    title('Segmentation mask')



    subplot(4,3,4)
	imagesc(max(imgStack{Surface_SegChannel},[],3))
	axis tight equal
    ylabel('Surface, max. projection')

	subplot(4,3,5)
	imagesc(max(segImg_surf,[],3))
	axis tight equal
    
    subplot(4,3,6)
    imagesc(max(SurfaceSegMask,[],3))
	axis tight equal



    subplot(4,3,7)
	imagesc(squeeze(imgStack{Condensate_SegChannel}(:,:,ceil(imgSize(3)./2))))
	axis tight equal
    ylabel('Condensate, xy-section')
	
	subplot(4,3,8)
	imagesc(squeeze(segImg_cond(:,:,ceil(imgSize(3)./2))))
	axis tight equal

	subplot(4,3,9)
	imagesc(squeeze(CondensateSegMask(:,:,ceil(imgSize(3)./2))))
	axis tight equal



    subplot(4,3,10)
	imagesc(max(imgStack{Condensate_SegChannel},[],3))
	axis tight equal
    ylabel('Condensate, max. projection')
	
	subplot(4,3,11)
	imagesc(max(segImg_cond,[],3))
	axis tight equal

    subplot (4,3,12)
    imagesc(max(CondensateSegMask,[],3))
	axis tight equal
    
    colormap(inferno)


    % Uncomment the following two lines if you want to check the extracted
    % images one by one
   % fprintf('File name: %s\n',thisFilePath)
	waitforbuttonpress
	
	% --- Connected component segmentation of surfaces
	comps_surf = bwconncomp(SurfaceSegMask,18);
	numPxls = cellfun(@(elmt)numel(elmt),comps_surf.PixelIdxList);
	minPixels = Surface_minVol./(pixelSize.^2)./zStepSize;
	comps_surf.NumObjects = sum(numPxls>=minPixels);
	comps_surf.PixelIdxList = comps_surf.PixelIdxList(numPxls>=minPixels);
	numPxls = cellfun(@(elmt)numel(elmt),comps_surf.PixelIdxList);
		
	numSurf = comps_surf.NumObjects;
	numSurf_vec(ff) = numSurf;

	% --- Connected component segmentation of condensates
	comps_cond = bwconncomp(CondensateSegMask,18);
	numPxls = cellfun(@(elmt)numel(elmt),comps_cond.PixelIdxList);
	minPixels = Condensate_minVol./(pixelSize.^2)./zStepSize;
	comps_cond.NumObjects = sum(numPxls>=minPixels);
	comps_cond.PixelIdxList = comps_cond.PixelIdxList(numPxls>=minPixels);
	numPxls = cellfun(@(elmt)numel(elmt),comps_cond.PixelIdxList);
		
	numCond = comps_cond.NumObjects;
    numCond_vec(ff) = numCond;

    props_surf = regionprops3(comps_surf,imgStack{Surface_SegChannel},...
        'Volume','VoxelValues','Solidity','VoxelIdxList',...
        'BoundingBox','Centroid');

    Surf_Volume_array = [props_surf.Volume].*pixelSize.^2.*zStepSize;
    if numSurf>0
        Surf_Intensity_array = cellfun(@(vals)median(vals),props_surf.VoxelValues);
        Surf_Centroid_array = props_surf.Centroid...
            .*[pixelSize,pixelSize,zStepSize];
    else
        Surf_Intensity_array = [];
        Surf_Centroid_array = [];
    end
    Surf_Solidity_array = [props_surf.Solidity];
   

    props_drop = regionprops3(comps_cond,imgStack{Condensate_SegChannel},...
        'Volume','VoxelValues','Solidity','VoxelIdxList',...
        'BoundingBox','Centroid');

    Cond_Volume_array = [props_drop.Volume].*pixelSize.^2.*zStepSize;
    if numCond>0
        Cond_Intensity_array = cellfun(@(vals)median(vals),props_drop.VoxelValues);
        Cond_Centroid_array = props_drop.Centroid...
            .*[pixelSize,pixelSize,zStepSize];
    else
        Cond_Intensity_array = [];
        Cond_Centroid_array = [];
    end
    Cond_Solidity_array = [props_drop.Solidity];
   

    cond_intCell{ff} = cell(1,numQuantChannels);
    drop_intCell{ff} = cell(1,numQuantChannels);
    for qq = 1:numQuantChannels
        quantImg = imgStack{quantChannels(qq)};
        quantProps = regionprops3(comps_surf,quantImg,...
            'MeanIntensity','VoxelIdxList','VoxelValues');
        Surface_intCell{ff}{qq} = [quantProps.MeanIntensity];
        quantProps = regionprops3(comps_cond,quantImg,...
            'MeanIntensity','VoxelIdxList','VoxelValues');
        Condensate_intCell{ff}{qq} = [quantProps.MeanIntensity];
    end

    Condensate_volCell{ff} = Cond_Volume_array;
    Condensate_solCell{ff} = Cond_Solidity_array;
    %Condensate_eloCell{ff} = vertcat(Condensate_elongation{:});
    %Condensate_imgCell{ff} = vertcat(Condensate_centralSlices_store{:});
    Condensate_centCell{ff} = Cond_Centroid_array;
    Condensate_intCell{ff} = cell(1,numQuantChannels);
        
    Surface_volCell{ff} = Surf_Volume_array;
    Surface_solCell{ff} = Surf_Solidity_array;
    %Surface_eloCell{ff} = vertcat(Surface_elongation{:});
    %Surface_imgCell{ff} = vertcat(Surface_centralSlices_store{:});
    Surface_centCell{ff} = Surf_Centroid_array;
    Surface_intCell{ff} = cell(1,numQuantChannels);


    % ---
    % Beginning nearest-neighbor distance calculation

    if numCond==0
        Condensate_Surface_distCell{ff} = [];
        Surface_Condensate_distCell{ff} = NaN(size(Surf_Solidity_array));
    elseif numSurf==0
        Condensate_Surface_distCell{ff} = NaN(size(Cond_Solidity_array));
        Surface_Condensate_distCell{ff} = [];
    else
        pwDistMatrix = pdist2(Cond_Centroid_array,Surf_Centroid_array,'euclidean');
        Condensate_Surface_distCell{ff} = min(pwDistMatrix,[],2);
        Surface_Condensate_distCell{ff} = min(pwDistMatrix,[],1)';
    end

    % End nearest-neighbor distance calculation
    % ---
    
    %continue

    validFileFlag(ff) = true;

           
end

% % Retain only files that returned nuclei
% 
% condNames = condNames(validFileFlag);
% surf_intCell = surf_intCell(validFileFlag);
% drop_intCell = drop_intCell(validFileFlag);
% nuc_stdCell = nuc_stdCell(validFileFlag);
% nuc_medianVolCell = nuc_medianVolCell(validFileFlag);
% perNuc_countCell = perNuc_countCell(validFileFlag);
% 
% Droplet_xyVoxelSizeCell = Droplet_xyVoxelSizeCell(validFileFlag);
% Droplet_zVoxelSizeCell = Droplet_zVoxelSizeCell(validFileFlag);
% Surface_xyVoxelSizeCell = Surface_xyVoxelSizeCell(validFileFlag);
% Surface_zVoxelSizeCell = Surface_zVoxelSizeCell(validFileFlag);
% 
% Droplet_volCell = Droplet_volCell(validFileFlag);
% Droplet_solCell = Droplet_solCell(validFileFlag);
% Droplet_eloCell = Droplet_eloCell(validFileFlag);
% Droplet_imgCell = Droplet_imgCell(validFileFlag);
% Droplet_centCell = Droplet_centCell(validFileFlag);
% Droplet_Surface_distCell = Droplet_Surface_distCell(validFileFlag);
% Droplet_intCell = Droplet_intCell(validFileFlag);
% Droplet_nucIntCell = Droplet_nucIntCell(validFileFlag);
% 
% Surface_volCell = Surface_volCell(validFileFlag);
% Surface_solCell = Surface_solCell(validFileFlag);
% Surface_eloCell = Surface_eloCell(validFileFlag);
% Surface_imgCell = Surface_imgCell(validFileFlag);
% Surface_centCell = Surface_centCell(validFileFlag);
% Surface_Droplet_distCell = Surface_Droplet_distCell(validFileFlag);
% Surface_intCell = Surface_intCell(validFileFlag);
% Surface_nucIntCell = Surface_nucIntCell(validFileFlag);


%% --- Plotting of object properties

% Before we can plot the analysis results, we need to bring them into the
% format of vectors of numbers. Below, we carry out the according
% reorganization operations.

% ---
% Reorganize object shape properties into vector arrays. In this case, it
% is relatively easy. We have a cell array with one element per cell
% nucleus. We then call up the content of each element using the :
% operator. Each element contains a vector, which then gets concatenated by
% the vertcat command.

Condensate_Volume = vertcat(Condensate_volCell{:});
Condensate_Solidity = vertcat(Condensate_solCell{:});
% Distance

Surface_Volume = vertcat(Surface_volCell{:});
Surface_Solidity = vertcat(Surface_solCell{:});
%Distance

% ---
% Reorganize the object intensities into vector arrays. For the object
% intensities, it is a little tricker. We need to resort ot using the
% cellfun operation, as we need to make a two-index call, calling into each
% element of the top cell array for a specific sub-element.

Droplet_DNA_intensity = cellfun(@(elmt)elmt{1}',Condensate_intCell,...
    'UniformOutput',false);
Droplet_DNA_intensity = [Droplet_DNA_intensity{:}];
Droplet_Droplet_intensity = cellfun(@(elmt)elmt{2}',Condensate_intCell,...
    'UniformOutput',false);
Droplet_Droplet_intensity = [Droplet_Droplet_intensity{:}];
Droplet_Surface_intensity = cellfun(@(elmt)elmt{3}',Condensate_intCell,...
    'UniformOutput',false);
Droplet_Surface_intensity = [Droplet_Surface_intensity{:}];

Surface_DNA_intensity = cellfun(@(elmt)elmt{1}',Surface_intCell,...
    'UniformOutput',false);
Surface_DNA_intensity = [Surface_DNA_intensity{:}];
Surface_Droplet_intensity = cellfun(@(elmt)elmt{2}',Surface_intCell,...
    'UniformOutput',false);
Surface_Droplet_intensity = [Surface_Droplet_intensity{:}];
Surface_Surface_intensity = cellfun(@(elmt)elmt{3}',Surface_intCell,...
    'UniformOutput',false);
Surface_Surface_intensity = [Surface_Surface_intensity{:}];

% Histogram

figure(1)
clf

subplot(1,1,1)
cla
Surface_Droplet_distCell_vertical_array = vertcat(Surface_Condensate_distCell{:})
hist(Surface_Droplet_distCell_vertical_array,100);

[binCounts,binCenters] = ...
    hist(Surface_Droplet_distCell_vertical_array,nBins);
hold on
plot(binCenters,binCounts/(0.01*sum(binCounts)),...
    'k-','LineWidth',1)
hold off
ylabel('Count')
xlabel('Distance [µm]')
set(gca,'Box','on')
title('Selective analysis')


figure(2)

Surface_S2P_intensity = cellfun(@(elmt)elmt{1}',Surface_intCell,...
    'UniformOutput',false);
Surface_S2P_intensity = [Surface_S2P_intensity{:}];
subplot(1,1,1)
plot(Surface_Droplet_distCell_vertical_array,Surface_S2P_intensity,...
    'k.','MarkerSize',1)
title('Pol II clusters')
xlabel('Distance [µm]')
ylabel('S2P intensity')
set(gca,'Box','on')





subplot(2,3,1)
plot(Droplet_DNA_intensity,Droplet_Surface_intensity,...
    'k.','MarkerSize',1)
title('Pol II clusters')
xlabel('Distance [µm]')
ylabel('Surface intensity')
set(gca,'XLim',[0.4,2.4],'YLim',[0.4,2.4])

subplot(2,3,2)
plot(Droplet_DNA_intensity,Droplet_Droplet_intensity,...
    'k.','MarkerSize',1)
title('Pol II clusters')
xlabel('DNA intensity')
ylabel('Distance [µm]')
set(gca,'XLim',[0.4,2.4],'YLim',[0.4,2.4])

subplot(2,3,3)
plot(Droplet_Droplet_intensity,Droplet_Surface_intensity,...
    'k.','MarkerSize',1)
title('Pol II clusters')
xlabel('Pol II Ser5P intensity')
ylabel('Surface intensity')
set(gca,'XLim',[0.4,2.4],'YLim',[0.4,2.4])




subplot(2,3,4)
plot(Surface_DNA_intensity,Surface_Surface_intensity,...
    'k.','MarkerSize',1)
title('Surface')
xlabel('DNA intensity')
ylabel('Surface intensity')
set(gca,'XLim',[0.4,2.4],'YLim',[0.4,2.4])

subplot(2,3,5)
plot(Surface_DNA_intensity,Surface_Droplet_intensity,...
    'k.','MarkerSize',1)
title('Surface')
xlabel('DNA intensity')
ylabel('Pol II Ser5P intensity')
set(gca,'XLim',[0.4,2.4],'YLim',[0.4,2.4])

subplot(2,3,6)
plot(Surface_Droplet_intensity,Surface_Surface_intensity,...
    'k.','MarkerSize',1)
title('Surface')
xlabel('Pol II Ser5P intensity')
ylabel('Surface intensity')
set(gca,'XLim',[0.4,2.4],'YLim',[0.4,2.4])


%% --- Selective analysis of high Surface intensity objects

% Based on the scatter plots of object intensities, it looks as if there is
% a distinct population of high Surface intensity objects. This population is
% present no matter whether we segment Surface directly, or actually
% Pol II clusters. This suggests that Pol II clusters with high Surface
% levels might be reflecting the same population of objects as the high
% intensity Surface. We can check on this a bit more closely. In
% concrete terms, we will compare the shapes of the objects obtained by
% both approaches of segmentation.

SurfaceIntCutoff = 1.25;
Droplet_highSurface_inds = Droplet_Surface_intensity>SurfaceIntCutoff;
Surface_highSurface_inds = Surface_Surface_intensity>SurfaceIntCutoff;

figure(1)

subplot(3,3,7)
cla
nBins = 10;
[binCounts,binCenters] = ...
    hist(Droplet_DNA_intensity(Droplet_highSurface_inds),nBins);
plot([1,1],[0,0.3],'k-','Color',[0.5,0.5,0.5])
hold on
plot(binCenters,binCounts./sum(binCounts),...
    'k-','LineWidth',1)
[binCounts,binCenters] = ...
    hist(Surface_DNA_intensity(Surface_highSurface_inds),nBins);
plot(binCenters,binCounts./sum(binCounts),...
    'r--','LineWidth',1)
[binCounts,binCenters] = ...
    hist(Droplet_DNA_intensity(~Droplet_highSurface_inds),nBins);
plot(binCenters,binCounts./sum(binCounts),...
    'b:','LineWidth',1)
hold off
legend('Pol II Ser5P (high Surface)','Surface (high Surface)',...
    'Pol II Ser5P (low Surface)')
xlabel('Normalized DNA intensity')
ylabel('Normalized count')
set(gca,'Box','on')
title('Selective analysis')


subplot(3,3,8)
cla
[binCounts,binCenters] = ...
    hist(Droplet_Elongation(Droplet_highSurface_inds),nBins);
hold on
plot(binCenters,binCounts./sum(binCounts),...
    'k-','LineWidth',1)
[binCounts,binCenters] = ...
    hist(Surface_Elongation(Surface_highSurface_inds),nBins);
plot(binCenters,binCounts./sum(binCounts),...
    'r--','LineWidth',1)
[binCounts,binCenters] = ...
    hist(Droplet_Elongation(~Droplet_highSurface_inds),nBins);
plot(binCenters,binCounts./sum(binCounts),...
    'b:','LineWidth',1)
hold off
legend('Pol II Ser5P (high Surface)','Surface (high Surface)',...
    'Pol II Ser5P (low Surface)')
xlabel('Elongation')
ylabel('Normalized count')
set(gca,'Box','on')
title('Selective analysis')


subplot(3,3,9)
cla
[binCounts,binCenters] = ...
    hist(Condensate_Solidity(Droplet_highSurface_inds),nBins);
hold on
plot(binCenters,binCounts./sum(binCounts),...
    'k-','LineWidth',1)
[binCounts,binCenters] = ...
    hist(Surface_Solidity(Surface_highSurface_inds),nBins);
plot(binCenters,binCounts./sum(binCounts),...
    'r--','LineWidth',1)
[binCounts,binCenters] = ...
    hist(Condensate_Solidity(~Droplet_highSurface_inds),nBins);
plot(binCenters,binCounts./sum(binCounts),...
    'b:','LineWidth',1)
hold off
legend('Pol II Ser5P (high Surface)','Surface (high Surface)',...
    'Pol II Ser5P (low Surface)','Location','Southeast')
xlabel('Solidity')
ylabel('Normalized count')
set(gca,'Box','on')
title('Selective analysis')
