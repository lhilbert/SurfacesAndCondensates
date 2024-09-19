clear all

% Specify the directory that contains the extracted files from the last
% step, where you extracted from the raw files obtained from the microscope
sourceDirectory = ...
    './ExtractedStacks/**/';

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
Condensate_segBlurSigma_BG_removal = 5.0; % in microns
% Number of standard deviations above background for the robust threshold
Condensate_seg_numStdDev = 6.0;

% Blurring and segmentation parameters for the Surface
Surface_segBlurSigma_object = 0.5; % in microns
Surface_segBlurSigma_BG_removal = 5.0; % in microns
Surface_seg_numStdDev = 6.0; % number of standard deviations in robust threshold

% Minimum volumes for condensate and surface objects
Condensate_minVol = 0.02; % cubic microns
Surface_minVol = 0.02; % cubic microns

surf_DBSCAN_epsilon = 3;
cond_DBSCAN_epsilon = 3;

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

parfor ff = 1:numFiles
	
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
%     fprintf('File name: %s, Condition: %s\n',...
%         thisFilePath,thisCondName)
% 	waitforbuttonpress
	
	% --- Connected component segmentation of surfaces
	comps_surf = bwconncomp(SurfaceSegMask,18);
	numPxls = cellfun(@(elmt)numel(elmt),comps_surf.PixelIdxList);
	minPixels = Surface_minVol./(pixelSize.^2)./zStepSize;
	comps_surf.NumObjects = sum(numPxls>=minPixels);
	comps_surf.PixelIdxList = comps_surf.PixelIdxList(numPxls>=minPixels);
	numPxls = cellfun(@(elmt)numel(elmt),comps_surf.PixelIdxList);

	% --- Connected component segmentation of condensates
	comps_cond = bwconncomp(CondensateSegMask,18);
	numPxls = cellfun(@(elmt)numel(elmt),comps_cond.PixelIdxList);
	minPixels = Condensate_minVol./(pixelSize.^2)./zStepSize;
	comps_cond.NumObjects = sum(numPxls>=minPixels);
	comps_cond.PixelIdxList = comps_cond.PixelIdxList(numPxls>=minPixels);
	numPxls = cellfun(@(elmt)numel(elmt),comps_cond.PixelIdxList);
		
    % --- connect nearby objects into clusters using DBSCAN, end

    if comps_surf.NumObjects>0 && surf_DBSCAN_epsilon > 0

        props_surf = regionprops3(comps_surf,imgStack{Surface_SegChannel},...
            'Centroid');

        centroid_coords = ...
            props_surf.Centroid.*[pixelSize,pixelSize,zStepSize];
        dbscan_inds = ...
            dbscan(centroid_coords,surf_DBSCAN_epsilon,1);

        unique_inds = unique(dbscan_inds);
        num_inds = numel(unique_inds);
        updated_comps = comps_surf;
        updated_comps.NumObjects = num_inds;
        updated_comps.PixelIdxList = cell(1,num_inds);
        for ii = 1:num_inds
            updated_comps.PixelIdxList{ii} = ...
                sort(vertcat(comps_surf.PixelIdxList{...
                dbscan_inds==unique_inds(ii)} ...
                ));
        end
        comps_surf = updated_comps;

        %Display hierarchical clustering
        LL = labelmatrix(comps_surf);

        centerPlaneInd = round(imgSize(3).*0.5);

        subplot(1,3,1)
        imagesc(squeeze(max(SurfaceSegMask,[],3)))
        axis tight equal
        title('Max cluster mask',FontSize=20);

        subplot(1,3,2)
        imagesc(squeeze(max(LL,[],3)))
        axis tight equal
        title('Max cluster comps',FontSize=20);

        subplot(1,3,3)
        imagesc(squeeze(LL(:,:,centerPlaneInd)))
        axis tight equal
        title('Mid cluster comps',FontSize=20);
        %set(gca,'Colormap',lines)

        %waitforbuttonpress

    end

    if comps_cond.NumObjects>0 && cond_DBSCAN_epsilon > 0

        props_cond = regionprops3(comps_cond,imgStack{Condensate_SegChannel},...
            'Centroid');

        centroid_coords = ...
            props_cond.Centroid.*[pixelSize,pixelSize,zStepSize];
        dbscan_inds = ...
            dbscan(centroid_coords,cond_DBSCAN_epsilon,1);

        unique_inds = unique(dbscan_inds);
        num_inds = numel(unique_inds);
        updated_comps = comps_surf;
        updated_comps.NumObjects = num_inds;
        updated_comps.PixelIdxList = cell(1,num_inds);
        for ii = 1:num_inds
            updated_comps.PixelIdxList{ii} = ...
                sort(vertcat(comps_cond.PixelIdxList{...
                dbscan_inds==unique_inds(ii)} ...
                ));
        end
        comps_cond = updated_comps;

        %Display hierarchical clustering
        LL = labelmatrix(comps_cond);

        centerPlaneInd = round(imgSize(3).*0.5);

        subplot(1,3,1)
        imagesc(squeeze(max(CondensateSegMask,[],3)))
        axis tight equal
        title('Max cluster mask',FontSize=20);

        subplot(1,3,2)
        imagesc(squeeze(max(LL,[],3)))
        axis tight equal
        title('Max cluster comps',FontSize=20);

        subplot(1,3,3)
        imagesc(squeeze(LL(:,:,centerPlaneInd)))
        axis tight equal
        title('Mid cluster comps',FontSize=20);
        %set(gca,'Colormap',lines)

        %waitforbuttonpress

    end


    % --- connect nearby objects into clusters using DBSCAN, end

    props_surf = regionprops3(comps_surf,imgStack{Surface_SegChannel},...
        'Volume','VoxelValues','Solidity','VoxelIdxList',...
        'BoundingBox','Centroid');

    numSurf = comps_surf.NumObjects;
    numSurf_vec(ff) = numSurf;

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

    numCond = comps_cond.NumObjects;
    numCond_vec(ff) = numCond;

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

    for qq = 1:numQuantChannels
        quantImg = imgStack{quantChannels(qq)};
        quantProps = regionprops3(comps_surf,quantImg,...
            'MeanIntensity','VoxelIdxList','VoxelValues');
        Surface_intCell{ff}{qq} = [quantProps.MeanIntensity];
        quantProps = regionprops3(comps_cond,quantImg,...
            'MeanIntensity','VoxelIdxList','VoxelValues');
        Condensate_intCell{ff}{qq} = [quantProps.MeanIntensity];
    end


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

%% Sort into conditions

uniqueCondNames = unique(condNames);
numPlotSets = numel(uniqueCondNames);
fileIndsCell = cell(1,numPlotSets);
numFiles_perCond = zeros(1,numPlotSets);
for cc = 1:numPlotSets
	fileIndsCell{cc} = cellfun(...
		@(elmt)strcmp(elmt,uniqueCondNames{cc}),condNames);
	numFiles_perCond(cc) = sum(fileIndsCell{cc});
end

sortedCondNames = cell(1,numPlotSets);
sortedNumFiles = zeros(1,numPlotSets);

sortedDropletPixelSize_xy = cell(1,numPlotSets);
sortedDropletPixelSize_z = cell(1,numPlotSets);
sortedSurfacePixelSize_xy = cell(1,numPlotSets);
sortedSurfacePixelSize_z = cell(1,numPlotSets);

sortedDropSurfDistCell = cell(1,numPlotSets);
sortedSurfDropDistCell = cell(1,numPlotSets);

sortedDropletNumCell = cell(1,numPlotSets);
sortedDropletVolCell = cell(1,numPlotSets);
sortedDropletSolCell = cell(1,numPlotSets);
sortedDropletEloCell = cell(1,numPlotSets);
sortedDropletCentralSliceCell = cell(1,numPlotSets);
sortedDropletCentroidsCell = cell(1,numPlotSets);
sortedDropletIntCell = cell(1,numQuantChannels);

sortedSurfaceNumCell = cell(1,numPlotSets);
sortedSurfaceVolCell = cell(1,numPlotSets);
sortedSurfaceSolCell = cell(1,numPlotSets);
sortedSurfaceEloCell = cell(1,numPlotSets);
sortedSurfaceCentralSliceCell = cell(1,numPlotSets);
sortedSurfaceCentroidsCell = cell(1,numPlotSets);
sortedSurfaceIntCell = cell(1,numQuantChannels);

for qq = 1:numQuantChannels
	sortedDropletIntCell{qq} = cell(1,numPlotSets);
	sortedSurfaceIntCell{qq} = cell(1,numPlotSets);
end

for cc = 1:numPlotSets
	
	sortedCondNames{cc} = ...
		condNames(fileIndsCell{cc});
	sortedCondNames{cc} = sortedCondNames{cc}{1};
	
	sortedNumFiles(cc) = sum(fileIndsCell{cc});

    sortedDropletNumCell{cc} = numCond_vec(fileIndsCell{cc});
    sortedSurfaceNumCell{cc} = numSurf_vec(fileIndsCell{cc});

    CS_dists = vertcat(Condensate_Surface_distCell{fileIndsCell{cc}});
    SC_dists = vertcat(Surface_Condensate_distCell{fileIndsCell{cc}});

    sortedDropSurfDistCell{cc} = CS_dists;
    sortedSurfDropDistCell{cc} = SC_dists;

	Drop_vols = vertcat(Condensate_volCell{fileIndsCell{cc}});
    Drop_sols = vertcat(Condensate_solCell{fileIndsCell{cc}});
	Drop_ints = Condensate_intCell(fileIndsCell{cc});
    
    sortedDropletVolCell{cc} = Drop_vols;
    sortedDropletSolCell{cc} = Drop_sols;
                    
	Surf_vols = vertcat(Surface_volCell{fileIndsCell{cc}});
    Surf_sols = vertcat(Surface_solCell{fileIndsCell{cc}});
	Surf_ints = Surface_intCell(fileIndsCell{cc});
    
    sortedSurfaceVolCell{cc} = Surf_vols;
    sortedSurfaceSolCell{cc} = Surf_sols;
	
	for qq = 1:numQuantChannels
        
		sortedDropletIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)Condensate_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedDropletIntCell{qq}{cc} = vertcat(...
            sortedDropletIntCell{qq}{cc}{:});

		sortedSurfaceIntCell{qq}{cc} = vertcat(arrayfun(...
			@(ind)Surface_intCell{ind}{qq},....
			find(fileIndsCell{cc}),...
			'UniformOutput',false));
		sortedSurfaceIntCell{qq}{cc} = vertcat(...
            sortedSurfaceIntCell{qq}{cc}{:});
        
	end
	
end

save('sortedResultsBundle')