clear all
	
%% Plot the example images

SourceFileCell = {...
	'./ExtractedStacks/Cond_1/Image_10.mat',...
	'./ExtractedStacks/Cond_2/Image_3.mat',...
    './ExtractedStacks/Cond_3/Image_4.mat',...
	'./ExtractedStacks/Cond_4/Image_2.mat',...
	'./ExtractedStacks/Cond_5/Image_5.mat',...
	'./ExtractedStacks/Cond_6/Image_6.mat',...
    };

Surf_blurRange = 0.1; % micrometers
Cond_blurRange = 0.1; % micrometers

imgSources = [1,2,3,4,5,6];
imgRanges = {...
	[37,37,72,72]+[-6,+6,-6,+6]};
imgRanges = {...
	[75,75,110,110]+[-20,+20,-15,+15],...
    [75,75,110,110]+[-20,+20,-15,+15],...
    [75,75,115,115]+[-20,+20,-15,+15],...
    [60,60,142,142]+[-20,+20,-15,+15],...
    [100,100,140,140]+[-20,+20,-15,+15],...
    [75,75,115,115]+[-20,+20,-15,+15],...
    };
zCoordinate = [35,30,26,5,35,10];
rotate_flag = [false,false,false,false,false,false];

plotTitles = {'+Surf','+Surf','+Surf','+Surf','+Surf','-Surf'};

scaleBar = 10.0; % in microns

% --- parameters end

numPlots = numel(imgRanges);

figure(1)
clf

for pp = 1:numPlots
	
	pp
	
	thisFilePath = SourceFileCell{pp};
	loadStruct = load(thisFilePath,...
		'imgStack','imgSize','pixelSize','zStepSize','condName');
	imgStack = loadStruct.imgStack;
	imgSize = loadStruct.imgSize;
	pixelSize = loadStruct.pixelSize;
	zStepSize = loadStruct.zStepSize;
    condName = loadStruct.condName;

	thisImg = imgStack;
	thisPixelSize = pixelSize;
	
	% --- Removal of DNA background
	Surf_img = double(thisImg{1}(:,:,zCoordinate(pp)));
	Cond_img = (double(thisImg{2}(:,:,zCoordinate(pp))));

	if rotate_flag(pp)
		Surf_img = Surf_img';
		Cond_img = Cond_img';
	end

	Surf_img = imgaussfilt(Surf_img,...
		Surf_blurRange./thisPixelSize);
	Cond_img = imgaussfilt(Cond_img,...
		Cond_blurRange./thisPixelSize);
	
	thisImgRange = round(imgRanges{pp}./thisPixelSize)+1;
	Surf_img = Surf_img(...
		thisImgRange(1):thisImgRange(2),...
		thisImgRange(3):thisImgRange(4));
	Cond_img = Cond_img(...
		thisImgRange(1):thisImgRange(2),...
		thisImgRange(3):thisImgRange(4));
	thisSize = size(Cond_img);
	
    Surf_lims = prctile(Surf_img(:),[0.01,99.99]);
    Cond_lims = prctile(Cond_img(:),[0.01,99.99]);

	subplot(3,numPlots,numPlots.*0+pp)
	imagesc([0,thisSize(2)].*thisPixelSize,...
		[0,thisSize(1)].*thisPixelSize,...
		Surf_img,Surf_lims)
	axis equal tight
	set(gca,'XTick',[],'YTick',[])
	set(gca,'YDir','normal')
	hold on
	plot([1,1+scaleBar],thisSize(1).*thisPixelSize-[1,1],...
		'w-','LineWidth',3) % Plot scale bar
	textObject = ...
		text(1+0.5.*scaleBar,thisSize(1).*thisPixelSize-2.0,...
		sprintf('%d \\mum',scaleBar));
	set(textObject,'Color',[1,1,1],'HorizontalAlignment','center',...
		'FontSize',12)
	ylabel('Surface')
    title(condName)
    xlabel(sprintf('Int. range [%d,%d]',...
        round(Surf_lims(1)),round(Surf_lims(2))))
	
	
	subplot(3,numPlots,numPlots.*1+pp)
	imagesc([0,thisSize(2)].*thisPixelSize,...
		[0,thisSize(1)].*thisPixelSize,...
		Cond_img,Cond_lims)
	axis equal tight
	set(gca,'XTick',[],'YTick',[])
	set(gca,'YDir','normal')
	ylabel('X-Motif')
    xlabel(sprintf('Int. range [%d,%d]',...
        round(Cond_lims(1)),round(Cond_lims(2))))
	
	% --- Color merge

    subplot(3,numPlots,numPlots.*2+pp)

	% Prepare the color channels for RGB plotting
	this_magenta_plot = ...
		(Cond_img-Cond_lims(1))./diff(Cond_lims);
	this_green_plot = ...
		(Surf_img-Surf_lims(1))./diff(Surf_lims);
	thisSize = size(Surf_img);
	
	redChannel = this_magenta_plot;
	blueChannel = this_magenta_plot;
	greenChannel = this_green_plot;
	rgb_img = zeros(thisSize(1),thisSize(2),3);
	rgb_img(:,:,1) = redChannel;
	rgb_img(:,:,2) = greenChannel;
	rgb_img(:,:,3) = blueChannel;
	
	image([0,thisSize(2)].*thisPixelSize,...
		[0,thisSize(1)].*thisPixelSize,...
		rgb_img)
	set(gca,'XTick',[],'YTick',[])
	set(gca,'YDir','normal')
	ylabel('Surface / X-Motif')
	axis equal tight
	
end