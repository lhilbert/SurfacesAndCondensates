clear all

% This is a new type of variable, it can take the values true or false.
% This logical variable, also called a Boolean, an be used as a kind of
% switch, to turn on / off certain functions to be carried out in a script.
% In this case, when you set it to "true", a central z section of the image
% data will be shown for each of the extracted image positions. This can be
% very helpful to check what is going on, that the color channels are
% chosen properly, and just to get a feeling for what we have in front of
% us. On the other hand, it requires that you interact with the images
% continuously during the extraction process, and you might not want to
% babysit your code - in that case, just change the value of this variable
% to false, and no data will be displaye during extraction, and you do not
% have to click on anything.
plotFlag = false;

% Here, you can specify in which directories on your hard drive the raw
% image data can be found. You can specify not only one directory, but
% several. In case you have several experimental conditions that you would
% like to compare, you can already sort your data on your hard drive: each
% condition should, ideally, be sorted into one directory. For
% demonstration purpose, I have only included a single directory in this
% first illustration of the extraction script.
%
% There is also a few details to note about the directory naming. the first
% thing, depending on which operating system (Mac vs Windows mostly) you
% are using, the file separator symbol might differ. On Mac, the forward
% slash / is used to separate directories. On Windows, the backslash \ is
% used, and you might have to fix this in the script if you are working on
% Windows.
%
% The second point to note is the use of the double point ..
% the double point means "go up one directory". So, if you are in the
% directory for Day 1, but the raw data are stored one directory level
% higher, you can use .. to go to the higher directory level, and then
% choose the directory that contains the raw data after the .. symbol. We
% use this here, and the best is to simply look at how it is used:
sourceDirectories = { ...
	'C:\Users\X\Desktop\UNI\HIWI\Analysis_titration_testdata\Control\0_0\',...	
    'C:\Users\X\Desktop\UNI\HIWI\Analysis_titration_testdata\Control\1_5\',...	
    'C:\Users\X\Desktop\UNI\HIWI\Analysis_titration_testdata\Control\5_0\',...
    'C:\Users\X\Desktop\UNI\HIWI\Analysis_titration_testdata\Surface\0_0\',...
    'C:\Users\X\Desktop\UNI\HIWI\Analysis_titration_testdata\Surface\1_5\',...
    'C:\Users\X\Desktop\UNI\HIWI\Analysis_titration_testdata\Surface\5_0\',...
            };

% The extracted files need to be placed on the hard drive, and you have to
% tell the script exactly where that should be. If you do not change this,
% a directory called 'ExtractedStacks' will be generated in the folder that
% contains this analysis script. That is not a bad choice, so unless you
% have a specific reason, you do not really need to change this directory.
extractTargetFolder = 'ExtractedStacks_Titration';

% While a lot is automated in this script, you still need to tell it which
% directory corresponds to which experimental condition. You can also
% assign several directories to the same condition, so you can pool data
% that come from different directories.
condInds = (1:6)';
% Each condition needs a name, and you need to enter it in this cell array.
% You can be very short, it just has to work for you to distinguish the
% different conditions lateron.
condLabels = {...
	'C0_0','C1_5','C5_0','S0_0','S1_5','S5_0'
	};
% This variable is an interesting setting if you want to first test the
% extraction script quickly. If you enter a number, for example 2, you can
% limit how many different positions are extracted from a given file on the
% hard drive. So, for example, if you have raw data with 10 or 20 positions
% per each file, that will take a long time and fill a lot of space on your
% hard drive. Instead, you can set a lower number here, and then first get
% a quicker look into the data set. If you want to extract all positions
% that are available, enter the value 'Inf', which stands for infinity.
maxNumSeries = Inf;

% The skip list tells the script to not extract data from specific
% directories. You can enter the number for the directory, so if you put
% the number 1 into this vector, the first directory in the directory list
% will not be extracted. This can be extremely useful if only additional
% data in one directory comes, or parts of a certain condition need to be
% removed because technical problems were discovered. The other directories
% do not need to be re-extracted, so this can save a lot of time.
% Leave empty array [] if all directories shoudl be processed.
skipList = []; % Directories to skip, for example if already done

% You also should choose which color channels you actually want to extract
% from the raw data. Often times, more color channels are recorded than
% should ultimately be analyzed. For example, our data have four color
% channels, but we only want to work with the first three. Accordingly, we
% set the color channel numbers 1, 2, 3. You can even use this option to
% reorder channels, so you could for example use the number 3,2,1 to
% reverse the channel order that is stored in the extracted data.
useChannel_inds = [1,2]; %Surface, Nanomotif
numChannels = numel(useChannel_inds);

% Here, you can specify the color ranges used to display the different
% colro channels. If you put [-Inf,Inf], this will automatically adjust the
% color range to the minimum and the maximum values found in a given image.
% Unless you have a strong reason to change this, just leave it as it is.
scaleChannels = {[-Inf,Inf],[-Inf,Inf],[-Inf,Inf]};

% --- From here on, the script does not need to be changed.
% We will go over a few main points of what is contained in the script.

% The script runs through one main for loop, which takes it through all the
% directories that you have specified. To that end, we first get the number
% of directories, and the loop over all these directories.

numDirs = numel(sourceDirectories);

for cc = 1:numDirs
    
    % The fprintf command can be used to generate character strings that
    % also use content from variables. In this case, we use it to report on
    % the progress of the script by printing to the console.
	fprintf('Extracting image data from directory %d of %d\n',cc,numDirs)
	
    % Here, we are using a case distinction with an if-else statement. The
    % idea is that we check whether the directory is on the skip list, and
    % only if it is not on the skip list, we carry out the extraction
    % procedure.
	if ismember(cc,skipList)
		
		fprintf('On skip list, skipping to next directory.\n')
		
	else
		
        % In the following two lines, we get the name of the directory from
        % the cell array, and use the function rdir (recursive directory
        % search for files) to find all the files in the directory and its
        % subdirectories that end on the .nd2 file ending.
        thisDir = sourceDirectories{cc};
		listing = rdir([thisDir,'*.nd2'],'~contains(name,''._'')');
		
        % For the current directory, we now will go and work through all
        % the files one after the other. So, again, we need to find out how
        % many files exist, and then loop over those files.
		numFiles = numel(listing);
				
		stackCounter = 0;
		mkdir(sprintf('./%s/Cond_%d',extractTargetFolder,cc))
		
        % This is the for loop that runs over the different files.
		for ff = 1:numFiles
			
			fprintf('Extracting images from file %d of %d\n',ff,numFiles)
			
            % Here, we make the filepath using the list of discovered files
            % from the recursive directory search.
			combined_filepath = ...
				fullfile(listing(ff).name);
			% This file path is now used to call the function nd2read,
            % which is basically a function containing the read-in routines
            % we developed in the last script. You can into that function
            % and take a look...
			[imgCell,imgSizeCell,...
				this_voxelSize_cell,seriesName_cell] = ...
				nd2read(combined_filepath,maxNumSeries);
			
            % The information from the file read process is stored in the
            % format of cell arrays. For the pixel xy and z size, this is
            % not the most convenient, and we would like to convert it into
            % normal vector type arrays. For this we use cellfun, but the
            % usage is a little heavy for the first day of the course, so
            % we will skip this for now.
			pixelSizeVec = cellfun(@(elmt)elmt(1),this_voxelSize_cell);
			zzStepVec = cellfun(@(elmt)elmt(3),this_voxelSize_cell);
			
            % The data read in process also finds the image stacks
            % corresponding to the multiple positions in the sample.
            % Accordingly, there will be several images in the imgCell
            % array. To sort the color channels for all of these images, we
            % again use a for loop that runs over these images, and uses
            % the requested channel order to reorder the color channels.
			numImgs = numel(imgCell);
			for kk = 1:numImgs
				imgCell{kk} = imgCell{kk}(useChannel_inds);
			end
			
            % The plotting here will be excuted in case the plot flag
            % varaiable was set to true. We will not go into the plotting
            % details here, because it will be part of your independent
            % coding exercise later.
			if plotFlag
				
				for kk = 1:numImgs
					
					thisPixelSize = pixelSizeVec(kk);
					thisSize = imgSizeCell{kk};
					
					numChannels = numel(imgCell{kk});
					
					for cc_2 = 1:numChannels
						
						subplot(1,numChannels,cc_2)
						
						imagesc([0,thisSize(2)].*thisPixelSize,...
							[0,thisSize(1)].*thisPixelSize,...
							imgCell{kk}{cc_2}(:,:,ceil(thisSize(3)./2)),...
							scaleChannels{1})
						axis equal tight
						
						xlabel('x [\mum]')
						
						colormap((gray))
						
					end
					
					disp(sprintf( ...
						'Directory number %d, Condition %s, File %d/%d, Image %d/%d.', ...
						cc,condLabels{condInds(cc)},ff,numFiles,kk,numImgs))
					
					waitforbuttonpress
					
				end
				
			end
			
			% --- saving of the retained images

			condName = condLabels{cc};
			condInd = condInds(cc);
			
			% Save all images into a separate file, but contained in a
            % folder that contains all images for this condition
			
			numImgs = numel(imgCell);
			for kk = 1:numImgs
                % We use a stack counter, which is increased by one for
                % every single image that is saved. This means, if for a
                % given directory several files are processed, for the
                % second file we just keep counting up in the stack
                % counter, it is not reset to zero for a new file. If that
                % sounds abstract and confusing, just look at the output in
                % the ExtractedStacks directory, you will see that all the
                % files have an increasing number at the end. The stack
                % counter is used to achieve this increasing number across
                % several raw image files from the same directory.
                stackCounter = stackCounter+1;

                % Here, we are simply pulling the information for this
                % particular image from the overall cell arrays.
                imgStack = imgCell{kk};
				imgSize = imgSizeCell{kk};
				pixelSize = pixelSizeVec(kk);
				zStepSize = zzStepVec(kk);
				
                % Now comes the actual saving process. We need to specify
                % the folder and file name, to which we would like to save.
                % That is achieved using the sprintf command, which can
                % create a character string based on something called
                % "regular exressions". Again, we will not go into details
                % today, but you can look up sprintf and regular
                % expressions if you are very interested.
                targetFile = sprintf('./%s/Cond_%d/Image_%d.mat',...
					extractTargetFolder,cc,stackCounter);
                % The target file is only one part of the save command. We
                % further need to specify which variables should be saved
                % into the file. This is relatively simple, you just list
                % the variables which should be saved to the file, and
                % MatLab will take care of the saving process.
                save(targetFile,...
					'imgStack','imgSize','pixelSize','zStepSize',...
					'condInd','condName')
			end
			
		end
		
	end
	
end