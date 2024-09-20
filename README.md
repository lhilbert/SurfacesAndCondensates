# SurfacesAndCondensates
 MatLab Image Analysis Scripts for Condensates that form on condensation surfaces

# Main steps of the analysis

## Data organization
Place all data in a folder, and make subfolders for each experimental condition, into which you sort the image data.

## Extraction and condition labeling of image data before processing
Enter the subdirectories into the file MultiPosition_Extraction_nd2.m (If you use other formats than .nd2, you will need a different importer function). This goes into the command assigning the varialbe sourceDirectories.
Enter into the variable extractTargetFolder into which target folder the extracted images should be saved.
Enter the number of conditions in the variable condInds, in the format condInds = {1:6}, if for example you want to have six different conditions.
Enter the short names of the conditions into the variable condLabels.
In the variable useChannel_inds, enter which channel contains the nanomotifs, and which the DNA surface strands.

## Running actual image analysis over the extracted data

In the file ObjectAnalysis_Titration.m set the directory to which you extracted the image data.
Also set the analysis parameters in the top section.
You can run the script in review and parameter adjustment mode if you comment out the following lines in the first for loop:
fprintf('File name: %s\n',thisFilePath)
waitforbuttonpress
Once the analysis is finished, the results of the analysis will be stored into a file sortedREsultBundle.mat

## Graphs and statistics of the analysis results

After the image analysis itself has completed, you can use further scripts to generate overview plots and statistical analysis.

List of such files:
- OverviewPlots.m
- Titration_ExampleImages.m
- ...

# Example data
Example image data for processing are available via Zenodo. We recommend downloading 3-5 of the .nd2 files from the following data repository

First experimental repeat:
[https://zenodo.org/uploads/13785805](https://zenodo.org/records/13788484)

Second experimental repeat:
[https://doi.org/10.5281/zenodo.13788846](https://doi.org/10.5281/zenodo.13788846)
