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



# Example data
Example image data for processing are available via Zenodo. We recommend downloading 3-5 of the .nd2 files from the following data repository

https://zenodo.org/uploads/13785805
