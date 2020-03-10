# EMsoft-utilities
Some useful code for interacting with EMsoft in MATLAB. These programs also allow you to use hexagonal grid data and index multiple phases, which cannot be done directly in EMsoft.

*Different versions of EMsoft have separate code here because of the differences in pattern orientation. Using the incorrect version can lead to incorrect indexing (and may not give an error). We have only tested EMsoft v4.0 and v4.3. Please let us know if one or neither of these versions seems to work for you.*

More information about EMsoft can be found at: https://github.com/EMsoft-org/EMsoft

Program descriptions:
- FindMapPoints_edaxang.m: Script that plots data from an EDAX .ang file and lets you click on points to extract info for that map point.
- GetPatternCenter_edaxang.m: Script that extracts pattern center from an EDAX .ang file and converts it to EMsoft coordinates.
- ComparePatterns.m: Compare 2 patterns by flashing between the two images.
- ConvertStorePatternsAsBinary.m: Stores all images in a folder into a .data file that can be read in by EMsoft.
- RunEMEBSD.m: Run EMsoft's EMEBSD program to simulate patterns and save as image files.
- RunEMEBSDDI.m: Run EMsoft's EMEBSDDI program to perform dictionary indexing and save .ang file.
- RunEMEBSDDI_refine.m: Run EMsoft's EMEBSDDI program to perform dictionary indexing, then immediately run EMsoft's EMFitOrientation program to refine dictionary indexing data and save .ang file.
- RunEMFitOrientation.m : Run EMsoft's EMFitOrientation program to refine dictionary indexing data and save .ang file.
- CombineAngFiles.m: Use to stitch together data from multiple .ang files (one phase each), selecting the phase with the highest NDP.
- TestErrorsTimings.m: Use this program to report errors and timings for EMsoft's EMEBSDDI program. Certain combinations of number of threads, numexptsingle/numdictsingle, and orientations causes errors (bug in EMsoft). Parameter values also greatly affect speed.

A walkthrough of how to perform a dictionary indexing run using EMsoft-utilities can be found in the file Instructions.pdf.

Some of these programs currently only work for EDAX .ang data. Please contact me if you would like the code adapted to data in another format. Feel free to email me at epang@mit.edu if you have any questions/difficulties/suggestions.
