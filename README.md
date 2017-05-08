# ImageStitching Instructions
To run this program, the easiest way is to open Matlab and set the "MATLAB" folder as the current directory.
Before running the program, the VLFeat Matlab toolbox must be installed (for the use of the VLFeat implementation of the SIFT algorithm)
The VLFeat Toolbox has been included inside the "MATLAB" folder, so to add the toolbox to the MATLAB search path you just need to run the
following command in the MATLAB command window: "run vlfeat-0.9.20\toolbox\vl_setup".

Note: There is sometimes a bug where running the above command will give an error where the file path is listed as incorrect. Should this occur
you will need to either try the full file path, or re-download the VLFeat binary package from "http://www.vlfeat.org/download.html",
extract the files, and run the command using it's full file path.

The "MATLAB" folder contains a number of .m files and folders needed to run the program.
The main program is contained in the "finalStitch.m" file.
The other files, "firstIterationStitch", "secondIterationStitch", and "thirdIterationStitch" are from previous versions of the artefact 
as described in the project report.



% pipeline as follows
% - read in images to mosaic
% - find sift for each and check there's enough matches, if not read in
% another image
% - use ransac on sift matches to get homograpy
% - refine homography using levenberg-marquardt
% - stitch images