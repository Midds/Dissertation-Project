# ImageStitching Instructions
To run this program, the easiest way is to open Matlab and set the "MATLAB" folder as the current directory.
Before running the program, the VLFeat Matlab toolbox must be installed (for the use of the VLFeat implementation of the SIFT algorithm)
The VLFeat Toolbox has been included inside the "MATLAB" folder, so to add the toolbox to the MATLAB search path you just need to run the
following command in the MATLAB command window: "run vlfeat-0.9.20\toolbox\vl_setup".

NOTE: There is sometimes a bug where running the above command will give an error where the file path is listed as incorrect. Should this occur
you will need to either try the full file path, or re-download the VLFeat binary package from "http://www.vlfeat.org/download.html",
extract the files, and run the command using it's full file path.

The "MATLAB" folder contains a number of .m files and folders needed to run the program.
The main program is contained in the "finalStitch.m" file.
The other files, "firstIterationStitch", "secondIterationStitch", and "thirdIterationStitch" are from previous versions of the artefact 
as described in the project report.

In "finalStitch.m", there are a number of variables at the beginning of the file that can be changed to alter what image is being stitched,
as well as how many images will be stitched. These should be self-explained in the comments in the file. However, the most important variables 
to note are:
		"numToStitch" - changes the number of images in a stitch
		"startImage" - corresponds to the image number in the chosen image folder
		"euclideanThresh" - a threshold variable used in the matching of SIFT features (a higher number is needed for images with lower matches)
To change the source of the images being stitched, the "filename" variable should be altered in the Reading images loop. For example choose:
	filename = sprintf('barret1/im%d.jpeg', n); % This will use images from the barret1 video 
	filename = sprintf('london/im%d.jpeg', n); % This will use images from the london poster video
	filename = sprintf('tiger/im%d.jpeg', n); % This will use images from the tiger poster video
In addition to changing the filename, line 32 allows the resizing of images. I recommend up sizing barret images as they are fairly small. Down sizing 
the "tiger" or "london" images to 0.5 is also highly recommended to greatly reduce the time required to make stitches, as these images are very large.
	


The "thirdIterationStitch" file is set-up very similarly to the "finalStitch" file. The major differences are the use of a different method of 
ransac (no levenberg-marquardt refinement). In addition to this, the "thirdIterationStitch" file contains a form of adaptive frame selection where the program will
attempt to keep increasing the distance between the frames being stitched, depending on the amount of matches that are found between each image pair.

-------------------- KNOWN BUGS --------------------

1. There is a bug in the file where certain images will not stitch. This usually happens when there are not enough matches found in the images. The barret2 
images are very hard to find suitable places to get matches from. Barret1 is more suitable although setting the "startImage" variable to certain values
will result in this bug occurring. Changing the "euclideanThresh" variable can help this issue. So can using the "thirdIterationStitch" file instead of "finalStitch".



	
% pipeline as follows
% - read in images to mosaic
% - find sift for each and check there's enough matches, if not read in
% another image
% - use ransac on sift matches to get homograpy
% - refine homography using levenberg-marquardt
% - stitch images