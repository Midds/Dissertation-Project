# ImageStitching Instructions
To view a short video of this program running, a short <1min video is provided in this folder named "Artefact running example"

To run this program, the easiest way is to open Matlab and set the "MATLAB" folder as the current directory.
Before running the program, the VLFeat Matlab toolbox must be installed (for the use of the VLFeat implementation of the SIFT algorithm)
The VLFeat Toolbox has been included inside the "MATLAB" folder, so to add the toolbox to the MATLAB search path you just need to run the
following command in the MATLAB command window: "run vlfeat-0.9.20\toolbox\vl_setup".

NOTE: There seems to be a bug where running the above command for the first time will give an error where the file path is listed as incorrect. Should this occur
you can either try the full file path (which will likely give the same error), or re-download the VLFeat binary package from "http://www.vlfeat.org/download.html",
extract the files, and run the command using it's full file path. Doing the second method of downloading this package, extracting and then running the command with 
the full file path seems to fix the issue 100% of the time.

The "MATLAB" folder contains a number of .m files and folders needed to run the program.
The main program is contained in the "finalStitch.m" file, in addition to this "thirdIterationStitch" provides a file with adaptive frame selection.
The other files, "firstIterationStitch" and "secondIterationStitch" are from previous versions of the artefact as described in the project report.
In addition to these files there are a number of function files used within the aforementioned files.

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
attempt to keep increasing the distance between the frames being stitched, depending on the amount of matches that are found between each image pair. The threshold 
for this adaptive selection can be altered by changing the "matchThresh" variable on line 16.

Like the "finalStitch.m" file, you can control the number of images to stitch and the starting image in the same way.
IMPORTANT: For "thirdIterationStitch", If changing the filename on line 24 - from which the program will get the images from. The filename on line 118 must be changed to match this.

-------------------------- KNOWN BUGS --------------------------

1. There is a bug in the file where certain images will not stitch and the program will crash. This usually happens when there are not enough matches found in the images. The barret2 
images are very hard to find suitable places to get matches from. Barret1 is more suitable although setting the "startImage" variable to certain values
will result in this bug occurring. Changing the "euclideanThresh" variable can help this issue. So can using the "thirdIterationStitch" file instead of "finalStitch".

2. There is a bug where the final mosaic will only sometimes be very distorted. This bug has no discernable cause, and running the exact same file again will often resolve the issue. 
However using the "clear;" command seems to help the issue although not always. If this happens consider clearing the workspace and running again. 

3. In some cases using too many images in one stitch will end up preventing the final stitch from working. An error will be thrown saying certain indices or matrix dimensions do not match.

4. Sometimes trying to stitch Barret images using "finalStitch" will result in the program running into an infinite loop around line 142. This is because the program cannot find enough
suitable matches. If this happens change the images to stitch, or change the "euclideanThresh" variable to a higher value (doesn't alwasy work and can lead to a worse/distorted final stitch).