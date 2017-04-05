close all;
%% Loading input video file
%pV = VideoReader('9poster4.mp4');
pV = VideoReader('tiger2.mp4');
%pV = VideoReader('Barretxx1.mpg');

%% PRE-PROCESSING
% Splitting video into a number of keyframes (to emulate lower framerate).
% This takes some compuational time, but is only necessary for this
% emulated endoscopy, so the time taken is irrelevant to final product.

% framesToGet will give the final number of frames used.
% keyFrameNo is used to decide framesToGet and is used in the for loop.
% A higher keyFrameNo means less frames in the final video.
totalFrames = pV.NumberOfFrames;
keyFrameNo = 1; % gets every x frame, change number to change x
framesToGet = round(totalFrames / keyFrameNo);

%% Looping through the video picking out frames at the specified interval
j = 1; % Count variable for loop
frameArray = {0:framesToGet}; % Create an array to hold the selected frames

for i = 1:framesToGet
    f = read(pV, j);
    frameArray{i} = f;
    j = j + keyFrameNo;
end

%% Display new video
for i = 1:framesToGet
    imshow(frameArray{i});
end

%% saving images to file - uncomment to save

% variables needed for imcrop
xmin = 100;
ymin = 35;
width = 220;
height = 170;
rect = [xmin ymin width height];
    
for i = 1:framesToGet
    % cropping for barret videos (to crop out the black border)
    % make sure to comment out for videos other than the barret videos
    %frameArray{i} = imcrop(frameArray{i},rect);
    % write to file
    %imwrite(frameArray{i}, sprintf('barret1/im%d.jpeg', i));
    imwrite(frameArray{i}, sprintf('tiger/im%d.jpeg', i));

end




