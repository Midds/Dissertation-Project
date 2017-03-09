
% Loading input video file
pV = VideoReader('9poster4.mp4');

% PRE-PROCESSING
% Splitting video into a number of keyframes (to emulate lower framerate).
% This takes some compuational time, but is only necessary for this
% emulated endoscopy, so the time taken is irrelevant to final product.

% framesToGet will give the final number of frames used.
% keyFrameNo is used to decide framesToGet and is used in the for loop.
% A higher keyFrameNo means less frames in the final video.
totalFrames = pV.NumberOfFrames;
keyFrameNo = 20;
framesToGet = round(totalFrames / keyFrameNo);

% Looping through the video picking out frames at the specified interval
j = 1; % Count variable for loop
frameArray = {0:framesToGet}; % Create an array to hold the selected frames

for i = 1:framesToGet
    f = read(pV, j);
    frameArray{i} = f;
    j = j + keyFrameNo;
end

% Display new video
for i = 1:framesToGet
    imshow(frameArray{i});
end




