% MATLAB script for Assessment Item-1
close all;

% Step-1: Load input image
InputImage = imread('AssignmentInput.jpg');
figure;
imshow(InputImage);
title('Step-1: Load input image');

% Step-2: Conversion of input image to greyscale
red = InputImage(:,:,1); 
green = InputImage(:,:,2);
blue = InputImage(:,:,3);
I = (0.299 * red) + (0.587 * green) + (0.114 * blue);
figure;
imshow(I);
title('Step-2: Conversion of input image to greyscale');

% Step-3: Noise removal
% Median filter performs neighbourhood operation on each pixel to change
% the value of the pixel to the median of it's neighbourhood.
mf = medfilt2(I, [3 3], 'symmetric');
figure;
imshow(mf);
title('Step-3: Median filter Noise removal - 3x3 mask');

% Step-4: Enhancing the image - justify choice

% adaptive histeq
% this is the chosen enhancement to make out of all those tested
mfadaptHE = adapthisteq(mf);
figure;
imshow(mfadaptHE);
title('Step-4: Enhanced with adaptive histogram equalisation');
%figure;
%imhist(mfadaptHE);
%title('Image enhanced with adaptive histogram equalisation');

% The other tested enhancements are commented out below should they be
% needed.

% Contrast stretching
%Im2D = 255*im2double(mf);
%mi = min(min(Im2D));
%ma = max(max(Im2D));
%contStretchedIm = imadjust(mf, [0.75; ma/255], [0;1]);
%figure;
%imshow(contStretchedIm);
%title('Step-4: Enhanced with contrast stretching');
%figure;
%imhist(contStretchedIm);
%title('Image enhanced with contrast stretching');

% Histogram equalisation
%mfhe = histeq(mf);
%figure;
%imshow(mfhe);
%title('Step-4: Enhanced with histogram equalisation');
%figure;
%imhist(mfhe);
%title('Image enhanced with histogram equalisation');

% Step-5: Segment the image into foreground and background
% Loops through image and assigns pixels to 1 or 0 (225 or 0) depending on 
% if they are above or below a certain intensity value.
bi = zeros(size(mfadaptHE));
threshold = (graythresh(mfadaptHE)) * 255;

for row = 1:size(bi, 1)
    for col = 1:size(bi, 2)
        if ((mfadaptHE(row,col)) > (threshold))
            bi(row,col) = 255;
        end
    end
end

figure;
imshow(~bi);
title('Step-5: Segmented image');

% Step-6: Use of morphological processing
se = strel('disk', 2);

biOpen = imopen(~bi, se);

% closing the already imopened image to close some holes in the starfish
biClose = imclose(biOpen, se);
figure;
imshow(biClose);
title('Step 6: Use of morphological processing - imopen + imclose');

% Step-7: Recognition of starfishes
% First labels all components in the image, then accesses the properties of
% each component using regionprops.
CC = bwconncomp(biClose);
stats = regionprops('table',CC, 'Area', 'Perimeter');
area = mean(stats.Area);
perim = mean(stats.Perimeter);

metric = 0;
E = zeros(size(biClose));  

for i = 1:CC.NumObjects
    area = stats.Area(i);
    perim = stats.Perimeter(i);   
    metric = ( (4 * pi * area) / (perim.^2) );
    
    % If the component falls within a certain range then it is a starfish.
    % This component is then written to the image (E).
    if ( (metric < 0.253) && (metric > 0.21) )  
       E(CC.PixelIdxList{i}) = 255;
    end
end

figure;
imshow(E);
title('Step 7:  Recognition of starfishes');

