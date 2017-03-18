% TO RUN, VLFEAT MUST FIRST BE INSTALLED ON THE MACHINE
% VLFEAT can be downloaded from http://www.vlfeat.org/download.html or http://www.vlfeat.org/index.html
% once downloaded and unpacked the command below must be ran on each Matlab restart
% run D:\Users\James\Documents\MATLAB\vlfeat-0.9.20/toolbox/vl_setup
% - with the pathway changed to match the vl_setup path
close all;

% images to merge
%imgs = imageDatastore({'im15.jpeg';'im16.jpeg'});
imgs = imageDatastore({'a.jpg';'b.jpg'});

figure;
montage(imgs.Files);
title('Montage');

%im1 = imread('im15.jpeg');
%im2 = imread('im16.jpeg');
im1 = imread('a.jpg');
im2 = imread('b.jpg');

% PREPROCESSING
% convert to greyscale
im1 = rgb2gray(im1);
im2 = rgb2gray(im2);

% make single
% vl_feat requires single precision greyscale image
im1 = single(im1);
im2 = single(im2);

disp('pre-processing done');

% FINDING SIFT FEATURES AND DESCRIPTORS
% vl_sift uses the vlfeat open source implementation of sift to find
% keypoints [k] and descriptors [d]
[k1,d1] = vl_sift(im1);
[k2,d2] = vl_sift(im2);
disp('sift features found');

% MATCH LOCAL DESCRIPTORS
% for each descriptor(key feature) in d1, vl_ubcmatch finds the closest descriptor in d2
% the index gets stored in matches and the distance between in scores
[matches, scores] = vl_ubcmatch(d1, d2);
disp('matched local descriptors');

%% Draw Some Matching Features
npts = 10;

figure(1), colormap gray; imagesc(im1);
figure(2), colormap gray; imagesc(im2);
for i = 1 : npts
    ind1 = matches(1,i);
    ind2 = matches(2,i);

    figure(1);
    plot1 = vl_plotsiftdescriptor(d1(:, ind1),k1(:, ind1));
    set(plot1, 'color', hsv2rgb([i / npts, 1, 1]));
    
    figure(2)
    plot2 = vl_plotsiftdescriptor(d2(:, ind2),k2(:, ind2));
    set(plot2, 'color', hsv2rgb([i / npts, 1, 1]));
    
end

disp('End');

figure;
imagesc(cat(2, im1, im2)) ;
hold on;
plot (k1(1,matches(1,:)), k1(2, matches(1,:)), 'b*');

hold on ;

xa = k1(1,matches(1,:)) ;
xb = k2(1,matches(2,:)) + size(im1,2) ;
ya = k1(2,matches(1,:)) ;
yb = k2(2,matches(2,:)) ;

h = line([xa ; xb], [ya ; yb]) ;
set(h,'linewidth', 1, 'color', 'b') ;


hold on;
k2(1,:) = k2(1,:) + size(im1,2) ;
plot (k2(1, matches(2,:)), k2(2, matches (2,:)), 'r*');
axis image off ;

disp('end2');

figure; clf ;
imagesc(cat(2, im1, im2)) ;
%%%%
xa = k1(1,matches(1,:)) ;
xb = k2(1,matches(2,:)) + size(im1,2) ;
ya = k1(2,matches(1,:)) ;
yb = k2(2,matches(2,:)) ;

hold on ;
h = line([xa ; xb], [ya ; yb]) ;
set(h,'linewidth', 1, 'color', 'b') ;

vl_plotframe(k1(:,matches(1,:))) ;
k2(1,:) = k2(1,:) + size(im1,2) ;
vl_plotframe(k2(:,matches(2,:))) ;
axis image off ;
