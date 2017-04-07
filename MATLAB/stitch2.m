%% TO RUN, VLFEAT MUST FIRST BE INSTALLED ON THE MACHINE
% VLFEAT can be downloaded from http://www.vlfeat.org/download.html or http://www.vlfeat.org/index.html
% once downloaded and unpacked the command below must be ran on each Matlab restart
% run D:\Users\James\Documents\GitHub\ImageStitching\MATLAB\vlfeat-0.9.20/toolbox/vl_setup
% - with the pathway changed to match the vl_setup path
close all;

%% Images to merge
% uncomment depending on which image you want to see
 imds = imageDatastore({'tiger/im168small.jpeg';'tiger/im169small.jpeg';'tiger/im170small.jpeg';'tiger/im171small.jpeg';'tiger/im172small.jpeg';'tiger/im173small.jpeg';'tiger/im174small.jpeg'});
% imds = imageDatastore({'tiger/tigerSmall10.jpeg';'tiger/tigerSmall11.jpeg'});
% imds = imageDatastore({'london/im15.jpeg';'london/im16.jpeg'});
% imds = imageDatastore({'a.jpg';'b.jpg'});
% imds = imageDatastore({'barret1/im67.jpeg';'barret1/im68.jpeg'});
% imds = imageDatastore({'mosaic2.jpeg';'mosaicResized.jpeg'});
% imds = imageDatastore({'barret1/im67.jpeg';'barret1/im68.jpeg'});

% imds = imageDatastore({'mosaic.jpeg';'testResized.jpeg'}); %mosaic = 67 + 68

figure;
montage(imds.Files);
title('Image Montage');
numImages = numel(imds.Files);

% Read the first image from the image set.
im1 = readimage(imds, 1);

% preprocessing for im1
im1 = im2single(im1);
% make grayscale
    if size(im1,3) > 1 
        im1g = rgb2gray(im1); 
    else
        im1g = im1;
    end

% finding sift kypoints for im1
% vl_sift uses the vlfeat open source implementation of sift to find
% keypoints(feature frame) and descriptors based on the sift algorithm (Lowe, 2004)
[keypointCurrent,descCurrent] = vl_sift(im1g);
disp('sift features found for image: 1');

%points = detectSURFFeatures(grayImage);
%[features, points] = extractFeatures(grayImage, points);



% Now loop through the remaining images
% the loop will return a homography matrix Hx for each successive image
% pair - these can then be multipled together to form the final mosaic
for n = 2:numImages
    
    % Store points and features for I(n-1).
    % d1,k1 will then be updated later in the loop to store n
    descPrev = descCurrent;
    keypointPrev = keypointCurrent;
    imPrev = im1g;
    
    %% PREPROCESSING
    % read the next image from the image datastore (starts at image 2 as first image alread read in)
    I = readimage(imds, n);
    % make single
    % vl_feat requires single precision greyscale image
    I = im2single(I);
    
    % making the image grayscale
    if size(I,3) > 1
        Ig = rgb2gray(I);
    else
        Ig = I; 
    end
    
    fprintf('pre-processing done for image: %d \n', n);
    fprintf('');

    %% FINDING SIFT FEATURES AND DESCRIPTORS
    
    [keypointCurrent,descCurrent] = vl_sift(Ig);

    fprintf('sift features found for image: %d\n', n);

    % MATCH LOCAL DESCRIPTORS
    % for each descriptor(key feature) in d1, vl_ubcmatch finds the closest descriptor in d2
    % the index gets stored in matches and the distance between in scores.
    % ubcmatch computes approximate matches between the images - so it's
    % fast but could be more accurate
    [matches, scores] = vl_ubcmatch(descPrev, descCurrent);

    numMatches = size(matches,2);
    fprintf('matched %d local descriptors for image: %d and %d\n', numMatches, n-1, n);

    %% Draw Some Matching Features
    % npts = 10;
    % 
    % figure, colormap gray; imagesc(im1);
    % figure, colormap gray; imagesc(im2);
    % for i = 1 : npts
    %     ind1 = matches(1,i);
    %     ind2 = matches(2,i);
    % 
    %     figure(1);
    %     plot1 = vl_plotsiftdescriptor(d1(:, ind1),k1(:, ind1));
    %     set(plot1, 'color', hsv2rgb([i / npts, 1, 1]));
    %     
    %     figure(2)
    %     plot2 = vl_plotsiftdescriptor(d2(:, ind2),k2(:, ind2));
    %     set(plot2, 'color', hsv2rgb([i / npts, 1, 1]));
    % end
    % disp('Matching features done');

%     %% DISPLAYING MATCHED KEYPOINTS BETWEEN IMAGES
%     figure;
%     % use im1 and im2 to display the points instead of Ig(these are the original images
%     % not greyscale versions)
%     imagesc(cat(2, im1, im2)) ;
%     hold on;
%     plot (keypointPrev(1,matches(1,:)), keypointPrev(2, matches(1,:)), 'b*');
% 
%     hold on ;
% 
%     x1 = keypointPrev(1,matches(1,:)) ;
%     x2 = keypointCurrent(1,matches(2,:)) + size(im1,2) ;
%     x3 = keypointPrev(2,matches(1,:)) ;
%     x4 = keypointCurrent(2,matches(2,:)) ;
% 
%     connectLine = line([x1 ; x2], [x3 ; x4]) ;
%     set(connectLine,'linewidth', 1, 'color', 'b') ;
% 
%     hold on;
%     temp = keypointCurrent; %k2 is used later on so temp is needed here
%     temp(1,:) = keypointCurrent(1,:) + size(im1,2) ;
%     plot (temp(1, matches(2,:)), temp(2, matches (2,:)), 'r*');
%     axis image off ;
%     title('Matched Features (SIFT)');
% 
%     fprintf('Matched features displayed\n');

    %% RANSAC
    % Ransac computes a homography matrix that can then be used to map 
    % the coordinates of one image to the coordinates of another image
    X1 = keypointPrev(1:2,matches(1,:)) ; X1(3,:) = 1 ;
    X2 = keypointCurrent(1:2,matches(2,:)) ; X2(3,:) = 1 ;

    clear H score ok ;
    for t = 1:100
      % estimate homograpyh
      subset = vl_colsubset(1:numMatches, 4) ;
      A = [] ;
      for i = subset
        A = cat(1, A, kron(X1(:,i)', vl_hat(X2(:,i)))) ;
      end
      [U,S,V] = svd(A) ;
      H{t} = reshape(V(:,9),3,3) ;

      % score homography
      X2_ = H{t} * X1 ;
      du = X2_(1,:)./X2_(3,:) - X2(1,:)./X2(3,:) ;
      dv = X2_(2,:)./X2_(3,:) - X2(2,:)./X2(3,:) ;
      ok{t} = (du.*du + dv.*dv) < 6*6 ;
      score(t) = sum(ok{t}) ;
    end

    [score, best] = max(score) ;
    H = H{best} ; % all matches
    ok = ok{best} ; % inliner matches

    fprintf('Saving H matrix for im%d and im%d\n',n-1, n);

    name = ['H', string(n-1), '_', string(n)];
    save(char(name) , 'H');
    
    %%
    %Showing inliner matches
    % again - for displaying points use the original im1 and im2 not Ig
    dh1 = max(size(Ig,1)-size(imPrev,1),0) ;
    dh2 = max(size(imPrev,1)-size(Ig,1),0) ;

    figure; clf ;
    subplot(2,1,1) ;
    imagesc([padarray(imPrev,dh1,'post') padarray(Ig,dh2,'post')]) ;
    o = size(imPrev,2) ;
    line([keypointPrev(1,matches(1,:));keypointCurrent(1,matches(2,:))+o], ...
         [keypointPrev(2,matches(1,:));keypointCurrent(2,matches(2,:))]) ;
    title(sprintf('%d tentative matches', numMatches)) ;
    axis image off ;

    subplot(2,1,2) ;
    imagesc([padarray(imPrev,dh1,'post') padarray(Ig,dh2,'post')]) ;
    o = size(imPrev,2) ;
    line([keypointPrev(1,matches(1,ok));keypointCurrent(1,matches(2,ok))+o], ...
         [keypointPrev(2,matches(1,ok));keypointCurrent(2,matches(2,ok))]) ;
    title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
                  sum(ok), ...
                  100*sum(ok)/numMatches, ...
                  numMatches)) ;
    axis image off ;

    drawnow ;

end
[M,N,C] = size(im3);

% Load estimated and refined homographies in previous steps.
% All the refined homographies were saved as mat files.
H12 = load('H1_2'); H12 = H12.H; % Homography of im1 to im2
H23 = load('H2_3'); H23 = H23.H; % Homography of im3 to im2
H34 = load('H3_4'); H34 = H34.H; % Homography of im1 to im2
H54 = load('H4_5'); H54 = H54.H; % Homography of im3 to im2
H65 = load('H5_6'); H65 = H65.H; % Homography of im1 to im2
H76 = load('H6_7'); H76 = H76.H; % Homography of im3 to im2
% im1 im2 im3 im4 im5 im6 im7 <- order of single images : im4 is in center.
H14 = H12*H23*H34;
H24 = H23*H34;
H64 = H65*H54;
H74 = H76*H65*H54;
%% Boundary Condition of Mosaiced Image %%
h14 = H14'; h14 = h14(:); % Change homograpy to a vector form.
h24 = H24'; h24 = h24(:);
h34 = H34'; h34 = h34(:);
h54 = H54'; h54 = h54(:);
h64 = H64'; h64 = h64(:);
h74 = H74'; h74 = h74(:);
c14 = fun(h14,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im1
c74 = fun(h74,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im7
x = [1,3,5,7];
y = [2,4,6,8];
xmin = round(min([c14(x);c74(x)]));
xmax = round(max([c14(x);c74(x)]));
ymin = round(min([c14(y);c74(y)]));
ymax = round(max([c14(y);c74(y)]));
%% Assign pixel values into the mosaiced image %%
img = zeros(ymax-ymin+1,xmax-xmin+1,C); % Initialize mosaiced image
img = mosaic(img,im1,H14,xmin,ymin); % Mosaicking im1
img = mosaic(img,im7,H74,xmin,ymin); % Mosaicking im7
img = mosaic(img,im2,H24,xmin,ymin); % Mosaicking im2
img = mosaic(img,im6,H64,xmin,ymin); % Mosaicking im6
img = mosaic(img,im3,H34,xmin,ymin); % Mosaicking im3
img = mosaic(img,im5,H54,xmin,ymin); % Mosaicking im5
img(2-ymin:M+1-ymin,2-xmin:N+1-xmin,:) = im4; % Mosaicking im4
figure; imshow(img); imwrite(img,'Mosaic_apt');

% % Load estimated and refined homographies in previous steps.
% % All the refined homographies were saved as mat files.
% H12 = load('H1_2'); H12 = H12.H; % Homography of im1 to im2
% H23 = load('H2_3'); H23 = H23.H; % Homography of im2 to im3
% H34 = load('H3_4'); H34 = H34.H; % Homography of im3 to im4
% H45 = load('H4_5'); H45 = H45.H; % Homography of im4 to im5
% 
% % im1 im2 im3 im4 im5 <- order of single images : im3 is in center.
% H13 = H12*H23;
% H35 = H45*H34;
% 
% %% Boundary Condition of Mosaiced Image %%
% h13 = H13'; h13 = h13(:); % Change homograpy to a vector form.
% h23 = H23'; h23 = h23(:);
% h34 = H34'; h34 = h34(:);
% h35 = H35'; h35 = h35(:);
% 
% c14 = fun(h13,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im1
% c74 = fun(h35,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im7
% x = [1,3,5,7];
% y = [2,4,6,8];
% xmin = round(min([c14(x);c74(x)]));
% xmax = round(max([c14(x);c74(x)]));
% ymin = round(min([c14(y);c74(y)]));
% ymax = round(max([c14(y);c74(y)]));
% %% Assign pixel values into the mosaiced image %%
% img = zeros(ymax-ymin+1,xmax-xmin+1,C); % Initialize mosaiced image
% img = mosaic(img,im1,H13,xmin,ymin); % Mosaicking im1
% img = mosaic(img,im5,H35,xmin,ymin); % Mosaicking im5
% img = mosaic(img,im2,H23,xmin,ymin); % Mosaicking im2
% img = mosaic(img,im4,H34,xmin,ymin); % Mosaicking im4
% 
% img(2-ymin:M+1-ymin,2-xmin:N+1-xmin,:) = im3; % Mosaicking im4
% figure; imshow(img); imwrite(img,'Mosaic_apt');


% %% MOSAICING
% % Works by transforming the first image onto the plane of the second,
% % before stiching them together. Images are stitched by transforming 
% % pixel coordinates of the first image to the pixel coordinates of the 
% % second image plane by multiplying with the homography matrix H
% 
% % box2 is a 2d matrix that holds the corner coordinates of the second image
% box2 = [1  size(im2,2) size(im2,2)  1           ;
%         1  1           size(im2,1)  size(im2,1) ;
%         1  1           1            1 ]         ;
% %box2_ is box2 transformed into the coordinates of the first image
% box2_ = inv(H) * box2 ;
% box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
% box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
% % ur and vr are the range of coordinates of where to project im2
% ur = min([1 box2_(1,:)]):max([size(im1,2) box2_(1,:)]) ;
% vr = min([1 box2_(2,:)]):max([size(im1,1) box2_(2,:)]) ;
% 
% [u,v] = meshgrid(ur,vr) ;
% im1_ = vl_imwbackward(im2double(im1),u,v) ;
% 
% z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
% u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
% v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
% im2_ = vl_imwbackward(im2double(im2),u_,v_) ;
% % mass is used to find out how many images cover each pixel 
% mass = ~isnan(im1_) + ~isnan(im2_) ;
% im1_(isnan(im1_)) = 0 ;
% im2_(isnan(im2_)) = 0 ;
% mosaic = (im1_ + im2_) ./ mass ;
% 
% fprintf('Mosaicing done\n');
% 
% figure ; clf ;
% imagesc(mosaic) ; axis image off ;
% title('Mosaic') ;
% 
% imwrite(mosaic, sprintf('mosaic.jpeg'));