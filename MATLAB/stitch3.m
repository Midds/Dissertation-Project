%% TO RUN, VLFEAT MUST FIRST BE INSTALLED ON THE MACHINE
% VLFEAT can be downloaded from http://www.vlfeat.org/download.html or http://www.vlfeat.org/index.html
% once downloaded and unpacked the command below must be ran on each Matlab restart
% run D:\Users\James\Documents\GitHub\ImageStitching\MATLAB\vlfeat-0.9.20/toolbox/vl_setup
% - with the pathway changed to match the vl_setup path
close all;

%% Images to merge
% uncomment depending on which image you want to see
% imds = imageDatastore({'tiger/im168small.jpeg';'tiger/im169small.jpeg';'tiger/im170small.jpeg';'tiger/im171small.jpeg';'tiger/im172small.jpeg';'tiger/im173small.jpeg';'tiger/im174small.jpeg'});
 imds = imageDatastore({'barret2/im170.jpeg';'barret2/im171.jpeg';'barret2/im172.jpeg';'barret2/im173.jpeg';'barret2/im174.jpeg';'barret2/im175.jpeg';'barret2/im176.jpeg'});

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
im1g = single(rgb2gray(im1));

% finding sift kypoints for im1
% vl_sift uses the vlfeat open source implementation of sift to find
% keypoints(feature frame) and descriptors based on the sift algorithm (Lowe, 2004)
[F2,D2] = vl_sift(im1g);
disp('sift features found for image: 1');


%% Preprocess %%
% Read images
% im2 = imread('tiger/im169small.jpeg');
% im3 = imread('tiger/im170small.jpeg');
% im4 = imread('tiger/im171small.jpeg');
% im5 = imread('tiger/im172small.jpeg');
% im6 = imread('tiger/im173small.jpeg');
% im7 = imread('tiger/im174small.jpeg');

 im1 = imresize(imread('barret2/im170.jpeg'), 2);
 im2 = imresize(imread('barret2/im171.jpeg'), 2);
 im3 = imresize(imread('barret2/im172.jpeg'), 2);
 im4 = imresize(imread('barret2/im173.jpeg'), 2);
 im5 = imresize(imread('barret2/im174.jpeg'), 2);
 im6 = imresize(imread('barret2/im175.jpeg'), 2);
 im7 = imresize(imread('barret2/im176.jpeg'), 2);

% Convert RGB image to grayscale
im2g = single(rgb2gray(im2));
im3g = single(rgb2gray(im3));
im4g = single(rgb2gray(im4));
im5g = single(rgb2gray(im5));
im6g = single(rgb2gray(im6));
im7g = single(rgb2gray(im7));
[R,C] = size(im1g);

% Now loop through the remaining images
% the loop will return a homography matrix Hx for each successive image
% pair - these can then be multipled together to form the final mosaic
for n = 2:numImages
  
    % Store points and features for I(n-1).
    % d1,k1 will then be updated later in the loop to store n
    F1 = F2;
    D1 = D2;
    
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
    

    %% FINDING SIFT FEATURES AND DESCRIPTORS
    

    
    
    
    %% RANSAC
    %% SIFT of im1 and im2 by VL_SIFT function %%
    [F2,D2] = vl_sift(Ig); % Each column of D is a discriptor
    fprintf('sift features found for image: %d\n', n);
    
    d = dist(D1',D2); % Distance between D1's column and D2's column
    [Y I] = min(d);
    count = 0; % Number of non-overlapped correspondences
    c1 = zeros(1,2); % Corresponding feature coordinates of im1
    c2 = zeros(1,2); % Corresponding feature coordinates of im2
%% Find correspondences between two images using Euclidean distance %%
img = [im1,im2];
for k = 1:length(Y)
    ind = 1; % Indicator to avoid overlapped correspondences
    for l = 1:length(I)
        if l~=k && I(l)==I(k)
            ind = 0;
            break;
        end
    end
    if ind && Y(k) < 35 % Threshold for Euclidean distance
        count = count + 1;
        c1(count,:) = round(F1(1:2,I(k)));
        c2(count,:) = round(F2(1:2,k));
    end
end
%% RANSAC algorithm %%
nc = 6; % Number of correspondences used to find a homography
N = fix(log(1-.99)/log(1-(1-.1)^nc)); % Number of trials by 10% rule
M = fix((1-.1)*count); % Minimum size for the inlier set
d_min = 1e100;
for o = 1:N
    lcv = 1; % Loop control variable
    while lcv % To avoid repeated selection
        r = randi(count,nc,1);
        r = sort(r);
        for k = 1:nc-1
            lcv = lcv*(r(k+1)-r(k));
        end
        lcv = ~lcv;
    end
    A = zeros(2*nc,9);
    for k = 1:nc
        A(2*k-1:2*k,:)=...
            [0,0,0,-[c1(r(k),:),1],c2(r(k),2)*[c1(r(k),:),1];
            [c1(r(k),:),1],0,0,0,-c2(r(k),1)*[c1(r(k),:),1]];
    end
    [U,D,V] = svd(A);
    h = V(:,9);
    H = [h(1),h(2),h(3);h(4),h(5),h(6);h(7),h(8),h(9)];
    
    d2 = zeros(count,1); % d^2(x_measured, x_true)
    for k = 1:count
        x_true = H*[c1(k,:),1]'; % x_true in HC
        temp = x_true/x_true(3);
        x_true = temp(1:2); % x_true in image plane
        d = c2(k,:)-x_true';
        d2(k) = d(1)^2+d(2)^2;
    end
    [Y I] = sort(d2);
    if sum(Y(1:M)) < d_min
        d_min = sum(Y(1:M));
        inliers = I(1:M);
        outliers = I(M+1:end);
    end
end
% Visualize the inliers and outliers
figure; image(img); truesize; hold on;
for k = inliers'
    plot([c1(k,1),C+c2(k,1)],[c1(k,2),c2(k,2)],'-og','linewidth',1);
end
for k = outliers'
    plot([c1(k,1),C+c2(k,1)],[c1(k,2),c2(k,2)],'-or','linewidth',1);
end
plot([C,C],[1,R],'-k'); hold off;
%% Linear Least Squares %%
A = zeros(2*M,9);
for k = 1:M
    A(2*k-1:2*k,:)=...
        [0,0,0,-[c1(inliers(k),:),1],c2(inliers(k),2)*[c1(inliers(k),:),1];
        [c1(inliers(k),:),1],0,0,0,-c2(inliers(k),1)*[c1(inliers(k),:),1]];
end
[U,D,V] = svd(A);
h1 = V(:,9); % Homography estimated by LLS with all inliers
%% Non-linear Least Square (Levenberg-Marquardt) %%
c1 = c1(inliers,:)';
c1 = c1(:);
c2 = c2(inliers,:)';
c2 = c2(:);
opt = optimset('Algorithm','levenberg-marquardt');
h2 = lsqcurvefit(@fun,h1,c1,c2,[],[],opt); % Refined homography by L.M.
H = [h2(1),h2(2),h2(3);h2(4),h2(5),h2(6);h2(7),h2(8),h2(9)];

    fprintf('Saving H matrix for im%d and im%d\n',n-1, n);

    name = ['Z', string(n-1), '-', string(n)];
    save(char(name) , 'H');
    
%     %%
%     %Showing inliner matches
%     % again - for displaying points use the original im1 and im2 not Ig
%     dh1 = max(size(Ig,1)-size(imPrev,1),0) ;
%     dh2 = max(size(imPrev,1)-size(Ig,1),0) ;
% 
%     figure; clf ;
%     subplot(2,1,1) ;
%     imagesc([padarray(imPrev,dh1,'post') padarray(Ig,dh2,'post')]) ;
%     o = size(imPrev,2) ;
%     line([keypointPrev(1,matches(1,:));keypointCurrent(1,matches(2,:))+o], ...
%          [keypointPrev(2,matches(1,:));keypointCurrent(2,matches(2,:))]) ;
%     title(sprintf('%d tentative matches', numMatches)) ;
%     axis image off ;
% 
%     subplot(2,1,2) ;
%     imagesc([padarray(imPrev,dh1,'post') padarray(Ig,dh2,'post')]) ;
%     o = size(imPrev,2) ;
%     line([keypointPrev(1,matches(1,ok));keypointCurrent(1,matches(2,ok))+o], ...
%          [keypointPrev(2,matches(1,ok));keypointCurrent(2,matches(2,ok))]) ;
%     title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
%                   sum(ok), ...
%                   100*sum(ok)/numMatches, ...
%                   numMatches)) ;
%     axis image off ;
% 
%     drawnow ;

end
% [M,N,C] = size(im3);
% 
% % Load estimated and refined homographies in previous steps.
% % All the refined homographies were saved as mat files.
% H12 = load('H1_2'); H12 = H12.H; % Homography of im1 to im2
% H23 = load('H2_3'); H23 = H23.H; % Homography of im3 to im2
% H34 = load('H3_4'); H34 = H34.H; % Homography of im1 to im2
% H54 = load('H4_5'); H54 = H54.H; % Homography of im3 to im2
% H65 = load('H5_6'); H65 = H65.H; % Homography of im1 to im2
% H76 = load('H6_7'); H76 = H76.H; % Homography of im3 to im2
% % im1 im2 im3 im4 im5 im6 im7 <- order of single images : im4 is in center.
% H14 = H12*H23*H34;
% H24 = H23*H34;
% H64 = H65*H54;
% H74 = H76*H65*H54;
% %% Boundary Condition of Mosaiced Image %%
% h14 = H14'; h14 = h14(:); % Change homograpy to a vector form.
% h24 = H24'; h24 = h24(:);
% h34 = H34'; h34 = h34(:);
% h54 = H54'; h54 = h54(:);
% h64 = H64'; h64 = h64(:);
% h74 = H74'; h74 = h74(:);
% c14 = fun(h14,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im1
% c74 = fun(h74,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im7
% x = [1,3,5,7];
% y = [2,4,6,8];
% xmin = round(min([c14(x);c74(x)]));
% xmax = round(max([c14(x);c74(x)]));
% ymin = round(min([c14(y);c74(y)]));
% ymax = round(max([c14(y);c74(y)]));
% %% Assign pixel values into the mosaiced image %%
% img = zeros(ymax-ymin+1,xmax-xmin+1,C); % Initialize mosaiced image
% img = mosaic(img,im1,H14,xmin,ymin); % Mosaicking im1
% img = mosaic(img,im7,H74,xmin,ymin); % Mosaicking im7
% img = mosaic(img,im2,H24,xmin,ymin); % Mosaicking im2
% img = mosaic(img,im6,H64,xmin,ymin); % Mosaicking im6
% img = mosaic(img,im3,H34,xmin,ymin); % Mosaicking im3
% img = mosaic(img,im5,H54,xmin,ymin); % Mosaicking im5
% img(2-ymin:M+1-ymin,2-xmin:N+1-xmin,:) = im4; % Mosaicking im4
% figure; imshow(img); imwrite(img,'Mosaic_apt');
% 
% % % Load estimated and refined homographies in previous steps.
% % % All the refined homographies were saved as mat files.
% % H12 = load('H1_2'); H12 = H12.H; % Homography of im1 to im2
% % H23 = load('H2_3'); H23 = H23.H; % Homography of im2 to im3
% % H34 = load('H3_4'); H34 = H34.H; % Homography of im3 to im4
% % H45 = load('H4_5'); H45 = H45.H; % Homography of im4 to im5
% % 
% % % im1 im2 im3 im4 im5 <- order of single images : im3 is in center.
% % H13 = H12*H23;
% % H35 = H45*H34;
% % 
% % %% Boundary Condition of Mosaiced Image %%
% % h13 = H13'; h13 = h13(:); % Change homograpy to a vector form.
% % h23 = H23'; h23 = h23(:);
% % h34 = H34'; h34 = h34(:);
% % h35 = H35'; h35 = h35(:);
% % 
% % c14 = fun(h13,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im1
% % c74 = fun(h35,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im7
% % x = [1,3,5,7];
% % y = [2,4,6,8];
% % xmin = round(min([c14(x);c74(x)]));
% % xmax = round(max([c14(x);c74(x)]));
% % ymin = round(min([c14(y);c74(y)]));
% % ymax = round(max([c14(y);c74(y)]));
% % %% Assign pixel values into the mosaiced image %%
% % img = zeros(ymax-ymin+1,xmax-xmin+1,C); % Initialize mosaiced image
% % img = mosaic(img,im1,H13,xmin,ymin); % Mosaicking im1
% % img = mosaic(img,im5,H35,xmin,ymin); % Mosaicking im5
% % img = mosaic(img,im2,H23,xmin,ymin); % Mosaicking im2
% % img = mosaic(img,im4,H34,xmin,ymin); % Mosaicking im4
% % 
% % img(2-ymin:M+1-ymin,2-xmin:N+1-xmin,:) = im3; % Mosaicking im4
% % figure; imshow(img); imwrite(img,'Mosaic_apt');
% 
% 
% % %% MOSAICING
% % % Works by transforming the first image onto the plane of the second,
% % % before stiching them together. Images are stitched by transforming 
% % % pixel coordinates of the first image to the pixel coordinates of the 
% % % second image plane by multiplying with the homography matrix H
% % 
% % % box2 is a 2d matrix that holds the corner coordinates of the second image
% % box2 = [1  size(im2,2) size(im2,2)  1           ;
% %         1  1           size(im2,1)  size(im2,1) ;
% %         1  1           1            1 ]         ;
% % %box2_ is box2 transformed into the coordinates of the first image
% % box2_ = inv(H) * box2 ;
% % box2_(1,:) = box2_(1,:) ./ box2_(3,:) ;
% % box2_(2,:) = box2_(2,:) ./ box2_(3,:) ;
% % % ur and vr are the range of coordinates of where to project im2
% % ur = min([1 box2_(1,:)]):max([size(im1,2) box2_(1,:)]) ;
% % vr = min([1 box2_(2,:)]):max([size(im1,1) box2_(2,:)]) ;
% % 
% % [u,v] = meshgrid(ur,vr) ;
% % im1_ = vl_imwbackward(im2double(im1),u,v) ;
% % 
% % z_ = H(3,1) * u + H(3,2) * v + H(3,3) ;
% % u_ = (H(1,1) * u + H(1,2) * v + H(1,3)) ./ z_ ;
% % v_ = (H(2,1) * u + H(2,2) * v + H(2,3)) ./ z_ ;
% % im2_ = vl_imwbackward(im2double(im2),u_,v_) ;
% % % mass is used to find out how many images cover each pixel 
% % mass = ~isnan(im1_) + ~isnan(im2_) ;
% % im1_(isnan(im1_)) = 0 ;
% % im2_(isnan(im2_)) = 0 ;
% % mosaic = (im1_ + im2_) ./ mass ;
% % 
% % fprintf('Mosaicing done\n');
% % 
% % figure ; clf ;
% % imagesc(mosaic) ; axis image off ;
% % title('Mosaic') ;
% % 
% % imwrite(mosaic, sprintf('mosaic.jpeg'));