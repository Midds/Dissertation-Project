%% TO RUN, VLFEAT MUST FIRST BE INSTALLED ON THE MACHINE
% VLFEAT can be downloaded from http://www.vlfeat.org/download.html or http://www.vlfeat.org/index.html
% once downloaded and unpacked the command below must be ran on each Matlab restart
% run D:\Users\James\Documents\GitHub\ImageStitching\MATLAB\vlfeat-0.9.20/toolbox/vl_setup
% - with the pathway changed to match the vl_setup path
close all;

% pipeline as follows
% - read in images to mosaic
% - find sift for each and check there's enough matches, if not read in
% another image
% - use ransac on sift matches to get homograpy
% - refine homography using levenberg-marquardt
% - stitch images

% select number of images to stitch
numToStitch = 7;
% select the starting image to stitch (whatever number it is in the file)
% eg if image name = 'im168.jpeg', then startImage = 168
startImage = 168;
% Creating an array to store the images
imArray = {};

for n = startImage:(startImage+numToStitch)-1
    filename = sprintf('tiger/im%dsmall.jpeg', n); % defining the filename
    im = imread(filename); % reading the image from the given filename
    imArray = [imArray im]; % adding the image to the image Array
end
 
figure;
newimage = cell2mat(imArray);
imshow(newimage);
title('Images to stitch');

% Read the first image from the image set.
im1 = imArray{1};

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
% features F2 and descriptors D2 based on the sift algorithm (Lowe, 2004)
[F2,D2] = vl_sift(im1g);
disp('sift features found for image: 1');

%points = detectSURFFeatures(grayImage);
%[features, points] = extractFeatures(grayImage, points);

% Now loop through the remaining images
% The nested loop will save a homography matrix (Hn_n+1) for each successive image
% pair up until the centre image of the image array. The image array is
% then reversed in order, and the nested loop will go again, looping until
% it hits the centre image. This is needed as the H matrices need to be in
% this order to successfully stitch them together later.
m = 2;
for j = 1:2
    for n = 2:(numToStitch/2) + 1
        % Store points and features for I(n-1).
    % d1,k1 will then be updated later in the loop to store n
    F1 = F2;
    D1 = D2;
    
    %% PREPROCESSING
    % read the next image from the image datastore (starts at image 2 as first image alread read in)
    I = imArray{n};
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
img = [imArray{m},imArray{n+1}];
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

        
        fprintf('Saving H matrix for im%d and im%d\n',m-1, m);
        
        name = ['H', string(m-1), '_', string(m)];
        save(char(name) , 'H');
        
        if (m > 3)
            m = m - 1;
        else
            m = m + 1;
        end
        

%         %%
%         %Showing inliner matches
%         % again - for displaying points use the original im1 and im2 not Ig
%         dh1 = max(size(Ig,1)-size(imPrev,1),0) ;
%         dh2 = max(size(imPrev,1)-size(Ig,1),0) ;
%         
%         figure; clf ;
%         subplot(2,1,1) ;
%         imagesc([padarray(imPrev,dh1,'post') padarray(Ig,dh2,'post')]) ;
%         o = size(imPrev,2) ;
%         line([F1(1,matches(1,:));F2(1,matches(2,:))+o], ...
%             [F1(2,matches(1,:));F2(2,matches(2,:))]) ;
%         title(sprintf('%d tentative matches', numMatches)) ;
%         axis image off ;
%         
%         subplot(2,1,2) ;
%         imagesc([padarray(imPrev,dh1,'post') padarray(Ig,dh2,'post')]) ;
%         o = size(imPrev,2) ;
%         line([F1(1,matches(1,ok));F2(1,matches(2,ok))+o], ...
%             [F1(2,matches(1,ok));F2(2,matches(2,ok))]) ;
%         title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
%             sum(ok), ...
%             100*sum(ok)/numMatches, ...
%             numMatches)) ;
%         axis image off ;
%         
%         drawnow ;
        
        %     %% REWORK
        %
        %     d = dist(D1',D2); % Distance between D1's column and D2's column
        %     [Y I] = min(d);
        %     count = 0; % Number of non-overlapped correspondences
        %     c1 = zeros(1,2); % Corresponding feature coordinates of im1
        %     c2 = zeros(1,2); % Corresponding feature coordinates of im2
        % %% Find correspondences between two images using Euclidean distance %%
        % img = [imArray{1},imArray{2}];
        % for k = 1:length(Y)
        %     ind = 1; % Indicator to avoid overlapped correspondences
        %     for l = 1:length(I)
        %         if l~=k && I(l)==I(k)
        %             ind = 0;
        %             break;
        %         end
        %     end
        %     if ind && Y(k) < 35 % Threshold for Euclidean distance
        %         count = count + 1;
        %         c1(count,:) = round(F1(1:2,I(k)));
        %         c2(count,:) = round(F2(1:2,k));
        %     end
        % end
        %
        %
        %     %% RANSAC algorithm %%
        %     nc = 6; % Number of correspondences used to find a homography
        %     N = fix(log(1-.99)/log(1-(1-.1)^nc)); % Number of trials by 10% rule
        %     M = fix((1-.1)*count); % Minimum size for the inlier set
        %     d_min = 1e100;
        %     for o = 1:N
        %         lcv = 1; % Loop control variable
        %         while lcv % To avoid repeated selection
        %         r = randi(count,nc,1);
        %         r = sort(r);
        %             for k = 1:nc-1
        %             lcv = lcv*(r(k+1)-r(k));
        %             end
        %             lcv = ~lcv;
        %         end
        %         A = zeros(2*nc,9);
        %         for k = 1:nc
        %             A(2*k-1:2*k,:)=...
        %                 [0,0,0,-[c1(r(k),:),1],c2(r(k),2)*[c1(r(k),:),1];
        %                 [c1(r(k),:),1],0,0,0,-c2(r(k),1)*[c1(r(k),:),1]];
        %         end
        %         [U,D,V] = svd(A);
        %         h = V(:,9);
        %         H = [h(1),h(2),h(3);h(4),h(5),h(6);h(7),h(8),h(9)];
        %
        %         d2 = zeros(count,1); % d^2(x_measured, x_true)
        %         for k = 1:count
        %             x_true = H*[c1(k,:),1]'; % x_true in HC
        %             temp = x_true/x_true(3);
        %             x_true = temp(1:2); % x_true in image plane
        %             d = c2(k,:)-x_true';
        %             d2(k) = d(1)^2+d(2)^2;
        %         end
        %         [Y I] = sort(d2);
        %         if sum(Y(1:M)) < d_min
        %             d_min = sum(Y(1:M));
        %             inliers = I(1:M);
        %             outliers = I(M+1:end);
        %         end
        %     end
        %
        %     % Visualize the inliers and outliers
        %     figure; image(img); truesize; hold on;
        %     for k = inliers'
        %         plot([c1(k,1),C+c2(k,1)],[c1(k,2),c2(k,2)],'-og','linewidth',1);
        %     end
        %     for k = outliers'
        %         plot([c1(k,1),C+c2(k,1)],[c1(k,2),c2(k,2)],'-or','linewidth',1);
        %     end
        %     plot([C,C],[1,R],'-k'); hold off;
    end
    % now flip the order of the image array
    % the inner for loop will now loop again, getting H values starting
    % from the end of the array and working towards the centre.
    imArray = fliplr(imArray);
    m=7; % m is used in the inner loop to keep track of which images to save
end


%  imArray = fliplr(imArray);
% 
%  im1 = imArray{1};
%  im2 = imArray{2};
%  im3 = imArray{3};
%  im4 = imArray{4};
%  im5 = imArray{5};
%  im6 = imArray{6};
%  im7 = imArray{7};

im1 = imread('tiger/im168small.jpeg');
im2 = imread('tiger/im169small.jpeg');
im3 = imread('tiger/im170small.jpeg');
im4 = imread('tiger/im171small.jpeg');
im5 = imread('tiger/im172small.jpeg');
im6 = imread('tiger/im173small.jpeg');
im7 = imread('tiger/im174small.jpeg');
[M,N,C] = size(im2);

fprintf('imreads done \n');

% Load estimated and refined homographies in previous steps.
% All the refined homographies were saved as mat files.
% H12 = load('H1_2'); H12 = H12.H; % Homography of im1 to im2
% H23 = load('H2_3'); H23 = H23.H; % Homography of im3 to im2
% H34 = load('H3_4'); H34 = H34.H; % Homography of im1 to im2
% H54 = load('H6_7'); H54 = H54.H; % Homography of im3 to im2
% H65 = load('H5_6'); H65 = H65.H; % Homography of im1 to im2
% H76 = load('H4_5'); H76 = H76.H; % Homography of im3 to im2

H12 = load('Test1_2'); H12 = H12.H; % Homography of im1 to im2
H23 = load('Test2_3'); H23 = H23.H; % Homography of im3 to im2
H34 = load('Test3_4'); H34 = H34.H; % Homography of im1 to im2
H54 = load('Test5_4'); H54 = H54.H; % Homography of im3 to im2
H65 = load('Test6_5'); H65 = H65.H; % Homography of im1 to im2
H76 = load('Test7_6'); H76 = H76.H; % Homography of im3 to im2
fprintf('Homography matrices loaded \n');

% im1 im2 im3 im4 im5 im6 im7 <- order of single images : im4 is in center.
H14 = H12*H23*H34;
H24 = H23*H34;
H64 = H65*H54;
H74 = H76*H65*H54;

fprintf('Homography matrices multiplied\n');

%% Boundary Condition of Mosaiced Image %%
h14 = H14'; h14 = h14(:); % Change homograpy to a vector form.
h24 = H24'; h24 = h24(:);
h34 = H34'; h34 = h34(:);
h54 = H54'; h54 = h54(:);
h64 = H64'; h64 = h64(:);
h74 = H74'; h74 = h74(:);
fprintf('Changed to vector format\n');

c14 = fun(h14,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im1
c74 = fun(h74,[1,1,N,1,1,M,N,M]); % Transformed boundaries of im7
fprintf('boundaries found\n');

x = [1,3,5,7];
y = [2,4,6,8];
xmin = round(min([c14(x);c74(x)]));
xmax = round(max([c14(x);c74(x)]));
ymin = round(min([c14(y);c74(y)]));
ymax = round(max([c14(y);c74(y)]));

fprintf('minmax boundary done\n');

%% Assign pixel values into the mosaiced image %%
img = zeros(ymax-ymin+1,xmax-xmin+1,C); % Initialize mosaiced image
fprintf('assigned zeros done\n');
img = mosaic(img,im1,H14,xmin,ymin); % Mosaicking im1
fprintf('mosaic im1\n');
img = mosaic(img,im7,H74,xmin,ymin); % Mosaicking im7
fprintf('mosaic im7\n');
img = mosaic(img,im2,H24,xmin,ymin); % Mosaicking im2
fprintf('mosaic im2\n');
img = mosaic(img,im6,H64,xmin,ymin); % Mosaicking im6
fprintf('mosaic im6\n');
img = mosaic(img,im3,H34,xmin,ymin); % Mosaicking im3
fprintf('mosaic im3\n');
img = mosaic(img,im5,H54,xmin,ymin); % Mosaicking im5
fprintf('mosaic im5\n');
img(2-ymin:M+1-ymin,2-xmin:N+1-xmin,:) = im4; % Mosaicking im4
fprintf('mosaic im4\n');
figure; imshow(img);