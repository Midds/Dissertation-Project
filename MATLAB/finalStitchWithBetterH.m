%% TO RUN, VLFEAT MUST FIRST BE INSTALLED ON THE MACHINE
% VLFEAT can be downloaded from http://www.vlfeat.org/download.html or http://www.vlfeat.org/index.html
% once downloaded and unpacked the command below must be ran on each Matlab restart
% run D:\Users\James\Documents\GitHub\ImageStitching\MATLAB\vlfeat-0.9.20/toolbox/vl_setup
% - with the pathway changed to match the vl_setup path
%close all;

% pipeline as follows
% - read in images to mosaic
% - find sift for each and check there's enough matches, if not read in
% another image
% - use ransac on sift matches to get homograpy
% - refine homography using levenberg-marquardt
% - stitch images

%line 111 change threshold to 25 for barett images

% select number of images to stitch
numToStitch = 7;
% select the starting image to stitch (whatever number it is in the file)
% eg if image name = 'im168.jpeg', then startImage = 168
startImage = 168;
% Creating an array to store the images
imArray = {};

for n = startImage:(startImage+numToStitch)-1
    filename = sprintf('barret1/im%d.jpeg', n); % defining the filename
    im = imread(filename); % reading the image from the given filename
    im = imresize(im,1.5); %- required for barret1 images (change the number to 2)
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
        Ig = rgb2gray(im1); 
        Ig = imadjust(Ig);

    else
        Ig = im1;
    end

% finding sift kypoints for im1
% vl_sift uses the vlfeat open source implementation of sift to find
% features F2 and descriptors D2 based on the sift algorithm (Lowe, 2004)
[F2,D2] = vl_sift(Ig);
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
        %Store points and features for I(n-1).
        % d1,k1 will then be updated later in the loop to store n
        D1 = D2;
        F1 = F2;
        
        imPrev = Ig;
        
        %% PREPROCESSING
        % read the next image from the image array (starts at image 2 as first image alread read in)
        I = imArray{n};
        % make single
        % vl_feat requires single precision greyscale image
        I = im2single(I);
        
        % making the image grayscale
        if size(I,3) > 1
            Ig = rgb2gray(I);
            Ig = imadjust(Ig);

        else
            Ig = I;
        end
        
        fprintf('pre-processing done for image: %d \n', n);
        fprintf('');
        
        %% FINDING SIFT FEATURES AND DESCRIPTORS
        
        [F2,D2] = vl_sift(Ig);
        
        fprintf('sift features found for image: %d\n', n);
        
        d = dist(D1',D2); % Distance between D1's column and D2's column
    [Y I] = min(d);
    count = 0; % Number of non-overlapped correspondences
    c1 = zeros(1,2); % Corresponding feature coordinates of im1
    c2 = zeros(1,2); % Corresponding feature coordinates of im2
    %% Find correspondences between two images using Euclidean distance %%
    img = [imArray{n-1},imArray{n}];
    for k = 1:length(Y)
        ind = 1; % Indicator to avoid overlapped correspondences
        for l = 1:length(I)
            if l~=k && I(l)==I(k)
                ind = 0;
                break;
            end
        end
        if ind && Y(k) < 55 % Threshold for Euclidean distance
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
        
        %name = ['H', string(m-1), '_', string(m)];
        name = ['homography/H', string(m-1)];
        save(char(name) , 'H');
        
        if (m > numToStitch/2)
            m = m - 1;
        else
            m = m + 1;
        end
               
%         %% Showing inliner matches
%         % again - for displaying points use the original im1 and im2 not Ig
%         dh1 = max(size(Ig,1)-size(imPrev,1),0) ;
%         dh2 = max(size(imPrev,1)-size(Ig,1),0) ;
%         
%         % matches before ransac
%         figure; clf ;
%         subplot(2,1,1) ;
%         imagesc([padarray(imPrev,dh1,'post') padarray(Ig,dh2,'post')]) ;
%         o = size(imPrev,2) ;
%         line([F1(1,matches(1,:));F2(1,matches(2,:))+o], ...
%             [F1(2,matches(1,:));F2(2,matches(2,:))]) ;
%         title(sprintf('%d tentative matches', numMatches)) ;
%         axis image off ;
%         
%         % matches after ransac
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
        
    end
    % now flip the order of the image array
    % the inner for loop will now loop again, getting H values starting
    % from the end of the array and working towards the centre.
    imArray = fliplr(imArray);
    
     % Read the first image from the image set.
    im1 = imArray{1};
    
    % preprocessing for im1
    im1 = im2single(im1);
    % make grayscale
    if size(im1,3) > 1
        Ig = rgb2gray(im1);
    else
        Ig = im1;
    end
    
    % finding sift kypoints for im1
    % vl_sift uses the vlfeat open source implementation of sift to find
    % features F2 and descriptors D2 based on the sift algorithm (Lowe, 2004)
    [F2,D2] = vl_sift(Ig);
    disp('sift features found for image: 1');
    
    
    m=numToStitch; % m is used in the inner loop to keep track of which images to save
end

% im1 = imread('london2/im168.jpeg');
% im2 = imread('london2/im169.jpeg');
% im3 = imread('london2/im170.jpeg');
% im4 = imread('london2/im171.jpeg');
% im5 = imread('london2/im172.jpeg');
% im6 = imread('london2/im173.jpeg');
% im7 = imread('london2/im174.jpeg');
[M,N,C] = size(imArray{2});

fprintf('imreads done \n');

% Load estimated and refined homographies in previous steps.
% All the refined homographies were saved as mat files.
% preallocating array to hold H
% HxxArray = cell(1,(numToStitch-1));
HArray{1, (numToStitch-1)} = [];
% the above loop will always save 1 less Homography than numImages, so this
% will loop from 1:numToStitch-1.
for n = 1:(numToStitch-1)
    name = ['homography/H', string(n)];
    HArray{n} = load(char(name)); HArray{n} = HArray{n}.H;
end

%  H12 = load('homography/H1           '); H12 = H12.H; % Homography of im1 to im2
%  H23 = load('homography/H2           '); H23 = H23.H; % Homography of im3 to im2
%  H34 = load('homography/H3           '); H34 = H34.H; % Homography of im1 to im2
%  H54 = load('homography/H4           '); H54 = H54.H; % Homography of im3 to im2
%  H65 = load('homography/H5           '); H65 = H65.H; % Homography of im1 to im2
%  H76 = load('homography/H6           '); H76 = H76.H; % Homography of im3 to im2

fprintf('Homography matrices loaded \n');


% % im1 im2 im3 im4 im5 im6 im7 <- order of single images : im4 is in center.
% H14 = H12*H23*H34;
% H24 = H23*H34;
% H64 = H65*H54;
% H74 = H76*H65*H54;

% there is always 3 less homographies to create than numImages
% 2 arrays are needed, one for the inwards, one for outwards
%HxArray{1, (numToStitch-3)/2} = []; % creates an empty (1 x numToStich-3) cell array 
HxArray = cell(1,(numToStitch-3)/2);
HxArray2 = cell(1,(numToStitch-3)/2);

%HxArray2{1, (numToStitch-3)/2} = []; % creates an empty (1 x numToStich-3) cell array 

m = numToStitch/2;
o = 1;
for n = 1:(numToStitch-3)/2
    HxArray{n} = 1;
    for i = o:m
        fprintf('Going once %d\n', i);
        HxArray{n} = HxArray{n} * HArray{i};
        %HxArray{n} = HArray{1} * HArray{2} * HArray{3};
    end
    %m = m - 1;   

    o = o + 1;


    fprintf('\nFirst loop ended');
        
end
HArray = fliplr(HArray);
o = 1;
for n = 1:(numToStitch-3)/2
    HxArray2{n} = 1;
    for i = o:m
        fprintf('Going once %d\n', i);
        HxArray2{n} = HxArray2{n} * HArray{i};
        %HxArray{n} = HArray{1} * HArray{2} * HArray{3};
    end
    %m = m - 1;   

    o = o + 1;


    fprintf('\nFirst loop ended');
        
end

% flip HxArray2 so it's in the right order
HxArray2 = fliplr(HxArray2);

% add the two unmultipled H matrices in the correct places
HxArray = [HxArray,HArray{(size(HArray, 2)/2)+1}];
HxArray2 = [HArray{(size(HArray, 2)/2)}, HxArray2];

% concatenating the 2 arrays into one
HxFinal = [HxArray, HxArray2];

% Each cell in the cell array holds the homograpy matrix of {n} to the
% centre image. Eg.. H14, H24, ...Hn4
% Example: HxFinal now looks like this, (assuming 7 images being mosaiced)
% HxFinal{1} = H12*H23*H34;
% HxFinal{2} = H23*H34;
% HxFinal{3} = H34;
% HxFinal{4} = H54;
% HxFinal{5} = H65*H54;
% HxFinal{6} = H76*H65*H54;


fprintf('Homography matrices multiplied\n');

% %% Boundary Condition of Mosaiced Image %%
% h14 = H14'; h14 = h14(:); % Change homograpy to a vector form.
% h24 = H24'; h24 = h24(:);
% h34 = H34'; h34 = h34(:);
% h54 = H54'; h54 = h54(:);
% h64 = H64'; h64 = h64(:);
% h74 = H74'; h74 = h74(:);

HxxArray{1, (numToStitch-1)} = [];
% HxArray =
%14
%24
%74
%64

% Change homograpy to a vector form.
% for n = 1:numToStitch-1
%    HxFinal{n} = HxFinal{n}'; 
%    HxFinal{n} = HxFinal{n}(:);     
% end
% fprintf('Changed to vector format\n');
tempHFirst = HxFinal{1}'; tempHFirst = tempHFirst(:); %Change first homograpy to a vector form.
tempHFinal = HxFinal{size(HxFinal, 2)}'; tempHFinal = tempHFinal(:); %Change last homograpy to a vector form.

c14 = fun(tempHFirst,[1,1,N,1,1,M,N,M]); % Transformed boundaries of the first image
c74 = fun(tempHFinal,[1,1,N,1,1,M,N,M]); % Transformed boundaries of the last image
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
direction = 'forward';
m = numToStitch;
n = 1;
for i = 1:numToStitch
    if (strcmp(direction ,'forward') == 1)
        %img = mosaic(img,imArray{n},HxFinal{n},xmin,ymin);
        img = mosaic(img,imArray{n},HxFinal{n},xmin,ymin);

        fprintf('mosaic im%d\n', n);
        direction = 'backward';
        n = n + 1;
    else
        img = mosaic(img,imArray{m},HxFinal{m-1},xmin,ymin);
        fprintf('mosaic im%d\n', m);
        direction = 'forward';
        m = m - 1;
    end  
end

img(2-ymin:M+1-ymin,2-xmin:N+1-xmin,:) = imArray{ ((size(imArray,2)/2)+0.5) }; % Mosaicking last image
fprintf('mosaic im4\n');
figure; imshow(img); title('final');

% %% Assign pixel values into the mosaiced image %%
% img = zeros(ymax-ymin+1,xmax-xmin+1,C); % Initialize mosaiced image
% fprintf('assigned zeros done\n');
% figure; imshow(img);
% img = mosaic(img,im1,H14,xmin,ymin); % Mosaicking im1
% fprintf('mosaic im1\n');
% figure; imshow(img); title('im1');
% img = mosaic(img,im7,H74,xmin,ymin); % Mosaicking im7
% fprintf('mosaic im7\n');
% figure; imshow(img); title('im7');
% img = mosaic(img,im2,H24,xmin,ymin); % Mosaicking im2
% fprintf('mosaic im2\n');
% figure; imshow(img); title('im2');
% img = mosaic(img,im6,H64,xmin,ymin); % Mosaicking im6
% fprintf('mosaic im6\n');
% figure; imshow(img); title('im6');
% img = mosaic(img,im3,H34,xmin,ymin); % Mosaicking im3
% fprintf('mosaic im3\n');
% figure; imshow(img); title('im3');
% img = mosaic(img,im5,H54,xmin,ymin); % Mosaicking im5
% fprintf('mosaic im5\n');
% figure; imshow(img); title('im5');
% 
% img(2-ymin:M+1-ymin,2-xmin:N+1-xmin,:) = im4; % Mosaicking im4
% fprintf('mosaic im4\n');
% figure; imshow(img); title('final');