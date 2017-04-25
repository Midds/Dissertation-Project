%% TO RUN, VLFEAT MUST FIRST BE INSTALLED ON THE MACHINE
% VLFEAT can be downloaded from http://www.vlfeat.org/download.html or http://www.vlfeat.org/index.html
% once downloaded and unpacked the command below must be ran on each Matlab restart
% run D:\Users\James\Documents\GitHub\ImageStitching\MATLAB\vlfeat-0.9.20/toolbox/vl_setup
% - with the pathway changed to match the vl_setup path
close all;

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
        else
            Ig = I;
        end
        
        fprintf('pre-processing done for image: %d \n', n);
        fprintf('');
        
        %% FINDING SIFT FEATURES AND DESCRIPTORS
        
        [F2,D2] = vl_sift(Ig);
        
        fprintf('sift features found for image: %d\n', n);
        
        % MATCH LOCAL DESCRIPTORS
        % for each descriptor(key feature) in d1, vl_ubcmatch finds the closest descriptor in d2
        % the index gets stored in matches and the distance between in scores.
        % ubcmatch computes approximate matches between the images - so it's
        % fast but could be more accurate
        [matches, scores] = vl_ubcmatch(D1, D2);
        
        numMatches = size(matches,2);
        fprintf('matched %d local descriptors for image: %d and %d\n', numMatches, m-1, m);
        
        
        %% RANSAC
        % Ransac computes a homography matrix that can then be used to map
        % the coordinates of one image to the coordinates of another image
        X1 = F1(1:2,matches(1,:)) ; X1(3,:) = 1 ;
        X2 = F2(1:2,matches(2,:)) ; X2(3,:) = 1 ;
        
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
        
        fprintf('Saving H matrix for im%d and im%d\n',m-1, m);
        
        %name = ['H', string(m-1), '_', string(m)];
        name = ['homography/H', string(m-1)];
        save(char(name) , 'H');
        
        if (m > numToStitch/2)
            m = m - 1;
        else
            m = m + 1;
        end
               
        %% Showing inliner matches
        % again - for displaying points use the original im1 and im2 not Ig
        dh1 = max(size(Ig,1)-size(imPrev,1),0) ;
        dh2 = max(size(imPrev,1)-size(Ig,1),0) ;
        
        % matches before ransac
        figure; clf ;
        subplot(2,1,1) ;
        imagesc([padarray(imPrev,dh1,'post') padarray(Ig,dh2,'post')]) ;
        o = size(imPrev,2) ;
        line([F1(1,matches(1,:));F2(1,matches(2,:))+o], ...
            [F1(2,matches(1,:));F2(2,matches(2,:))]) ;
        title(sprintf('%d tentative matches', numMatches)) ;
        axis image off ;
        
        % matches after ransac
        subplot(2,1,2) ;
        imagesc([padarray(imPrev,dh1,'post') padarray(Ig,dh2,'post')]) ;
        o = size(imPrev,2) ;
        line([F1(1,matches(1,ok));F2(1,matches(2,ok))+o], ...
            [F1(2,matches(1,ok));F2(2,matches(2,ok))]) ;
        title(sprintf('%d (%.2f%%) inliner matches out of %d', ...
            sum(ok), ...
            100*sum(ok)/numMatches, ...
            numMatches)) ;
        axis image off ;
        
        drawnow ;
        
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

 H12 = load('homography/H1           '); H12 = H12.H; % Homography of im1 to im2
 H23 = load('homography/H2           '); H23 = H23.H; % Homography of im3 to im2
 H34 = load('homography/H3           '); H34 = H34.H; % Homography of im1 to im2
 H54 = load('homography/H4           '); H54 = H54.H; % Homography of im3 to im2
 H65 = load('homography/H5           '); H65 = H65.H; % Homography of im1 to im2
 H76 = load('homography/H6           '); H76 = H76.H; % Homography of im3 to im2

fprintf('Homography matrices loaded \n');


% im1 im2 im3 im4 im5 im6 im7 <- order of single images : im4 is in center.
H14 = H12*H23*H34;
H24 = H23*H34;
H64 = H65*H54;
H74 = H76*H65*H54;

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

%% Boundary Condition of Mosaiced Image %%
h14 = H14'; h14 = h14(:); % Change homograpy to a vector form.
h24 = H24'; h24 = h24(:);
h34 = H34'; h34 = h34(:);
h54 = H54'; h54 = h54(:);
h64 = H64'; h64 = h64(:);
h74 = H74'; h74 = h74(:);

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