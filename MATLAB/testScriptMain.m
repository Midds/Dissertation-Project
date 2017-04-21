clear; close;
%% Preprocess %%
% Read images
% im1 = imread('tiger/im168small.jpeg');
% im1 = imread('tiger/im169small.jpeg');
% im1 = imread('tiger/im170small.jpeg');
% im2 = imread('tiger/im171small.jpeg');
% im2 = imread('tiger/im172small.jpeg');
 im2 = imread('tiger/im173small.jpeg');
 im1 = imread('tiger/im174small.jpeg');

% im1 = imread('barret1/im14.jpeg');
% barret14 doesnt work

% im1 = imread('london2/im370.jpeg');
% im1 = imread('london2/im371.jpeg');
% im1 = imread('london2/im372.jpeg');
% im2 = imread('london2/im373.jpeg');
% im2 = imread('london2/im374.jpeg');
% im2 = imread('london2/im375.jpeg');
% im1 = imread('london2/im376.jpeg');

% im1 = imresize(imread('barret2/im170.jpeg'), 2);
% im2 = imresize(imread('barret2/im171.jpeg'), 2);
% im1 = imread('london2/im372.jpeg');
% im2 = imread('london2/im373.jpeg');
% im2 = imread('london2/im374.jpeg');
% im2 = imread('london2/im375.jpeg');
% im1 = imread('london2/im376.jpeg');



% Convert RGB image to grayscale
im1g = single(rgb2gray(im1));
im2g = single(rgb2gray(im2));
[R,C] = size(im1g);
%% SIFT of im1 and im2 by VL_SIFT function %%
[F1,D1] = vl_sift(im1g); % Each column of F is a feature frame
[F2,D2] = vl_sift(im2g); % Each column of D is a discriptor
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
for n = 1:N
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

fprintf('Saving H matrix');
        
        name = ['Test7_6'];
        save(char(name) , 'H');