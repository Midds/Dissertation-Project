close all;
%% Preprocess %%
% Load single images.
% im1 = imread('tiger/im168small.jpeg');
% im2 = imread('tiger/im169small.jpeg');
% im3 = imread('tiger/im170small.jpeg');
% im4 = imread('tiger/im171small.jpeg');
% im5 = imread('tiger/im172small.jpeg');
% im6 = imread('tiger/im173small.jpeg');
% im7 = imread('tiger/im174small.jpeg');

 im1 = imread('london2/im370.jpeg');
 im2 = imread('london2/im371.jpeg');
 im3 = imread('london2/im372.jpeg');
 im4 = imread('london2/im373.jpeg');
 im5 = imread('london2/im374.jpeg');
 im6 = imread('london2/im375.jpeg');
 im7 = imread('london2/im376.jpeg');
[M,N,C] = size(im2);

fprintf('imreads done \n');

% Load estimated and refined homographies in previous steps.
% All the refined homographies were saved as mat files.
H12 = load('L1-2'); H12 = H12.H; % Homography of im1 to im2
H23 = load('L2-3'); H23 = H23.H; % Homography of im3 to im2
H34 = load('L3-4'); H34 = H34.H; % Homography of im1 to im2
H54 = load('L5-4'); H54 = H54.H; % Homography of im3 to im2
H65 = load('L6-5'); H65 = H65.H; % Homography of im1 to im2
H76 = load('L7-6'); H76 = H76.H; % Homography of im3 to im2

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