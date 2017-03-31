close all;
%% File for cropping barratxx1 and barratxx2 videos
I = imread('barret1/im1.jpeg');

xmin = 100;
ymin = 35;
width = 220;
height = 170;
rect = [xmin ymin width height];

I2 = imcrop(I,rect);

imshow(I2);

