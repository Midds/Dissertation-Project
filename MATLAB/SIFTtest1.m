close all;
imds = imageDatastore({'im15.jpeg';'im16.jpeg';'im17.jpeg';});

figure;
montage(imds.Files);
title('Montage');

