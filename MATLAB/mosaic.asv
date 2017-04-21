function [ IMG ] = mosaic( IM, im , H, x_min, y_min )
% This function maps image 'im' in a plane to a new image
% 'img' in the mosaiced plane transformed by a homography H
% Arguments:
% im - single image
% H - homography
% IM - current mosaiced image
% x_min - shift to x-direction
% y_min - shift to y-direction
% IMG - updated mosaiced image

% Size of single image
[M,N,C] = size(im);
% Size of mosaiced image
[height,width,C] = size(IM);
IMG = IM;
% Assign pixel values
x_new = [0 0 1]'; % Homogeneous coordinate in new plane
for m = 1 : height
    x_new(2) = m + y_min - 1;
    for n = 1 : width
        x_new(1) = n + x_min - 1;
        for c = 1 : C
            x_org = H \ x_new;
            x = x_org(1) / x_org(3);
            fx = x - fix(x);
            y = x_org(2) / x_org(3);
            fy = y - fix(y);
            if (1 <= x && x <= N && 1 <= y && y <= M)
                % Use bilinear interpolation
                IMG(m,n,c) = (1 - fx) * (1 - fy) * im(fix(y),fix(x),c) +...
                    (1 - fx) * fy * im(ceil(y),fix(x),c) +...
                    fx * (1 - fy) * im(fix(y),ceil(x),c) +...
                    fx * fy * im(ceil(y),ceil(x),c);
            end
        end
    end
end
IMG = uint8(IMG);
end