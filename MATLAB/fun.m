function [ F ] = fun( h,x )
% This is a user-defined function to express the homography mapping
% and 'fun' will be used as a parameter of the non-linear least square
% fitting function lsqcurvefit().
% Arguments :
% h = [h11,h12,h13,h21,h22,h23,h31,h32,h33]'
% F = Transformed coordinate by h in the range plane
% x = Coordinate in the domain plane
L = length(x); % Length of xdata
F = zeros(L,1); % Initialize function values
for k = 1:2:L
    F(k) = (h(1)*x(k)+h(2)*x(k+1)+h(3))/(h(7)*x(k)+h(8)*x(k+1)+h(9));
    F(k+1) = (h(4)*x(k)+h(5)*x(k+1)+h(6))/(h(7)*x(k)+h(8)*x(k+1)+h(9));
end
end