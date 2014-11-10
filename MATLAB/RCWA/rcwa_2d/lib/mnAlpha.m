function [m,n] = mnAlpha(a,maxx)
% function [m,n] = mnAlpha(a,maxx)
% Inverzni prevod k Alpha

mMAX = maxx(1);
nMAX = maxx(2);
n = floor((a-1)/(2*mMAX+1)) - nMAX;
m = rem(a-1, 2*mMAX+1) - mMAX;
