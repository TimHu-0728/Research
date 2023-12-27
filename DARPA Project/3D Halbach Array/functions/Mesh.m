function [X,Y,Z] = Mesh(xi,xf,yi,yf,zi,zf,gridsz)
%MESH Summary of this function goes here
%   Detailed explanation goes here
x       = linspace(xi,xf,gridsz);
y       = linspace(yi,yf,gridsz);
z       = linspace(zi,zf,gridsz);
[X,Y,Z] = meshgrid(x,y,z);
end

