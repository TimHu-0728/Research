function [Bx,By,Bz] = B_multiple(x,y,z,X,Y,Z,Mm,dimensions,r_mag,orient_mag,interp_method)
%B_MULTIPLE 此处显示有关此函数的摘要
%   r_mag, orient_mag are n by 3 array, where n is the number of magnet
%   x, y, z are 3d arrays 
%   Bx, By, Bz are 3d arrays

%%
% Caulculate B_base
B_base = B_1magnet(x,y,z,Mm,dimensions);

%
nX = size(X,1);
nY = size(X,2);
nZ = size(X,3);
n_mag = size(r_mag,1);

% Computations
Bx_multi = zeros(nX,nY,nZ,n_mag);
By_multi = zeros(nX,nY,nZ,n_mag);
Bz_multi = zeros(nX,nY,nZ,n_mag);

% Calculate the magnetic field components Bx, By, Bz of each coil with shifted coil origin
parfor i=1:n_mag
    [Bx_multi(:,:,:,i),By_multi(:,:,:,i),Bz_multi(:,:,:,i)]= B_rotate(x,y,z,X-r_mag(i,1),Y-r_mag(i,2),Z-r_mag(i,3),B_base(:,:,:,1),B_base(:,:,:,2),B_base(:,:,:,3),orient_mag(i,:),interp_method);
end

% Get total magnetic field
Bx=sum(Bx_multi,4);
By=sum(By_multi,4);
Bz=sum(Bz_multi,4);
end

