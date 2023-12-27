function [Bx_new,By_new,Bz_new] = B_rotate(x,y,z,X,Y,Z,Bx,By,Bz,orient_mag,interp_method)
%B_MULTIPLE 此处显示有关此函数的摘要
%   此处显示详细说明
%%
R_rev  = transpose(eul2rotm(orient_mag, 'ZYX'));
R_forw = eul2rotm(orient_mag, 'ZYX');

% Calculate the position vector before the rotation r_old (m)
XYZ4d  = cat(4,X,Y,Z);
XYZ4d_old  = reshape((R_rev*reshape(XYZ4d,[],3)')',size(XYZ4d));

X_old  = XYZ4d_old(:,:,:,1);
Y_old  = XYZ4d_old(:,:,:,2);
Z_old  = XYZ4d_old(:,:,:,3);

Bx_old = interp3(x,y,z,Bx,X_old,Y_old,Z_old,interp_method);
By_old = interp3(x,y,z,By,X_old,Y_old,Z_old,interp_method);
Bz_old = interp3(x,y,z,Bz,X_old,Y_old,Z_old,interp_method);

B4d_old = cat(4,Bx_old,By_old,Bz_old);
B4d_new = reshape((R_forw*reshape(B4d_old,[],3)')',size(B4d_old));

Bx_new  = B4d_new(:,:,:,1);
By_new  = B4d_new(:,:,:,2);
Bz_new  = B4d_new(:,:,:,3);
end