function B = B_1magnet(X,Y,Z,Mm,dimensions)
% X, Y, Z are 3d arrays
% B is 4d matrix, 4th dimension is Bx, By, Bz [T]

%%
% extract magnet dimensions
a      = dimensions(1);
b      = dimensions(2);
h      = dimensions(3);

% generate meshpoints inside magnets
gridsz = 22;                            % Increase this to get more accurate result, but more computation time
x_mag  = linspace(-a/2,a/2,gridsz);
y_mag  = linspace(-b/2,b/2,gridsz);
z_mag  = linspace(-h/2,h/2,gridsz);
x_mag  = x_mag(2:end-1);
y_mag  = y_mag(2:end-1);
z_mag  = z_mag(2:end-1);
[X_mag,Y_mag,Z_mag] = meshgrid(x_mag,y_mag,z_mag);

% Number of inspecting point (assume inspecting point has same number in all three dimensions)
xn = size(X,2);
yn = size(Y,1);
zn = size(Z,3);
A  = zeros([size(X) 3]);
parfor i = 1:xn
    for j = 1:yn
        for k = 1:zn
           A(i,j,k,:) = A_1magnet(X(i,j,k),Y(i,j,k),Z(i,j,k),Mm,X_mag,Y_mag,Z_mag,x_mag,y_mag,z_mag);
        end
    end
end
% figure
% quiver3(X,Y,Z,A(:,:,:,1),A(:,:,:,2),A(:,:,:,3))

% B = curl (A)
[Bx,By,Bz,~] = curl(X,Y,Z,A(:,:,:,1),A(:,:,:,2),A(:,:,:,3));
% figure
% quiver3(X,Y,Z,Bx,By,Bz)
B            = cat(4,Bx,By,Bz);
end