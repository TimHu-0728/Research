function A = A_1magnet(X,Y,Z,Mm,X_mag,Y_mag,Z_mag,x_mag,y_mag,z_mag)
%A_1MAGNET Simulate A field of 1 magent, upward , at the origin
%   [x, y, z] is vector, A = [Ax, Ay, Az] is a vector, Mm [A/m]
%   dimensions = [a,b,h] [m]

mu0    = 4*pi*1e-7;
% position vector from source points to inspecting point
rx     = -X_mag+X;
ry     = -Y_mag+Y;
rz     = -Z_mag+Z;
r_norm = sqrt(rx.^2+ry.^2+rz.^2);
r_hat_over_r_square  = cat(4,rx./(r_norm.^3),ry./(r_norm.^3),rz./(r_norm.^3));

% Magnetic dipole moment per unit volume
Mx     = 0*ones(size(rx));
My     = 0*ones(size(ry));
Mz     = Mm*ones(size(rz));
M      = cat(4,Mx,My,Mz);

% Calculate integrand
integrand = mu0/(4*pi)*cross(M,r_hat_over_r_square,4);

% Numerical integration to get A
A = squeeze(trapz(z_mag,trapz(y_mag,trapz(x_mag,integrand,2),1),3));
end

