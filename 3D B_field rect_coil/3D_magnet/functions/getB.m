function B=getB(x,y,z,a,b,N,I,h,mu,psi,theta,phi) 
% Function calculate components of Magnetic Field Bx, By, Bx with input parameters:
% Position r = [x,y,z] (m),
% Coil length a (m), width b (m), total height h (m), number of stack N, current I (A), mu (N*A^(-2))
% Coil orientation angles, psi,theta,phi (rad)
%
% B=getB(x,y,z,a,b,N,I,h,mu,psi,theta,phi) 
%
% INPUTS:
%
% OUTPUTS:
%

% Computations
% Calculate the position vector before the rotation r_old (m)
r_old= quatrotate(eul2quat(-[psi theta phi]),[x,y,z]);

% Evaluate the magnetic field components Bx, By, Bz at r_old
B_old=B_coil(r_old(1),r_old(2),r_old(3),a,b,mu,I,N,h);

% Rotate the magnetic field B evaluated at r_old to the new position
B=-quatrotate(eul2quat([psi theta phi]),B_old')';
end


