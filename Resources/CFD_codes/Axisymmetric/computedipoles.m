function [m1, z1, m2, z2] = computedipoles(M, Mz, I, Mm, r15, r17, z14, z16)
% Computes the magnetic dipole moments and positions for the magnetic
% system composed of (1) a coil/magnet, and (2) the ferrofluid volume
%
% INPUT: 
% M        = Integral of the magnetization in the ferrofluid volume (Am2)
% Mz       = Integral of magnetization * z in the ferrofluid volume (Am3)
% I        = Current intensity flowing through the coil (0 if magnet) (A)
% Mm       = Magnetization of the magnet (A/m)
% r15      = Radius of the internal wall of the coil/magnet (m)
% r17      = Radius of the external wall of the coil/magnet (m)
% z14      = Height of the lower boundary of the coil/magnet (m)
% z16      = Height of the upper boundary of the coil/magnet (m)
%
% OUTPUT: 
% m1     : Magnetic dipole moment 1 (Am2)
% m2     : Magnetic dipole moment 2 (Am2)
% z1     : Magnetic dipole height 1 (m)
% z2     : Magnetic dipole height 2 (m)

% Preliminary computations
R    = mean([r15, r17]);
Z    = mean([z14, z16]);

% Dipole of coil/magnet
if I ~= 0 % Coil
    m2 = I*pi*R^2;
    z2 = Z;
else % Magnet magnetized vertically at Mm A/m
    m2 = 2*pi*R * (r17-r15) * (z16-z14) * Mm;
    z2 = Z;
end

% Dipole of ferrofluid volume
z1   = Mz/(M+1e-13); 
m1   = M; 