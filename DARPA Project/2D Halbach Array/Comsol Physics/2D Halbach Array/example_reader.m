%% Example import
% This script shows how to use Matlab's importfile tool. 
%
% Álvaro Romero-Calvo, V1.0 (17/12/2023)

% Parameters
filename  = './B_0.5_0.5_0.5.txt';
side      = 0.5*0.0254;             % Magnet side length (m)

% File loading (B field at y = 0)
[x, ~, z, bx, ~, bz] = importfile(filename);

% Matrix form
xu        = unique(x);
zu        = unique(z);
X         = zeros(length(xu),length(zu));
Z         = zeros(length(xu),length(zu));
Bx        = zeros(length(xu),length(zu));
Bz        = zeros(length(xu),length(zu));
for i = 1:length(xu)
    for j = 1:length(zu)
        X(i,j)  = xu(i);
        Z(i,j)  = zu(j);
        Bx(i,j) = bx(x==xu(i) & z == zu(j));
        Bz(i,j) = bz(x==xu(i) & z == zu(j));
    end
end

% Intermediate variables
B         = sqrt(Bx.^2 + Bz.^2);
% [x,z]     = meshgrid(0,0.00635);
% B_top     = interp2(X,Z,B,x,z,'spline')
% B_side    = interp2(X,Z,B,x,z,"spline")
%% Plot
figure, hold on
contourf(X,Z,B)
% quiver(X,Z,Bx./B,Bz./B, 'r')
patch([-side, side, side, -side]/2, [-side, -side, side, side]/2,'k')
xlabel('x (m)')
ylabel('z (m)')
c = colorbar; 
ylabel(c,'B (T)')
axis equal
