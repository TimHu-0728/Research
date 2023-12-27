function plot_magnet(r_mag,orient_mag,dimensions,plot_mag)
%PLOT_MAGNET 此处显示有关此函数的摘要
%   此处显示详细说明

% Create a rotation matrix from the Euler angles
R = eul2rotm(orient_mag, 'ZYX');
direction = R * [0;0;1];
% Define the vertices of the magnet at the origin
a = dimensions(1);
b = dimensions(2);
h = dimensions(3);
vertices = [a/2 b/2 h/2;
            -a/2 b/2 h/2;
            -a/2 -b/2 h/2;
            a/2 -b/2 h/2;
            a/2 b/2 -h/2;
            -a/2 b/2 -h/2;
            -a/2 -b/2 -h/2;
            a/2 -b/2 -h/2];

%Apply rotation and translation to the vertices
vertices = (R * vertices' + r_mag').';

% Define the faces of the cuboid
faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];
c     = [112,128,144];
c     = c./norm(c);
%Plot the cuboid
hold on;
if plot_mag
    for i = 1:size(faces, 1)
        X = vertices(faces(i,:), 1);
        Y = vertices(faces(i,:), 2);
        Z = vertices(faces(i,:), 3);
        fill3(X, Y, Z, c); % 'b' sets the color blue for all faces
    end
end
% Plot orientation arrow (e.g., along the local z-axis of the cuboid)
quiver3(r_mag(1), r_mag(2), r_mag(3), direction(1), direction(2), direction(3), 0.8*h, 'k', 'LineWidth', 1.5,'MaxHeadSize',1.5);
end



