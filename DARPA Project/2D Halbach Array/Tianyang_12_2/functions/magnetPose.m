function [r_mag,orient_mag] = magnetPose(n_mag,range)
%MAGNETPOSE 此处显示有关此函数的摘要
%   此处显示详细说明

x = linspace(-range,range,n_mag);
y = linspace(-range,range,n_mag);
[X,Y] = meshgrid(x,y);
Z = -5*ones(n_mag*n_mag,1);

r_mag      = zeros(n_mag*n_mag,3);
orient_mag = zeros(n_mag*n_mag,3);
idx        = 1;
for i = 1:(n_mag*n_mag)
    r_mag(i,:)      = [X(i) Y(i) Z(i)];
    orient_mag(i,:) = deg2rad([0 0 90*(idx-1)]);
    if idx == n_mag
        idx = idx - n_mag;
    end
    idx = idx + 1;
end
end

