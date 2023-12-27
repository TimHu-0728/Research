%% Convergence
clear all

%% Parameters
% Physics
N        = 200;
I        = 20;

% Geometry
RminD    = 0.0805;  % Minimum coil radius
RmaxD    = 0.1115;  % Maximum coil radius
ZminD    = -0.0075; % Minimum coil height
ZmaxD    = +0.0205; % Maximum coil height

% Domain
[r,z]    = meshgrid(linspace(0,0.15,200),linspace(-0.1,0.15,400));

% Meshing
nr       = 15;
nz       = nr;

%% Computations

% Cylindrical thread vector potential
phi      = zeros(size(r));
for i = 1:nr
    for j = 1:nz
        ra    = RminD + (RmaxD-RminD)/nr * (1/2 + i-1);
        za    = ZminD + (ZmaxD-ZminD)/nz * (1/2 + j-1);
        k     = sqrt(4 * ra * r ./ ((ra+r).^2 + (z-za).^2));
        [K,E] = ellipke(k.^2);
        phi   = phi + N*I*sqrt(ra.*r)./(pi*k*nr*nz) .* ((1-k.^2/2).*K - E);
    end
end

%% Magnetic field
Aphi    = phi ./ r;
[dAphidr,dAphidz] = gradient(Aphi,r(1,2)-r(1,1),z(2,1)-z(1,1));
Hr                = -dAphidz;
Hz                = Aphi./r + dAphidr;
H                 = sqrt(Hr.^2+Hz.^2);

%% Represent
figure
hold on
plot3(r(:),z(:),phi(:),'.')
plot([RminD,RmaxD,RmaxD,RminD,RminD],[ZminD,ZminD,ZmaxD,ZmaxD,ZminD],'k-','linewidth',1)
xlabel('r [m]')
ylabel('z [m]')

figure
hold on
contourf(r,z,H,0:2500:45000)
h = colorbar; ylabel(h,'$H$ [A/m]', 'Interpreter', 'latex','FontSize', 14)
plot([RminD,RmaxD,RmaxD,RminD,RminD],[ZminD,ZminD,ZmaxD,ZmaxD,ZminD],'k-','linewidth',1)
axis equal
xlabel('r [m]')
ylabel('z [m]')
