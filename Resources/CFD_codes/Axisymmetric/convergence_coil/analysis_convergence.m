%% Convergence
clear all

%% Parameters
% Physics
N        = 200;
I        = 5;

% Geometry
RminD    = 0.0805;  % Minimum coil radius
RmaxD    = 0.1115;  % Maximum coil radius
ZminD    = -0.0075; % Minimum coil height
ZmaxD    = +0.0205; % Maximum coil height

% Domain
[r,z]    = meshgrid(linspace(0,0.6,100),linspace(-0.6,0.6,200));

% Meshing
nr       = 10;
nz       = nr;

%% Computations
figure
for nr = 1:2:20
    nz = nr;
    % Cylindrical thread vector potential
    if nr > 1
        phi0     = phi;
    end
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
    
    %% Represent
    subplot(2,5,(nr-1)/2+1)
    hold on
    if nr == 1
        contourf(r,z,phi,20)
        title('Potential single loop')
    else
        contourf(r,z,log10(abs((phi-phi0)./phi0)*100),20)
        title(['Variation (log10(%)) for n_r=',num2str(nr)])
    end
    colorbar
    plot([RminD,RmaxD,RmaxD,RminD,RminD],[ZminD,ZminD,ZmaxD,ZmaxD,ZminD],'k-','linewidth',1)
end