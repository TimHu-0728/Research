function [f_vol_y_g, Roughness, xi, keps] = find_GandE(lambda, hratio, Ms,RhovsMs_cfit)
    %% list of our four magnets' widths and heights (b) - modify this to add/remove magnets
    b      = lambda;
    
    %% Constants
     
    g           = 9.81;                           %m/s^2
    mu0         = 4*pi*1e-7;                      %m*kg/(s^2 A^2)
    sigma       = 0.025;                          %N/m, for water based ferrofluid
    M0          = 1.48/mu0;                       %N42 (1.32) ->N52 (1.48) grade magnet, gauss is in kg/(A*s^2), this is A/m
    delta       = 0.001;                          %m, thickness of the ferrofluid layer
    
    x           = 0;                              %here x and y represent the simulation point
    y           = delta;
    
    N           = 1000;                           %number of points per axis
    h           = hratio*lambda;                  %range of h values we want
    
    chi0        = Ms / 28500;                     %EMG700 in its undiluted state has a magnetic susceptability of 1
    rho         = RhovsMs_cfit(Ms);      %kg/m^3; EMG 700 base is 1290, and is linearly modified with water
    
    %% magnet of choice:
    
    %% Surface roughness and force
    
    k           = 2.*pi./lambda;                       %wave number, m^-1
    H0          = (2^1.5 / pi) * M0 .* exp(-k.*h) .* (1 - exp(-k.*b));%Calculating magnetic H field at the surface of the ferrofluid
    beta        = (1 - exp(-5 .* k .* b))/(1 - exp(-1 .* k .* b));
    xi          = Ms.* exp(k.*delta)./H0;              %Calculating xi for verification
    f_vol_term1 = -mu0.*k.*Ms.*H0.*exp(-k.*y);
    f_vol_term2 = (1 + xi .* exp(2.*k.*(y-delta)).*cos(2.*k.*x));
    f_vol_y     = f_vol_term1 .* f_vol_term2;          %magnetic force per unit volume
    
    f_vol_y_g   = f_vol_y./(rho.*g);                   %magnetic force per unit mass - in units of g
    
    Roughness   = abs(xi./(4*k) .* (1 - 0.4 .* beta .* exp(-4 .* k .* (h + delta))) ./ (1 + (4.*sigma.*k + rho.*g./k)./(mu0.*Ms.*H0.*exp(-k.*delta)) - beta .* exp(-4 .* k .* (h + delta)))     );%Mid-wavefront error
    keps        = k.* eps;
end