clear
clc
clear
close all

%% 
% Interpolation data for EMG water-based series ferrofluid (Ferrotec Inc.)
% URL: https://ferrofluid.ferrotec.com/products/ferrofluid-emg/water/
% EMG   304, 308 408 507  508  509 605  607  700   705  707  708
Chi0_data = [5.03 0.5  0.5  1.63 0.88 0.5  3.02 1.63 12.57 4.04 1.51 0.63];     % Initial Susceptibility
Ms_data   = [27.5 6.6  6.6  11   6.6  3.3  22   11   35.5  22   11   6.6] * 795.7747;      % Saturation Magnetization [A/m]
phi_data  = [4.5  1.2  1.2  2    1.2  0.6  3.9  2    5.8   3.9  2    1.2];      % Magnetic Particle concentration [%]
rho_data  = [1240 1060 1070 1120 1070 1030 1180 1100 1290  1190 1100 1080];     % Density @ 25 C [kg/m^3]

% Find fit curves
[RhovsMs_cfit, gof1] = createFit_RhovsMs(Ms_data, rho_data);
[Chi0vsRho_cfit, gof2] = createFit_Chi0vsRho(rho_data, Chi0_data);


%%
stack_plots = 1;%If true, stacks all of the plots into one. If false, breaks them into four.

%% list of our four magnets' widths and heights (b) - modify this to add/remove magnets
width_list  = [1, 4, 8, 1, 4, 8, 1, 4, 8] * 0.0254 / 8;
b_list      = [1, 1, 1, 4, 4, 4, 8, 8, 8] * 0.0254 / 8;

%% Constants
 
g           = 9.81;                           %m/s^2
mu0         = 4*pi*1e-7;                      %m*kg/(s^2 A^2)
sigma       = 0.025;                          %N/m, for water based ferrofluid
M0          = 1.48/mu0;                       %N42 (1.32) ->N52 (1.48) grade magnet, gauss is in kg/(A*s^2), this is A/m
delta       = 0.001;                          %m, thickness of the ferrofluid layer

x           = 0;                              %here x and y represent the simulation point
y           = delta;

N           = 1000;                           %number of points per axis
Ms_1D       = 28500*10.^linspace(-2.5, 0, N); %range of Ms values we want
h_1D        = 0.001*linspace(0,20,N);         %range of h values we want

[Ms, h]     = meshgrid(Ms_1D, h_1D);          %2D grid of Ms and h for use in further calculations

rho         = reshape(RhovsMs_cfit(Ms),size(Ms));               %kg/m^3; EMG 700 base is 1290, and is linearly modified with water
chi0        = reshape(Chi0vsRho_cfit(rho),size(Ms));            %EMG700 in its undiluted state has a magnetic susceptability of 1

%% magnet of choice:

t = cell(length(b_list));
for i=[1:length(b_list)]                      %generating figure titles
    t{i}=string(width_list(i) / 0.0254) + " x " + string(b_list(i) / 0.0254) + " in magnet";
end

for i=[1:length(b_list)]                      %different figures for each magnet
    width       = width_list(i);              %width of magnet
    b           = b_list(i);                  %height of magnet - both this and the width are from the width and height lists

    %% Surface roughness and force
    
    lambda      = 4.*width+0.004;                            %frequency length
    k           = 2.*pi./lambda;                       %wave number, m^-1
    H0          = (2^1.5 / pi) * M0 .* exp(-k.*h) .* (1 - exp(-k.*b));%Calculating magnetic H field at the surface of the ferrofluid
    beta        = (1 - exp(-5 .* k .* b))/(1 - exp(-1 .* k .* b));
    xi          = Ms.* exp(k.*delta)./H0;              %Calculating xi for verification
    f_vol_term1 = -mu0.*k.*Ms.*H0.*exp(-k.*y);
    f_vol_term2 = (1 + xi .* exp(2.*k.*(y-delta)).*cos(2.*k.*x));
    f_vol_y     = f_vol_term1 .* f_vol_term2;          %magnetic force per unit volume
    
    f_vol_y_g   = f_vol_y./(rho.*g);                   %magnetic force per unit mass - in units of g
    
    Roughness   = abs(xi./(4*k) .* (1 - 0.4 .* beta .* exp(-4 .* k .* (h + delta))) ./ (1 + (4.*sigma.*k + rho.*g./k)./(mu0.*Ms.*H0.*exp(-k.*delta)) - beta .* exp(-4 .* k .* (h + delta)))     );%Mid-wavefront error
    fraction    = 1 - (4.*sigma.*k + rho.*g./k)./(mu0.*Ms.*H0.*exp(-k.*delta));%Fraction within surface roughness - this is zero when Roughness is zero, but varies less logatithmically

    %% Rosenswieg instability
    Hx = H0;                                                                                   %worst case scenario for mallenson configuration - H0 is entirely in the normal direction
    gamma   = 3.*chi0./Ms;
    r0      = sqrt(((H0+Ms.*(coth(gamma.*H0)-1./(gamma.*H0))).*(1+Ms.*(1./(H0.^2.*gamma)-gamma.*csch(H0.*gamma).^2)))./H0);%(geometric mean of the chord and tangent permeabilities)
    Mc      = sqrt(2./mu0.*(1+1./r0).*(sqrt(g.*rho.*sigma)+mu0.*(1.-r0).^2/2./(1.+r0).*Hx.^2));% Critical magnetization[A/m]
    M       = Ms.*(coth(gamma.*H0)-1./gamma./H0);                                              % Magnetization of ferrofluid [A/m]
    inst    = M - Mc;                                                                          %Represents the instability - if this is negative, we are good

    allcon  = max((-10 - f_vol_y_g),0) .* max((1 - xi * 10),0) .* max((0 - inst), 0);          %This is all of our constraints except the Roughness, and is >0 when we are feasible
    
    %bpoint1analy = %proposed analytical expression for the minimum Ms

    %% Plotting
    figure(i)                                                                              %generating the figure for this magnet
    hold on
    if stack_plots
        title(t(i))
        bpoint = megacontour(log10(Ms_1D), h_1D * 1000, f_vol_y_g, -10, -10.5,'#1b9e77');        %force constraint and ideal point
        megacontour(log10(Ms_1D), h_1D * 1000, xi, 0.01, 0.0095,'#823901');                     %stringent xi constraint
        megacontour(log10(Ms_1D), h_1D * 1000, xi, 0.1, 0.095,'#d95f02');                       %weak xi constraint
        megacontour(log10(Ms_1D), h_1D * 1000, inst, 0, -250,'#404040');                       %Rosenswieg constraint

        contourf(log10(Ms_1D), h_1D * 1000, allcon ./ allcon + allcon, [0 0], "FaceAlpha",0.25)%paint all feasible points a translucent color
        colormap("summer")                                                                %paint them green
        plot(bpoint(1), bpoint(2), '.k', 'MarkerSize',20)                                      %best point
        text(bpoint(1) - 0.1, bpoint(2) + 0.35,{'Best point'})
        xlabel("log_{10}(Saturation Magnetization (A/m))")
        ylabel('Ferrofluid height (mm)')
        legend('Force > 10g constraint', '', 'ξ < 0.01 constraint', '', 'ξ < 0.1 constraint', '', 'Rosenswieg instability constraint')
        
        clines = reshape([1; 3] * 10.^ linspace(-10, 4, 15), [1, 30]);%contour lines

        [c11, c12] = contour(log10(Ms), h * 1000, Roughness * 1e6, clines,'Color', '#7570b3','ShowText','on');  %plot of log Roughness
        c12.DisplayName = 'Mid-wavefront error (μm)';                                                      %adding this contour to the legend
        c12.LineWidth = 1.5;                                                                              %making the lines visible
        clabel(c11, c12, 'Color', [0.459 0.439 0.702]);                                                   %coloring the labels
        [c21, c22] = contour(log10(Ms), h * 1000, 0-f_vol_y_g, clines,'Color', '#1b9e77','ShowText','on');%plot of log force
        c22.DisplayName = 'Magnetic force (g)';                                                           %adding this contour to the legend
        c22.LineWidth = 1.5;                                                                              %making the lines visible
        clabel(c21, c22, 'Color', [0.106 0.620 0.467]);                                                   %coloring the labels
        [c31, c32] = contour(log10(Ms), h * 1000, xi, clines,'Color', '#d95f02','ShowText','on');         %plot of log xi
        c32.DisplayName = 'ξ';                                                                            %adding this contour to the legend
        c32.LineWidth = 1.5;                                                                              %making the lines visible
        clabel(c31, c32, 'Color', [0.851 0.373 0.008]);                                                   %coloring the labels
    else
        sgtitle(t(i))                                                                          %showing the viewer which magnet this figure is for
    
        ax4 = subplot(2,2,4);                                                                  %plot of all of the constraints - ax4 is for imposing a different colormap
        bpoint = megacontour(log10(Ms_1D), h_1D * 1000, f_vol_y_g, -10, -10.5,'#1b9e77');      %force constraint and ideal point
        megacontour(log10(Ms_1D), h_1D * 1000, xi, 0.01, 0.0095,'#823901');                    %stringent xi constraint
        megacontour(log10(Ms_1D), h_1D * 1000, xi, 0.1, 0.095,'#d95f02');                      %weak xi constraint
        megacontour(log10(Ms_1D), h_1D * 1000, inst, 0, -250,'#404040');                       %Rosenswieg constraint
    
        contourf(log10(Ms_1D), h_1D * 1000, allcon ./ allcon + allcon, [0 0], "FaceAlpha",0.25)%paint all feasible points a translucent color
        colormap(ax4, "summer")                                                                %paint them green
        plot(bpoint(1), bpoint(2), '.k', 'MarkerSize',20)                                      %best point
        text(bpoint(1) - 0.35, bpoint(2) + 1,{'Best point'})
    
        axis([min(log10(Ms_1D)) max(log10(Ms_1D)) min(h_1D * 1000) max(h_1D * 1000)])          %making sure the axes in this subplot matches those of the others
        xlabel("log_{10}(Saturation Magnetization (A/m))")
        ylabel('Ferrofluid height (mm)')
        legend('Force > 10g constraint', '', 'ξ < 0.01 constraint', '', 'ξ < 0.1 constraint', '', 'Rosenswieg instability constraint')
        title('Plot of constraints')
    
        subplot(2,2,1)%plot of log Roughness
        contourf(log10(Ms), h * 1000, log10(Roughness),'ShowText','on');
        colorbar;
        xlabel("log_{10}(Saturation Magnetization (A/m))")
        ylabel('Ferrofluid height (mm)')
        title('log_{10}(Mid-wavefront error (m))')
    
        subplot(2,2,2)%plot of log force
        contourf(log10(Ms), h * 1000, log10(0-f_vol_y_g),'ShowText','on');
        colorbar;
        xlabel("log_{10}(Saturation Magnetization (A/m))")
        ylabel('Ferrofluid height (mm)')
        title('log_{10}(Magnetic force (g))')
    
        subplot(2,2,3)%plot of log xi
        contourf(log10(Ms), h * 1000, log10(xi),'ShowText','on');
        colorbar;
        xlabel("log_{10}(Saturation Magnetization (A/m))")
        ylabel('Ferrofluid height (mm)')
        title('log_{10}(ξ)')
    end
    hold off
end