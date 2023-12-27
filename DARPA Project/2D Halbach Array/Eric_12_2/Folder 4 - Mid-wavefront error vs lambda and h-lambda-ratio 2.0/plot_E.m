clear
clc
clearvars
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

%% Generating variables and values

N         = 201;                                     %defining number of coordinates on each axis
lam_1D    = 0.0254 * 4 * 2 .^ linspace(-8.1, 4.1, N);%1D lambda
hr_1D     = linspace(0, 1, N);                       %1D h/lambda
[lam, hr] = meshgrid(lam_1D, hr_1D);                 %2D input variables

error     = lam * 0;                                 %Because I do a Newton method to find Ms, I have to use for loops.
Ms        = lam * 0;                                 %lam*0 is the easiest way I could think of to get a zero matrix of the correct size
xi        = Ms;
keps      = Ms;
for i = 1:N
    for j = 1:N
        [error(i, j), Ms(i, j), xi(i, j), keps(i, j)] = find_E(lam(i, j), hr(i, j),RhovsMs_cfit);%finding error and optimum Ms
    end
end

%% Generating plots

figure(1)
hold on
megacontour(log10(lam_1D), hr_1D, log10(Ms), 3.9, 3.88, '#1b9e77');
megacontour(log10(lam_1D), hr_1D, xi, 0.1, 0.085, '#7f7f7f');
legend('Rosenswieg instability constraint', '', 'ξ constraint', '')
clines = reshape([1; 3] * 10.^ linspace(-10, 4, 15), [1, 30]);%contour lines

[c11, c12] = contour(log10(lam), hr, error, clines,'ShowText','on','Color', '#7570b3');
c12.DisplayName = 'Mid-wavefront error (m)';%adding this contour to the legend
c12.LineWidth = 1.5;
clabel(c11, c12, 'Color', [0.459 0.439 0.702]);

[c21, c22] = contour(log10(lam), hr, 1./(lam / 4 / 0.0254), 2 .^ linspace(-10,10,21),'ShowText','on','Color', 'k');
c22.DisplayName = 'Magnets per inch';%adding this contour to the legend
c22.LineWidth = 1.5;

[c31, c32] = contour(log10(lam), hr, Ms, clines,'ShowText','on','Color', '#1b9e77');
c32.DisplayName = 'Saturation Magnetization (A/m)';%adding this contour to the legend
c32.LineWidth = 1.5;
clabel(c31, c32, 'Color', [0.106 0.620 0.467]);

[c41, c42] = contour(log10(lam), hr, xi, 10 .^ [-3 -2 -1 0 2 4 6 8],'ShowText','on','Color', '#d95f02');
c42.DisplayName = 'ξ';%adding this contour to the legend
c42.LineWidth = 1.5;
clabel(c41, c42, 'Color', [0.851 0.373 0.008]); 
% clabel(c41, c42, 'Color', [0.906 0.161 0.541]);%'#e7298a'

% [c51, c52] = contour(log10(lam), hr, log10(keps),'ShowText','on','Color', '#66a61e');
% c52.DisplayName = 'log_{10}(kε)';%adding this contour to the legend
% c52.LineWidth = 1.5;
% clabel(c51, c52, 'Color', [0.400 0.651 0.118]);

xlabel("log_{10}(Magnet wavelength λ (m))")
ylabel('Ferrofluid height h/λ (unitless)')
title('log_{10}(Mid-wavefront error)')
hold off

%% Computing specific values

'1/8 in magnet, h/λ = 1/8'
[gm, em] = find_E(0.0254 / 8 * 4, 1/8,RhovsMs_cfit);
gm * 1e9
em

'1/16 in magnet, h/λ = 1/8'
[gm, em] = find_E(0.0254 / 16 * 4, 1/8,RhovsMs_cfit);
gm * 1e9
em

'1/32 in magnet, h/λ = 1/8'
[gm, em] = find_E(0.0254 / 32 * 4, 1/8,RhovsMs_cfit);
gm * 1e9
em