clc
clear
close all
addpath('functions/')
addpath("Data/")
%% Import Demagnetization curve data of N52 magnets
% K&J Magnetics Demagnetization curve data for N52 magents in CGS unit
% URL: https://www.kjmagnetics.com/bhcurves.asp

% Transfer to SI unit
% Temperature: 20°C
H20  = [-11.30 -11.08 -10.98 -10.84 -10.57 -10.18 0.00 -10.03 -10.49 -10.81 -10.97 -11.08 -11.10]*1e3*1e3/4/pi;    % [A/m]
B20  = [0.00 13.02 13.45 13.74 13.89 14.00 14.38 3.99 3.41 2.93 2.47 1.93 0.00]*1e3*1e-4;                          % [T]

% Temperature: 40°C
H40 = [-9.33 -9.11 -8.96 -8.73 -8.29 -7.65 0.00 -8.04 -8.52 -8.84 -9.01 -9.13 -9.19]*1e3*1e3/4/pi;
B40 = [0.00 12.84 13.29 13.60 13.77 13.89 14.18 5.83 5.20 4.63 4.00 3.23 0.00]*1e3*1e-4;

% Temperature: 60°C
H60 = [-7.53 -7.32 -7.17 -6.93 -6.48 -5.81 0.00 -6.28 -6.74 -7.06 -7.22 -7.33 -7.42]*1e3*1e3/4/pi;
B60 = [0.00 12.63 13.10 13.42 13.60 13.72 13.94 7.41 6.79 6.22 5.57 4.77 0.00]*1e3*1e-4;

% Temperature: 80°C
H80 = [-5.90 -5.70 -5.60 -5.44 -5.13 -4.68 0.00 -4.75 -5.17 -5.45 -5.60 -5.70 -5.80]*1e3*1e3/4/pi;
B80 = [0.00 12.40 12.88 13.21 13.39 13.51 13.67 8.73 8.19 7.70 7.18 6.54 0.00]*1e3*1e-4;

% Temperature: 100°C
H100 = [-4.62 -4.42 -4.32 -4.17 -3.87 -3.45 0.00 -3.50 -3.90 -4.18 -4.32 -4.42 -4.53]*1e3*1e3/4/pi;
B100 = [0.00 12.02 12.55 12.91 13.11 13.24 13.36 9.73 9.20 8.73 8.23 7.62 0.00]*1e3*1e-4;

% Temperature: 120°C
H120 = [-3.50 -3.31 -3.21 -3.05 -2.76 -2.33 0.00 -2.45 -2.83 -3.08 -3.22 -3.31 -3.43]*1e3*1e3/4/pi;
B120 = [0.00 11.53 12.14 12.56 12.79 12.94 13.03 10.49 9.99 9.55 9.08 8.50 0.00]*1e3*1e-4;

% Temperature: 140°C
H140 = [-2.55 -2.38 -2.27 -2.10 -1.79 -1.33 0.00 -1.60 -1.94 -2.17 -2.30 -2.38 -2.50]*1e3*1e3/4/pi;
B140 = [0.00 10.93 11.67 12.17 12.44 12.61 12.67 11.01 10.56 10.16 9.72 9.17 0.00]*1e3*1e-4;

data = [H20;B20;H40;B40;H60;B60;H80;B80;H100;B100;H120;B120;H140;B140];
data_intrinsic = data(:,1:7); 

mu0   = 4 * pi * 1e-7;


H_20 = data_intrinsic(1,:);                                              % [A/m]
M_20 = data_intrinsic(2,:)/mu0;                                          % [A/m]
H_40 = data_intrinsic(3,:); 
M_40 = data_intrinsic(4,:)/mu0;
H_60 = data_intrinsic(5,:); 
M_60 = data_intrinsic(6,:)/mu0;
H_80 = data_intrinsic(7,:); 
M_80 = data_intrinsic(8,:)/mu0;
H_100 = data_intrinsic(9,:); 
M_100 = data_intrinsic(10,:)/mu0;
H_120 = data_intrinsic(11,:); 
M_120 = data_intrinsic(12,:)/mu0;
H_140 = data_intrinsic(13,:); 
M_140 = data_intrinsic(14,:)/mu0;


% Change this value for percentage demagnetization can be tolerant [User Define]
demag_percent = 0.96;

% Create fit curves for H vs M at different temperature using interpolation 
[HM_20,HM_40,HM_60,HM_80,HM_100,HM_120,HM_140,~,~,~,~,~,~,~] = createFit_HM(M_20, H_20, M_40, H_40, M_60, H_60, M_80, H_80, M_100, H_100, M_120, H_120, M_140, H_140);

m20 = linspace(0,M_20(end),3000)';
h20 = HM_20(m20);
m40 = linspace(0,M_40(end),3000)';
h40 = HM_40(m40);
m60 = linspace(0,M_60(end),3000)';
h60 = HM_60(m60);
m80 = linspace(0,M_80(end),3000)';
h80 = HM_80(m80);
m100 = linspace(0,M_100(end),3000)';
h100 = HM_100(m100);
m120 = linspace(0,M_120(end),3000)';
h120 = HM_120(m120);
m140 = linspace(0,M_140(end),3000)';
h140 = HM_140(m140);

figure
plot(m20,h20,m40,h40,m60,h60,m80,h80,m100,h100,m120,h120,m140,h140)
hold on



% Extracting the H field that demagnetize the magnets to 96 percent of Mr
H_demag = [HM_20(demag_percent*M_20(end));
           HM_40(demag_percent*M_40(end));
           HM_60(demag_percent*M_60(end));
           HM_80(demag_percent*M_80(end));
           HM_100(demag_percent*M_100(end));
           HM_120(demag_percent*M_120(end));
           HM_140(demag_percent*M_140(end))];

M_demag = [demag_percent*M_20(end);
           demag_percent*M_40(end);
           demag_percent*M_60(end);
           demag_percent*M_80(end);
           demag_percent*M_100(end);
           demag_percent*M_120(end);
           demag_percent*M_140(end)];

scatter(M_demag,H_demag, 60)

% Temperatures for interpolation
Temp       = [20;40;60;80;100;120;140];

% Create fit curve for T vs H_demag using interpolation
[TH, ~] = createFit_TH(H_demag, Temp);

h_demag = linspace(H_demag(1),H_demag(end),3000)';
T       = TH(h_demag);
H_max   = -0.762425861295835/mu0;
T_max   = TH(H_max);

figure  
plot(h_demag,T)
hold on
scatter(H_demag,Temp)
scatter(H_max,T_max)


% Write data to a folder
mkdir Data

writematrix([m20 h20],'Data/HM_20C.csv')
writematrix([m40 h40],'Data/HM_40C.csv')
writematrix([m60 h60],'Data/HM_60C.csv')
writematrix([m80 h80],'Data/HM_80C.csv')
writematrix([m100 h100],'Data/HM_100C.csv')
writematrix([m120 h120],'Data/HM_120C.csv')
writematrix([m140 h140],'Data/HM_140C.csv')
writematrix([M_demag H_demag],'Data/96% Demagnetization.csv')
writematrix([h_demag,T],'Data/TH line.csv')
writematrix([H_demag Temp],'Data/TH data.csv')
writematrix([H_max T_max],'Data/TH max data.csv')
