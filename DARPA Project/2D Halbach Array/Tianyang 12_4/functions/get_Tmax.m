function T_max = get_Tmax(B_max,data,demag_percent)
%GET_TMAX Calculate the maximum Temperature T_max (°C) such that the magnets, 
%         surrounded by B_max, will not demagnetize to 96% of its original 
%         value Mr (remanence magnetization) based on the demagnetization 
%         curve of N52 Neodymuim magnets. 
% ---------------------------------------------------------------------------------
%         Input : B_max      [T]
%
% ---------------------------------------------------------------------------------
%         Output : T_max      [°C]

%% 
mu0 = 4 * pi * 1e-7;

% K&J Magnetics Demagnetization curve data for N52 magents in CGS unit
% URL: https://www.kjmagnetics.com/bhcurves.asp
H20 = data(1,:);                                              % [A/m]
M20 = data(2,:)/mu0;                                          % [A/m]
H40 = data(3,:); 
M40 = data(4,:)/mu0;
H60 = data(5,:); 
M60 = data(6,:)/mu0;
H80 = data(7,:); 
M80 = data(8,:)/mu0;
H100 = data(9,:); 
M100 = data(10,:)/mu0;
H120 = data(11,:); 
M120 = data(12,:)/mu0;
H140 = data(13,:); 
M140 = data(14,:)/mu0;

% Change this value for percentage demagnetization can be tolerant [User Define]
% demag_percent = 0.96;

% Create fit curves for H vs M at different temperature using interpolation 
[HM_20,HM_40,HM_60,HM_80,HM_100,HM_120,HM_140,~,~,~,~,~,~,~] = createFit_HM(M20, H20, M40, H40, M60, H60, M80, H80, M100, H100, M120, H120, M140, H140);

% Extracting the H field that demagnetize the magnets to 98 percent of Mr
H_demag = [HM_20(demag_percent*M20(end));
           HM_40(demag_percent*M40(end));
           HM_60(demag_percent*M60(end));
           HM_80(demag_percent*M80(end));
           HM_100(demag_percent*M100(end));
           HM_120(demag_percent*M120(end));
           HM_140(demag_percent*M140(end))];

M_demag = [demag_percent*M20(end);
            demag_percent*M40(end);
            demag_percent*M60(end);
            demag_percent*M80(end);
            demag_percent*M100(end);
            demag_percent*M120(end);
            demag_percent*M140(end)];

% Plot the 96% demagnetization points
scatter(M_demag,H_demag,60,'magenta','filled')
legend('', '20 °C', '', '40 °C','', '60 °C','', '80 °C','', '100 °C','', '120 °C','', '140 °C','96% Demagnetization','Location', 'Eastoutside', 'Interpreter', 'none',fontsize = 16  );
hold off

% Temperatures for interpolation
T       = [20;40;60;80;100;120;140];

% Create fit curve for T vs H_demag using interpolation
[TH, ~] = createFit_TH(H_demag, T);

% Use fit curve for T vs H_demag to find maximum temperature T in (°C)
H_max = B_max/mu0;
T_max = TH(H_max);

% Plot T_max and B_max on the fit curve
scatter(H_max,T_max,80,'green','filled')
legend( 'Data', 'Fit Curve','T_max', 'Location', 'Eastoutside', 'Interpreter', 'none',fontsize = 16 );
hold off
end

