% Spherical Mirror Sensitivity Analysis
%
% 
% Citation: Cannavo' F., Sensitivity analysis for volcanic source modeling quality assessment and model selection, Computers & Geosciences, Vol. 44, July 2012, Pages 52-59, ISSN 0098-3004, http://dx.doi.org/10.1016/j.cageo.2012.03.008.

clear;
clc
close all
addpath('./GSAT')

tic
% create a new project 
pro = pro_Create();


% add to the project 10 input variables, named X*, distributed in the 
% range [ ] and indicate that the variables will be sampled following a 
% Sobol quasi-random set 
pro = pro_AddInput(pro, @()pdf_Sobol([800 2000]), 'X1');              % Ferrofluid material density (Ferrotec Inc.) [kg/m^3]
pro = pro_AddInput(pro, @()pdf_Sobol([0 100000]), 'X2');             % Saturation magnetization of Ferrofluid mateirial (Ferrotec Inc.) [A/m]
pro = pro_AddInput(pro, @()pdf_Sobol([0 13]), 'X3');                % Initial Magenetic susceptibility of Ferrofluid (Ferrotec Inc.)
pro = pro_AddInput(pro, @()pdf_Sobol([sqrt(9.806/30) 2*sqrt(9.806/30)]), 'X4');  % Angular velocity w [rad/s]
pro = pro_AddInput(pro, @()pdf_Sobol([0 exp(30)]), 'X5');           % Current inside coil 1 [A]
pro = pro_AddInput(pro, @()pdf_Sobol([0 exp(30)]), 'X6');           % Current inside coil 2 [A]
pro = pro_AddInput(pro, @()pdf_Sobol([0 exp(30)]), 'X7');           % Current inside coil 3 [A]
pro = pro_AddInput(pro, @()pdf_Sobol([0 exp(30)]), 'X8');           % Current inside coil 4 [A]
pro = pro_AddInput(pro, @()pdf_Sobol([0 20]), 'X9');                % R Positions of coil 1 [m]
pro = pro_AddInput(pro, @()pdf_Sobol([-6 4.5]), 'X10');              % Z Positions of coil 1 [m]
pro = pro_AddInput(pro, @()pdf_Sobol([0 20]), 'X11');               % R Positions of coil 2 [m]
pro = pro_AddInput(pro, @()pdf_Sobol([-6 4.5]), 'X12');              % Z Positions of coil 2 [m]
pro = pro_AddInput(pro, @()pdf_Sobol([0 20]), 'X13');               % R Positions of coil 3 [m]
pro = pro_AddInput(pro, @()pdf_Sobol([-6 4.5]), 'X14');              % Z Positions of coil 3 [m]
pro = pro_AddInput(pro, @()pdf_Sobol([0 20]), 'X15');               % R Positions of coil 4 [m]
pro = pro_AddInput(pro, @()pdf_Sobol([-6 4.5]), 'X16');              % Z Positions of coil 4 [m]


% set the model, and name it as 'model', to the project 
pro = pro_SetModel(pro, @(x)mymodel(x), 'model');

% set the number of samples for the quasi-random Monte Carlo simulation
pro.N = 50000;

% initialize the project by calculating the model at the sample points
pro = GSA_Init(pro);


% calculate the first order global sensitivity coefficients by using FAST
% algorithm
Sfast = GSA_FAST_GetSi(pro);
Sfast_norm = Sfast./sum(Sfast);

% Plot
explode = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
labels  = {'rho','Ms','chi0','w','I1','I2','I3','I4','X1','Y1','X2','Y2','X3','Y3','X4','Y4',};
pie(Sfast_norm,explode,labels)
title('Sensitivity coefficients for each variable (First Order)','FontSize',18)
toc