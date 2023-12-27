%% Studies magnetically induced viscosity

% Parameters
mu0 = 4*pi*1e-7;
Hn  = 0:100:1e6; %A/m
eta = 1.445/1000; %Pa*s
phi = 0.58 / 100; % Volume fraction of solids
beta = pi/2;
k   = 1.38*1e-23; %Nm/K
T   = 298; %K
d   = 10*1e-9; % Particle diameter (m)
V   = 4/3 * pi*(d/2)^3; % Particle volume(m-3)

% Magnetization
aM    = 459.70*2/pi; % Parameters of the magnetization curve
bM    = 2747.15*2/pi;
cM    = 5.7304*10^(-6);
dM    = 1.0267*10^(-4);
eM    = 0;
M     = (aM*atan(cM*Hn)+bM*atan(dM*Hn)+eM*Hn);

% zeta and tau
zeta = 1.5 * eta * phi;
tau  = 3*V*eta/k/T;

% Viscosity
Deta = zeta * (mu0 * M .* Hn * tau)./(4*zeta + mu0*M.*Hn*tau) * sin(beta)^2;

% Represent
figure
plot(Hn, eta+Deta)