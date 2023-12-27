%% Solver de Miguel Angel Herrada
% V3 (14/03/2021)
clear all

%% Initial setup
% Lets clean before start and decide where to locate jacobians functions
% and auxiliary elements

% Path locations
path_jacobian  = './eigen/jacobians2D/';
path_jacobian1 = './subroutinesMatlab/';
path_mesher    = './meshgenerator/';
path_jacmesh   = './meshgenerator/jacobians2D';

% Path operations
restoredefaultpath
addpath(path_jacobian)
addpath(path_jacobian1)
addpath(path_mesher)
addpath(path_jacmesh)

%% Parameters
% Geometry
% Block A (Ferrofluid)
R     = 0.055; % Radius of the container
H     = 0.05;  % Initial height of the liquid
nzA   = 101;    % Points in r
nrA   = 101;    % Points in z

% Block B (Air over ferrofluid)
H1    = 0.1; % Height of the dry container
nzB   = nzA; % Points in the r
nrB   = 41;  % Points in the z

% Block C (Air over container)
H2    = 0.15; % Height of the domain
nzC   = nzA; % Points in z
nrC   = 21;  % Points in r

% Block E (Air under container)
zout  = -0.1;  % Lower end of the domain
nzE   = nzA; % Points in z
nrE   = 61;  % Points in r

% Block D+ (Coil)
RminD    = 0.0805;  % Minimum coil radius
RmaxD    = 0.1115;  % Maximum coil radius
ZminD    = -0.0075; % Minimum coil height
ZmaxD    = +0.0205; % Maximum coil height
Rout     = 0.15;   % External radius of the domain
d_nz     = 14;  % Density of points per cm of D domain
nz1      = 3*floor((100*d_nz * (RminD-R))/2)+1;     % Pointers before the coil (in the radial direction)
nzbobina = floor((100*d_nz * (RmaxD-RminD))/2)+1; % Pointers inside the coil
nz2      = floor((100*d_nz * (Rout-RmaxD))/2)+1;  % Pointers after the coil

% Physical properties
rho_fe   = 1012;        % Density of the ferrofluid (kg/m3)
eta_fe   = 1.445/1000;  % Dynamic viscosity (Ns/m2 or Pa·s)
gamma_fe = 0.06239;     % Surface tension (N/m)
theta_fe = 113*pi/180;  % contact angle (rad)
Ib_fe    = 0;           % Coils current intensity (A)
Nturns   = 200;         % Number of turns of the coils
Mmagnet_fe = 0;         % Magnetization vector of magnet (A/m)
aM_fe    = 459.70*2/pi; % Parameters of the magnetization curve
bM_fe    = 2747.15*2/pi;
cM_fe    = 5.7304*10^(-6);
dM_fe    = 1.0267*10^(-4);
eM_fe    = 0;
zeta_fe  = 1.2571e-05; % Coefficients of magnetic viscosity
tau_fe   = 5.5194e-07;
grav_fe  = 0;

% Time simulation parameters
dt      = 100;
Nsteps  = 1;
miniter = 2;
maxerror= 0.001;
warning off

% Block compilation variables
leq      = 0; % 1 to symbolic functions and construction of the matrix, otherwise
lb       = 1; % Specifies block to be compiled (1:A, 2:B, 3:C, 4:D, 5:E)
lstart   = 0; % 0 to start meshing from scratch, 1 to start with name_prev_sol solution
lcompute = 1; % Compute mesh
name_prev_sol = './soluciones/JAM-09-May-2021_Ib0_zeta1,2571e-05tau5,5194e-07nrA101nzA101nrB41nrC21nrE61nzE101_t0.27.mat';
filesave = 0; % Save (1) or don't save (0) solution

%% Details of the block, variables and connections
%
%
% Fix the number of blocks (subdomains). In each block it has to be defined
% the number of variables and the number of symbolic variables.
%
%%
% number of parameters

Np = 31;

% To construct the numerical scheme we need to use the list of variables of
% each block (*Important: In the same order as defined in the corresponding
% block*).

list_block = {'A' 'B' 'C' 'D' 'E'}; %see scheme

list_var_A = {'wA', 'uA', 'pA','phiA','fA' 'gA' 'CA' 'etaA' 'psiA', 'VA' 'gaxisA'}; %phiA is the streamfunction for H CA is magnetic susceptibility
list_var_B = {'phiB','fB' 'gB' 'psiB'}; %phiB is the streamfunction for H
list_var_C = {'phiC','fC','psiC'}; %phiC is the streamfunction for H
list_var_E = {'phiE' ,'fE','psiE' }; %magentic streamfuction CE is magnetic susceptibility

list_var_D = {'phiD','psiD' 'FD' 'GD' 'MzD'};%phiB is the streamfunction for H ellipitical mesh



%conextions between blocks
list_con_A = {'AB' 'AD'};
list_con_B = {'BA' 'BC''BD'};
list_con_C = {'CB' 'CD'};
list_con_D = {'DA' 'DB' 'DC' 'DE'};
list_con_E = {'ED' 'EA' };



[basura nbl ] = size(list_block);


% The list of derivatives. We have assumend that the number and distribution
% is the same in all the block. This information is relevant in _Matrix_
%
list_der_A = {'', '0r', '0z', '0rr', '0zz', '0rz', '0time'};
list_der_B = {'', '0r', '0z', '0rr', '0zz', '0rz', '0time'};
list_der_C = {'', '0r', '0z', '0rr', '0zz', '0rz', '0time'};
list_der_D = {'', 'r0', 'z0', 'rr0', 'zz0', 'rz0' ,'t0'}; %mapping derivatives asumming y=FA(r0,z0,t0,time), x=GA(r0,z0,t0,time), z=t0
list_der_E = {'', '0r', '0z', '0rr', '0zz', '0rz', '0time'};
list_der=list_der_A;
for i=1:nbl
    bl = list_block{i};
    order = ['[bas NV' bl '] = size(list_var_' bl ');'];
    evalc(order);
    %     order = ['[bas ND' bl '] = size(list_der);' ];
    %     evalc(order);
    
    order = ['[bas ND' bl '] = size(list_der_' bl ');' ];
    evalc(order);
end


%% Compute symbolic functions and construction of matrix.m
if (leq==1)
    
    for i=1:Np
        pa(i,1)=sym(['pa_',num2str(i)],'real');
    end
    
    %% Form the array  of dimensionless parameters
    %phisiclal properties of the liquid
    rho   = pa(1); %density of the liquid
    eta0  = pa(2); %viscosity of the liquid
    gamma = pa(3); %surface tension
    theta = pa(4); %contact angle
    grav     = pa(31); %gravity
    
    %geometrical factors
    H1    = pa(5); % High of the box
    R     = pa(6); %Radius of the box
    H2    = pa(7); %upper high
    RminD = pa(25); % Minimum radius of coil
    RmaxD = pa(26);
    ZminD = pa(27); % Minimum height of coil
    ZmaxD = pa(28);
    Rout  = pa(29); % Outer radius of domain
    zout  = pa(30); % Minimum height of domain   
    
    %parameter of the ferrofluid (Magnetization curve)
    aM    = pa(8); %
    bM    = pa(9); %
    cM    = pa(10);%
    dM    = pa(11);%
    eM    = pa(12);%
    
    %coi parameters
    Ib    = pa(13); %intensity
    Sb    = pa(14);
    
    %parameters need to compute BC1678
    m1    = pa(15); % Magnetic dipole moment 1 (Am2)
    mz1   = pa(16);% Magnetic dipole height 1 (m)
    m2    = pa(17);% Magnetic dipole moment 2 (Am2)
    mz2   = pa(18);% Magnetic dipole moment 2 (Am2)
    
    %M
    Mmagnet = pa(19);
    
    %Parameters for varible viscosity
    zeta=pa(20);
    tau=pa(21);
    
    % Fluid volume
    VF=pa(22);
    
    Fcoil=pa(23);
    Gcoil=pa(24);
    
    %% Blocks evaluation
    %
    % display
    
    
    %Be carefull compile  each block separatly
    %block to compile
    switch lb
        case 1
            blockAnew2D;
        case 2
            blockBnew2D;
        case 3
            blockCnew2D;
        case 4
            blockDnew2D;
        case 5
            blockEnew2D;
    end
    stop
    matrixgen
end

if (lstart==0)
    %creating mesh in A and defining collocation matrices
    
    %% Block A: Concentrated meshing with 4th order in r and z
    [z0A,dz0A,dzz0A] = finitasdoblegood(nzA,1);
    [r0A,dr0A,drr0A] = finitasdoblegood(nzA,1);
    
    ntA              = nrA*nzA;
    r1A              = repmat(r0A', [1 nzA]);
    z1A              = repmat(z0A, [nrA 1]);
    
    %initial variables
    fA               = R*z0A;
    fA               = repmat(fA, [nrA 1]);
    gA               = H*ones(nrA,nzA);%+0.1*cos(2*pi*fA);
    gaxisA           = 0*gA;
    
    mass
    VF               = hnominal*R^2/2;
    
    %
    pref             = 1;
    f0               = z0A*R;
    g0               = ones(1,nzA)*hnominal;
    
    Bond             = 0;
    options          = optimset('Display','iter');
    x0               = [f0,g0,pref];
    
    % Solve interface
    %x                = fsolve(@(x) Fshapep(x,nzA,z0A,dz0A,dzz0A,R,VF,theta_fe),x0,options);
    x                = fsolve(@(x) Fshapepgravity(x,nzA,z0A,dz0A,dzz0A,R,VF,theta_fe,grav_fe,rho_fe),x0,options);
    nz               = nzA;
    f                = x(1:nz);
    h                = x(nz+1:2*nz);
    fA               = repmat(f, [nrA 1]);
    gA               = repmat(h, [nrA 1]);
    %  pref=x(2*nz+1);
    % plot(f,h,'g-')
    %
    %      pref=1;
    %      f0=z0A*R;
    %      g0=ones(1,nzA)*hnominal;
    %
    %  Bond=0;
    %   options=optimset('Display','iter');
    %   V=g0*0;
    %    x0=[f0,g0,V,pref];
    % %
    %  x=fsolve(@(x) Fshapep1(x,nzA,z0A,dz0A,dzz0A,Bond,R,VF,theta_fe),x0,options);
    %  nz=nzA;
    %     f=x(1:nz);
    %  h=x(nz+1:2*nz);
    % hold on
    %  pref=x(2*nz+1);
    %  plot(f,h,'x')
    %  stop
    % %
    
    %axial velocity
    gaxisA   = zeros(nrA,nzA);
    wA       = zeros(nrA,nzA);
    pA       = zeros(nrA,nzA);
    uA       = zeros(nrA,nzA);
    
    %streamfunction
    psiA     = 0*fA;
    phiA     = 0*fA;
    etaA     = eta_fe*ones(nrA,nzA);
    
    %magnetization
    CA       = 0*fA;
    
    %volume
    VA       = CA;
    
    %initial mesh in Block A
    rA       = fA;
    zA       = r1A.*gA;
    
    %maping points
    r1A      = reshape(r1A,nrA*nzA,1);
    z1A      = reshape(z1A,nrA*nzA,1);
    
    ntB      = nrB*nzB;
    
    %% Block B: Concentrated meshing with 4th order in r, and liquid-concentrated mesh in z
    [r0B,dr0B,drr0B] = Chevitanh4th_i(nrB,0,1,3);                %%%%%%%% 10/03/2021
    [z0B,dz0B,dzz0B] = finitasdoblegood(nzB,1);                 %%%%%%%% 18/09/2020
    
    r1B  = repmat(r0B', [1 nzB]);                        %Manual
    z1B  = repmat(z0B, [nrB 1]);                         %Manual
    fB   = repmat(f, [nrB 1]);
    gB   = repmat(h, [nrB 1]);
    
    %initial mesh....
    rB   = fB;
    zB   = gB+r1B.*(H1-gB);
    
    %creating Linea10B(bottom wall i=1 1<j<nzA)
    %converting in a vector
    
    psiB = 0*ones(nrB,nzB);
    phiB = 0*ones(nrB,nzB);
    
    %mapping points
    r1B  = reshape(r1B,nrB*nzB,1);
    z1B  = reshape(z1B,nrB*nzB,1);
    
    %% Block C: Concentrated meshing with 4th order in r, and uniform 4th order mesh in z
    [z0C,dz0C,dzz0C] = finitasdoblegood(nzC,1);                %%%%%%%% 18/09/2020
    [r0C,dr0C,drr0C] = finitas4ordentotal(nrC,1);
    
    ntC  = nrC*nzC;
    r1C  = repmat(r0C', [1 nzC]);
    z1C  = repmat(z0C, [nrC 1]);
    fC   = repmat(f, [nrC 1]);
    psiC = 0*ones(nrC,nzC);
    phiC = 0*ones(nrC,nzC);
    
    %initial mesh....
    rC   = fC;
    zC   = H1+r1C.*(H2-H1);
    
    %% Block E: Concentrated meshing with 4th order in r, and single-side concentrated 4th order mesh in z
    [z0E,dz0E,dzz0E]=finitasdoblegood(nzE,1);                %%%%%%%% 18/09/2020
    %[r0E,dr0E,drr0E]=finitas4ordentotal(nrE,1);
    [r0E,dr0E,drr0E]=Chevitanh4th(nrE,0,1,3);
    
    ntE  = nrE*nzE;
    r1E  = repmat(r0E', [1 nzE]);
    z1E  = repmat(z0E, [nrE 1]);
    fE   = repmat(f, [nrE 1]);
    psiE = 0*ones(nrE,nzE);
    phiE = 0*ones(nrE,nzE);
    
    %initial mesh....
    rE   = fE;
    zE   = zout-r1E.*zout;
    
    %% Pointers and creation of elliptical 4th order mesh in D
    gettingLinesandVertices
    ntD   = nrD*nzD;
    
    % Coil subdomain
    Fcoil = reshape(FD,nrD,nzD);
    Gcoil = reshape(GD,nrD,nzD);
    
    % Initial values
    psiD  = 0*ones(nrD,nzD);
    phiD  = 0*ones(nrD,nzD);
    MzD   = 0*ones(nrD,nzD);
    MzD(CoilD) = Mmagnet_fe;
    MzD(CoilD_contour) = Mmagnet_fe/2;
    r1D   = repmat(r0D', [1 nzD]);
    z1D   = repmat(z0D, [nrD 1]);
    
    %% Plot mesh
    plottingmesh
    
    %% Converting Matrices of block A, B, C and E
    [dd0rA,dd0rrA,dd0zA,dd0zzA,deltaA,dd0rzA]= matricesconversor(nrA,nzA,dr0A,drr0A,dz0A,dzz0A);
    [dd0rB,dd0rrB,dd0zB,dd0zzB,deltaB,dd0rzB]= matricesconversor(nrB,nzB,dr0B,drr0B,dz0B,dzz0B);
    [dd0rC,dd0rrC,dd0zC,dd0zzC,deltaC,dd0rzC]= matricesconversor(nrC,nzC,dr0C,drr0C,dz0C,dzz0C);
    [dd0rD,dd0rrD,dd0zD,dd0zzD,deltaD,dd0rzD]= matricesconversor(nrD,nzD,dr0D,drr0D,dz0D,dzz0D);
    [dd0rE,dd0rrE,dd0zE,dd0zzE,deltaE,dd0rzE]= matricesconversor(nrE,nzE,dr0E,drr0E,dz0E,dzz0E);
    
    % Physical Parameters
    rho   = rho_fe;   % Density of the liquid
    eta0  = eta_fe;   % Viscosity of the liquid
    gamma = gamma_fe; % Surface tension
    theta = theta_fe; % Contact angle
    zeta  = zeta_fe;  % Magnetic viscosity coefficients
    tau   = tau_fe;
    grav  = grav_fe;
    
    %geometrical factors
    %     H1=pa(5); % Height of the box
    %     R=pa(6);  % Radius of the box
    %     H2=pa(7); % Upper height
    
    %parameter of the ferrofluid (Magnetization curve)
    aM   = aM_fe;
    bM   = bM_fe;
    cM   = cM_fe;
    dM   = dM_fe;
    eM   = eM_fe;
    
    %coil parameters
    Ib   = Ib_fe * Nturns; %intensity
    Sb   = (max(RmaxD)-min(RminD))*(max(ZmaxD)-min(ZminD));
    
    %Magnet magnetization vector (if any)
    Mmagnet = Mmagnet_fe;
    
    %inital M and Mz=0
    M    = 0;
    Mz   = 0;
    
    %getting magnetic fieldsARREGLAR*
    [m1, mz1, m2, mz2] = computedipoles(M, Mz, Ib, Mmagnet, RminD, RmaxD, ZminD,ZmaxD);
    
    % parameters
    pa=[rho;eta0;gamma;theta;H1;R;H2;aM;bM;cM;dM;eM;Ib;Sb;m1;mz1;m2;mz2;Mmagnet;zeta;tau;VF;0;0;RminD;RmaxD;ZminD;ZmaxD;Rout;zout;grav];
    
    % % % x0=eval([order ']']);
    order='[';
    for j=1:nbl
        bl = list_block{j};
        NVAR = eval(['NV' bl]);
        for i=1:NVAR
            lv = eval(['list_var_' bl '{' num2str(i) '}']);
            order = [order 'reshape(' lv ',nt' bl ',1);'];
        end
    end
    x0=eval([order ']']);
    x0m = x0;
    x0mm = x0m;
    
else
    
    load(name_prev_sol, '-regexp','^(?!Ib_fe$|Nturns$|dt$|Nsteps$|miniter$|maxerror$|tiempo$|theta_fe$|gamma_fe$|hmin$|eta_fe$|zeta_fe$|tau_fe$|Mmagnet_fe$|filesave$|ksave$).');
    
    % EXTRACTION OF PHYSICAL PARAMETERS IN VECTOR OF PARAMETERS
    % Physical properties of the liquid
    rho   = pa(1); %density of the liquid
    eta0  = pa(2); %viscosity of the liquid
    gamma = pa(3); %surface tension
    theta = pa(4); %contact angle
    grav     = pa(31);%gravity
    
    % Geometrical factors
    H1    = pa(5); % High of the box
    R     = pa(6); %Radius of the box
    H2    = pa(7); %upper high
    RminD = pa(25); % Minimum radius of coil
    RmaxD = pa(26);
    ZminD = pa(27); % Minimum height of coil
    ZmaxD = pa(28);
    Rout  = pa(29); % Outer radius of domain
    zout  = pa(30); % Minimum height of domain  
    
    % Parameter of the ferrofluid (Magnetization curve)
    aM    = pa(8); %
    bM    = pa(9); %
    cM    = pa(10);%
    dM    = pa(11);%
    eM    = pa(12);%
    
    % Coil parameters
    Ib    = pa(13); %intensity
    Sb    = pa(14);
    
    % Parameters needed to compute BC1678
    m1    = pa(15); % Magnetic dipole moment 1 (Am2)
    mz1   = pa(16); % Magnetic dipole height 1 (m)
    m2    = pa(17); % Magnetic dipole moment 2 (Am2)
    mz2   = pa(18); % Magnetic dipole moment 2 (Am2)
    
    %M
    Mmagnet = pa(19);
    
    %Parameters for varible viscosity
    zeta  = pa(20);
    tau   = pa(21);
    
    % Fluid volume
    VF    = pa(22);
    
    % UPDATE OF KEY PARAMETERS IN VECTOR OF PARAMETERS
    % New viscosity
    eta0  = eta_fe;
    pa(2) = eta0;
    
    % New surface tension
    gamma = gamma_fe;
    pa(3) = gamma;
    
    % New contact angle
    theta = theta_fe;
    pa(4) = theta;
    
    % New magnet magnetization
    Mmagnet = Mmagnet_fe;
    pa(19)  = Mmagnet;
    
    % New intensity
    Ib     = Ib_fe*Nturns;
    pa(13) = Ib;
    
    % New dipole BCs
    M      = 0;
    Mz     = 0;
    [m1, mz1, m2, mz2] = computedipoles(M, Mz, Ib, Mmagnet, RminD, RmaxD, ZminD,ZmaxD);
    
    pa(15) = m1;
    pa(16) = mz1;
    pa(17) = m2;
    pa(18) = mz2;
    
    % New magnetic viscosity
    zeta   = zeta_fe;
    pa(20) = zeta;
    tau    = tau_fe;
    pa(21) = tau;
    
end


xotoredeable % changes the unknown format: from row to an intelligible one.


%% Time integration
% The generic time integration is a second order backward,
%
% $$ \dot{y}(t^n) = -\frac{2 \Delta t_1 +\Delta t_2}{\Delta t_1(\Delta t_1
% +\Delta t_2)} y^n + \frac{\Delta t_1 +\Delta t_2}{\Delta t_1 \Delta t_1
% } y^{n-1} + \frac{\Delta t_1}{\Delta t_2 (\Delta t_1
% +\Delta t_2)} y^{n-2}  $$
%
% In case of equal timestep it reduces to,
%
% $$ \dot{y}(t^n) = \frac{-3 y^n + 4 y^{n-1} - y^{n-2}}{2 \Delta t}  $$
%
% In this case since we search a steady case (It is not a time evolution
% problem) we set $\Delta t$ to infinity.
%

dt1=dt;
ksave=0;
mass;
for ll=1:Nsteps
    %
    % _tiempo_ is an array where computed instants are saved,
    %matlab
    if (ll==1)
        tiempo(ll)=0.28;
    else
        tiempo(ll)=tiempo(ll-1)+dt;
    end
    
    error=1e9;
    iter=0;
    
    % Display Iteration start message
    fprintf('Iteration %d of %d started.\n',ll, Nsteps)
    
    % Iterations proceed while the larger error is above 1e-9
    
    while (error > maxerror || iter < miniter)
        %%
        % Increase the iteration counter
        iter=iter+1;
        %%
        % Construct $\mathbf{A}$ and $\mathbf{B}$
        
        matrixABCDEfull
        
        % M1=reshape(MA,nrA,nzA);
        % M1z=reshape(MAz,nrA,nzA);
        % M1(:,1)=(4*M1(:,2)-M1(:,3))/3;
        % M1z(:,1)=(4*M1z(:,+2)-M1z(:,3))/3;
        %
        % MA=reshape(M1,ntA,1);
        % MAz=reshape(M1z,ntA,1);
        
        %          c1=isnan(a)
        %          c2=isnan(b)
        %          c3=isinf(a)
        %          c3=isinf(b)
        xotoredeable
        
        % Solve the linear system of equation to compute the correction $\delta \mathbf{x}$.
        dxa=a\b;
        error=max(abs(dxa));
        
        if(error > 10^11)
            cprintf('Errors', 'Solution diverges. Try a smaller time step size or a progressive increase of the external parameters\n')
            return
        end
        %%
        % Apply the correction
        n=length(dxa);
        %         ni=NVA*ntA+NVB*ntB+NVC*ntC;
        %         x0(ni+1:ni+n)=x0(ni+1:ni+n)+dxa;
        x0(1:n)=x0(1:n)+dxa;
        xotoredeable
        
        %puntos fisicos
        
        zA=reshape(r1A,nrA,nzA).*(gA-gaxisA);
        rA=fA;
        
        zB=gB+reshape(r1B,nrB,nzB).*(H1-gB);
        rB=fB;
        
        rC=fC;
        zC=H1+reshape(r1C,nrC,nzC).*(H2-H1);
        
        rD=FD;
        zD=GD;
        
        zE=zout-reshape(r1E,nrE,nzE).*zout;
        rE=fE;
        
        %plottingmesh
        
        
        
        %           plot(z0A,CA)
        %         stop
        %correcting mass
        % mass
        
        % % % x0=eval([order ']']);
        order='[';
        for j=1:nbl
            bl = list_block{j};
            NVAR = eval(['NV' bl]);
            for i=1:NVAR
                lv = eval(['list_var_' bl '{' num2str(i) '}']);
                order = [order 'reshape(' lv ',nt' bl ',1);'];
            end
        end
        x0=eval([order ']']);
        
        MA       = reshape(MA,nrA,nzA);
        MAz      = reshape(MAz,nrA,nzA);
        
        MA(:,1)  = MA(:,2);
        MAz(:,1) = MAz(:,2);
        dfAs     = fA*dz0A';
        
        %computing M and Mz
        M=0;
        Mz=0;
        for j=2:nzA
            for i=2:nrA
                app = gA(i,j)*fA(i,j)*MA(i,j)*dfAs(i,j);
                amp = gA(i-1,j)*fA(i-1,j)*MA(i-1,j)*dfAs(i-1,j);
                apm = gA(i,j-1)*fA(i,j-1)*MA(i,j-1)*dfAs(i,j-1);
                amm = gA(i-1,j-1)*fA(i-1,j-1)*MA(i-1,j-1)*dfAs(i-1,j-1);
                M   = M+0.25*(r0A(i)-r0A(i-1))*(z0A(j)-z0A(j-1))*(app+amp+apm+amm);
                app = gA(i,j)*fA(i,j)*MAz(i,j)*dfAs(i,j);
                amp = gA(i-1,j)*fA(i-1,j)*MAz(i-1,j)*dfAs(i-1,j);
                apm = gA(i,j-1)*fA(i,j-1)*MAz(i,j-1)*dfAs(i,j-1);
                amm = gA(i-1,j-1)*fA(i-1,j-1)*MAz(i-1,j-1)*dfAs(i-1,j-1);
                Mz  = Mz+0.25*(r0A(i)-r0A(i-1))*(z0A(j)-z0A(j-1))*(app+amp+apm+amm);
            end
        end
        M  = M*2*pi;
        Mz = Mz*2*pi;
        %end computing M Mz
        
        %getting magnetic fields
        [m1, mz1, m2, mz2] = computedipoles(M, Mz, Ib, Mmagnet, RminD, RmaxD, ZminD,ZmaxD);
        
        % parameters
        pa=[rho;eta0;gamma;theta;H1;R;H2;aM;bM;cM;dM;eM;Ib;Sb;m1;mz1;m2;mz2;Mmagnet;zeta;tau;VF;0;0;RminD;RmaxD;ZminD;ZmaxD;Rout;zout;grav];
       
        % Modify meshes in A+B+C+D
        % A
        rA   = reshape(fA,ntA,1);
        zA   = reshape(r1A,ntA,1).*(reshape(gA-gaxisA,ntA,1));
        
        % B
        rB   = reshape(fB,ntB,1);
        zB   = reshape(gB,ntB,1)+reshape(r1B,ntB,1).*(H1-reshape(gB,ntB,1));
        
        % C
        rC   = reshape(fC,ntC,1);
        zC   = H1+r1C.*(H2-H1);
        
        % Display message
        fprintf('Subiteration %d: error = %f, hnominal = %f\n', iter,error,hnominal)
        
        %         velocities(za,ra,zb,rb,zc,rc,zd,rd,ze,re,psiA,psiB,psiC,psiD,psiE)
        
    end    
    
    %cambiamos x0
    xotoredeable
    mass;
    
    % % % x0=eval([order ']']);
    order='[';
    for j=1:nbl
        bl = list_block{j};
        NVAR = eval(['NV' bl]);
        for i=1:NVAR
            lv = eval(['list_var_' bl '{' num2str(i) '}']);
            order = [order 'reshape(' lv ',nt' bl ',1);'];
        end
    end
    x0=eval([order ']']);
    x0mm = x0m;
    x0m  = x0;
    dt2  = dt1;
    dt1  = dt;
    ksave=ksave+1;
    
    %% Represent
%     visualizationHerrada(rA, rB, rC, FD, rE, zA, zB, zC, GD, zE, nrA, nrB,...
%         nrC, nrD, nrE, nzA, nzB, nzC, nzD, nzE, phiA, phiB, phiC, phiD, phiE,psiD, HrA, HrB,...
%         HrC, HrD, HrE, HzA, HzB, HzC, HzD, HzE, MAz0, MAr0, gA, fA, Forcer,...
%         Forcez, pA, etaA, Ib_fe, gaxisA, dd0rD,dd0rrD, dd0rzD,...
%         dd0zD, dd0zzD,ZminD,ZmaxD,RminD,RmaxD)
    
    %%  Save solution
    if ksave==1 && filesave ==1
        clear a b
        % Case name
        finame = ['soluciones/DoubleCoil-DYT-',date,'_Ib' po2com(Ib) '_zeta' po2com(zeta),...
            'tau' po2com(tau) 'nrA' num2str(nrA)     'nzA' num2str(nzA)  'nrB' ...
            num2str(nrB) 'nrC' num2str(nrC) 'nrE' num2str(nrE)    'nzE' num2str(nzE) '_t' num2str(tiempo(ll)) '.mat'];
        
        % Save file
        hmin(ll)=gA(nrA,1)-gaxisA(nrA,1);
        save(finame);
        
        % Update index
        ksave=0;
        
        %plot(tiempo(ll), gA(nzA,1),'k.')
        %plot(fA(nrA,:),gA(nzA,:),'linewidth',0.5)
        pause(0.001)
        
        % Display command line information
        hold on, plot(fA(nrA,:),gA(nrA,:),'c','linewidth',1)
        drawnow
        fprintf('Iteration %d of %d completed. Results saved as %s.mat\n-----\n',ll, Nsteps, finame)
    end
end


%% STREAM

% % % x0=eval([order ']']);
pA=zeros(nrA,nzA);
order='[';
for j=1:nbl
    bl = list_block{j};
    NVAR = eval(['NV' bl]);
    for i=1:NVAR
        lv = eval(['list_var_' bl '{' num2str(i) '}']);
        order = [order 'reshape(' lv ',nt' bl ',1);'];
    end
end
x0=eval([order ']']);
x0m = x0;
x0mm = x0m;

xotoredeable
r1A=repmat(r0A', [1 nzA]);

ZA=(gA-gaxisA).*r1A;

%contourf(fA,ZA,wA)
%     contourf(gA,fA,wA)

matrixAstream
xotoredeable
psi0=a\b; %MEJORAR EL TIEMPO DE INVERSÓN MENDIANTE MÉTODOS ITERATIVOS...
psi0=reshape(psi0,nrA,nzA);
figure
plot(fA(nrA,:),gA(nrA,:),'r-','LineWidth',2)
hold on
% stop
psi3a=zeros(nrA,nzA);
for j=1:nzA
    for i=1:nrA
        if (psi0(i,j)<0)
            psi3a(i,j)=-(-psi0(i,j))^(1/3);
        else
            psi3a(i,j)=(psi0(i,j))^(1/3);
        end
    end
end

c1=min(min(psi3a));
c2=max(max(psi3a));
% c1=-0.015;
% c2=0.015;




ns=24
for i=1:ns+1
    %va(i)=(i/15)^2*c1;
    vsa(i)=((i-1)/ns)*c1;
    vsb(i)=((i-1)/ns)*c2;
end

contour(fA,ZA,psi3a ,vsa,'c-')
hold on
contour(fA,ZA,psi3a,vsb,'b-')

