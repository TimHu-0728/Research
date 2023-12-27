function [omegar,omegai] = eigenfcn(Ib_fe)
%% Solver de Miguel Angel Herrada
% V4 (20/04/2021)

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

% File
name_sol = ['soluciones/DYT-05-May-2021_Ib',num2str(Ib_fe*200),...
    '_zeta1,2571e-05tau5,5194e-07nrA101nzA101nrB41nrC21nrE61nzE101.mat'];

%% Details of the block, variables and connections
%
%
% Fix the number of blocks (subdomains). In each block it has to be defined
% the number of variables and the number of symbolic variables.
%
%%
% number of parameters

Np = 30;

% To construct the numerical scheme we need to use the list of variables of
% each block (*Important: In the same order as defined in the corresponding
% block*).

list_block = {'A' 'B' 'C' 'D' 'E'}; %see scheme

list_var_A = {'wA', 'uA', 'pA','phiA','fA' 'gA' 'CA' 'etaA' 'psiA', 'VA' 'gaxisA'}; %phiA is the streamfunction for H CA is magnetic susceptibility
list_var_B = {'phiB','fB' 'gB' 'psiB'}; %phiB is the streamfunction for H
list_var_C = {'phiC','fC','psiC'}; %phiC is the streamfunction for H
list_var_D = {'phiD','psiD' 'FD' 'GD' 'MzD'};%phiB is the streamfunction for H ellipitical mesh
list_var_E = {'phiE' ,'fE','psiE' }; %magentic streamfuction CE is magnetic susceptibility

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

omega0=1;
npeigen=36;
%%%%%% ->
load(name_sol);

%% Form the array  of dimensionless parameters
%phisicla properties of the liquid
rho   = pa(1); %density of the liquid
eta0  = pa(2); %viscosity of the liquid
gamma = pa(3); %surface tension
theta = pa(4); %contact angle

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
x0m=x0;
x0mm=x0;
%hold on


%%
% The grid can be refined/coarsened. If _interpol_ = 1 the existing solution
% will be interpolated to the new grid.
xotoredeable%% Eigen value calculation
Fcoil=reshape(FD,nrD*nzD,1);
Gcoil=reshape(GD,nrD*nzD,1);

% Plot mesh
% plottingmesh
% hold on
%  rx=reshape(rA,nrA,nzA);
% %    plot(tiempo,hmin,'-')
% %    stop
%     rx=reshape(rA,nrA,nzA);
%         zx=reshape(zA,nrA,nzA);
%         Forcerx=reshape( Forcer,nrA,nzA);
%         Forcezx=reshape( Forcez,nrA,nzA);
%          %quiver(rx,zx, Forcerx, Forcezx)
%
%            contourf(rx,zx,wA)
%    stop
%         stop
%         hold on
% %         plot(fA(nrA,:),Forcezx(nrA,:),'r-')
%         stop
%  contourf(rx,zx,(Forcerx.^2+Forcezx.^2).^0.5)

% stop
%     plot(fA(nrA,:),wA(nrA,:),'-')
%      stop
%  stop
% We can automate a litte bit the search of eigensvalue bu means of a for
% loop indicating initial frequency , final frequency and number of points
% between both.
%
Nom = 1;
omegai0 =5
omegaf0=5;
deltaom = (omegaf0-omegai0)/(Nom-1);

if(isnan(deltaom))
    deltaom = 0.;
end

for iom=1:Nom
    %omega0=omegai0+deltaom*(iom-1);
    omega0=10-0.2*1i
    %%
    % We compute matrix a and b with _matrix_eigen_
    hold on
    matrixABCDE_eigen
    
    %%
    % Calling later to the eigen solver subroutine. eigensolver returns, _d_,
    % the sixteen larger eigenvalues and the corresponding 16 eigenvector _V1_
    % that are stored for a later post-processing at Vf and df
    
    [V1,d]=eigensolver(a,b,omega0,npeigen);
    
    Vf(:,:,iom)=V1;
    df(:,:,iom)=d;
end

ns=length(d(:,:));

%% Seaech for the most unstable, i.e. that with larger imaginary part
%

gi=-10000;

for is=1:Nom
    omegar(:,is)=real(diag(df(:,:,is)));
    omegai(:,is)=imag(diag(df(:,:,is)));
    
    for ll=1:ns
        %if(omegai(ll,is)>gi)&&(omegar(ll,is) > 1e-4)
        if(omegai(ll,is)>gi)&&(omegar(ll,is)>-0.001) &&(omegai(ll,is)<100)
            lt=ll;
            it=is;
            gi=omegai(ll,is);
            gr=omegar(ll,is);
        end
    end
    
end
