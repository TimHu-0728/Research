function MainmeshD
%% Solver de Miguel Angel Herrada
% V2.2 (04/11/2020)
% clear all

%% Initial setup
% Lets clean before start and decide where to locate jacobians functions
% and auxiliary elements

% Path of jacobians
path_jacobian  = './jacobians2D/';

% Compile? 
leq=0;

%
%%
% number of parameters

Np = 6;

% To construct the numerical scheme we need to use the list of variables of
% each block (*Important: In the same order as defined in the corresponding
% block*).

list_block = { 'D' }; %see scheme


list_var_D = { 'FD' 'GD'};%phiB is the streamfunction for H ellipitical mesh







[basura nbl ] = size(list_block);


% The list of derivatives. We have assumend that the number and distribution
% is the same in all the block. This information is relevant in _Matrix_
%

list_der_D = {'', 'r0', 'z0', 'rr0', 'zz0', 'rz0' ,'t0'}; %mapping derivatives asumming y=FA(r0,z0,t0,time), x=GA(r0,z0,t0,time), z=t0

list_der=list_der_D;
for i=1:nbl
    bl = list_block{i};
    order = ['[bas NV' bl '] = size(list_var_' bl ');'];
    evalc(order);
    %     order = ['[bas ND' bl '] = size(list_der);' ];
    %     evalc(order);
    
    order = ['[bas ND' bl '] = size(list_der_' bl ');' ];
    evalc(order);
end

list_con_D=[];
%% Compute symbolic functions and construction of matrix.m
if (leq==1)
    
    for i=1:Np
        pa(i,1)=sym(['pa_',num2str(i)],'real');
    end
    
    %% Form the array  of dimensionless parameters
    %phisicla properties of the liquid
    Rout=pa(1);
    H2=pa(2);
    zout=pa(3);
    Fcoil=pa(4);
    Gcoil=pa(5);
    alphaf=pa(6);
    
    %% Blocks evaluation
    %
    % display
    
    
    %Be carefull compile  each block separatly
    
    blockDnew2D_mg;
    stop
end


miniter=100;


%Reading initial mesh

load('initmesh')
ntD=nrD*nzD;

[dd0rD,dd0rrD,dd0zD,dd0zzD,deltaD,dd0rzD]= matricesconversor_mg(nrD,nzD,dr0D,drr0D,dz0D,dzz0D);

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


hold on
%      %ploting mesh E
%  for j=1:nzD
%  plot(FD(:,j),GD(:,j),'o')
%  end
%  for i=1:nrD
%  plot(FD(i,:),GD(i,:),'o')
%  end


%  plot(FD(:,imin),GD(:,imax),'x')
%   plot(FD(nx,:),GD(nx,:),'co')

Fcoil=reshape(FD,ntD,1);
Gcoil=reshape(GD,ntD,1);
Nsteps=1;
maxerror=10^7;

pa=zeros(6,1);
pa(1)=Rout;
pa(2)=H2;
pa(3)=zout;
alphaf=0;


for liter=1:2 % liter=1 mallado no eliptioc liter=2 mallado eliptico
    
    
    iter=0;
    
    error=1e9;
    
    %     if (liter==3)
    %     [z0D,dz0D,dzz0D]=Chevitanh2th(nzD,0,1,1.5);
    %      [dd0rD,dd0rrD,dd0zD,dd0zzD,deltaD,dd0rzD]= matricesconversor_mg(nrD,nzD,dr0D,drr0D,dz0D,dzz0D);
    %     end
    
    % Iterations proceed while the larger error is above 1e-9
    
    while (error > 1e-2 && iter < 300)
        %%
        % Increase the iteration counter
        iter=iter+1;
        %%
        % Construct $\mathbf{A}$ and $\mathbf{B}$
        
        
        %         stop
        matrixD_mg
        
        xotoredeable_mg
        
        
        
        % Solve the linear system of equation to compute the correction $\delta \mathbf{x}$.
        dxa=a\b;
        error=max(abs(dxa))
        
        
        
        
        x0=x0+dxa;
        xotoredeable_mg
        
        
    end
    
    
end


% hold on
% %      %ploting mesh E
% for j=1:nzD
%     plot(FD(:,j),GD(:,j),'b-')
% end
% for i=1:nrD
%     plot(FD(i,:),GD(i,:),'b-')
% end


save('./meshgenerator/finalmesh','FD','GD', 'r0D','z0D', 'dr0D','drr0D','dz0D','dzz0D',...
    'nrD','nzD', 'ndD', 'CoilD','Linea10D','Linea7D','Linea8D','Linea20D',...
    'Linea13D','Linea21D','Linea19D', 'VertexDBC','VertexDEA','VertexDAB',...
    'zout','H2','R','Rout' )
