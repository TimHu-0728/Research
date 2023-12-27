

    
    load('hola')
    
    
    
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
 name='soluciones/dt2casezeta0,00242tau1e-06Ib200nrA71nzA71nrB81nrC71nrE51nzE51'
 load(name)
       %% Form the array  of dimensionless parameters
    %phisicla properties of the liquid
    rho   = pa(1); %density of the liquid
    eta0   = pa(2); %viscosity of the liquid
    gamma = pa(3); %surface tension
    theta = pa(4); %contact angle
    %geometrical factors
    H1    = pa(5); % High of the box
    R     = pa(6); %Radius of the box
    H2    = pa(7); %upper high
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
    zeta=pa(20);
    tau=pa(21);
    dt=0.01
    clear tiempo
    clear hmin
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
figure, hold on
dt1=dt;
ksave=0;
mass  
for ll=1:Nsteps
    %
    % _tiempo_ is an array where computed instants are saved,
    %matlab
    if (ll==1)
        tiempo(ll)=0;
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
        
          matrixABCDE
%         stop
        [a,b,MA,MAz,HrA,HzA,HrB,HzB,HrC,HzC,HrD,HzD,HrE,HzE,MAr0,MAz0,Forcer,Forcez]=...
            matrixABCDEfull(x0,x0m,x0mm,dt,dt1,pa,ntA,ntB,ntC,ntD,ntE,...
            NVA,NDA,NVB,NDB, NVC,NDC, NVD,NDD,NVE,NDE,...
            z1A,r1A,z1B,r1B,z1C,r1C,zD,rD,rE,zE,...
            dd0rA,dd0zA,dd0rrA,dd0zzA,dd0rzA...
            ,dd0rB,dd0zB,dd0rrB,dd0zzB,dd0rzB...
            ,dd0rC,dd0zC,dd0rrC,dd0zzC,dd0rzC...
            ,dd0rD,dd0zD,dd0rrD,dd0zzD,dd0rzD...
            ,dd0rE,dd0zE,dd0rrE,dd0zzE,dd0rzE...
            , ndA,ndB,ndC,ndE,nid...
            ,Linea11A,Linea11B ...
            ,Linea12C,Linea12B...
            ,Linea13D,Linea13C...
            ,Linea9Ab,Linea9Ar...
            ,Linea9Db,Linea9Dr...
            ,Linea14D,Linea15D,Linea16D,Linea17D...
            ,Linea14E,Linea15E,Linea16E,Linea17E...
            ,Linea10D,Linea10B);
        
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
        x0=x0+dxa;
        %x0(1:n)=x0(1:n)+dxa;
        xotoredeable
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
        [m1, mz1, m2, mz2] = computedipoles(M, Mz, Ib, Mmagnet, rE(Linea15E(1)), rE(Linea17E(1)), zE(Linea14E(1)),zE(Linea16E(1)));
        
        % parameters
   pa=[rho;eta0;gamma;theta;H1;R;H2;aM;bM;cM;dM;eM;Ib;Sb;m1;mz1;m2;mz2;Mmagnet;zeta;tau];
     
        
        % Modify meshes in A+B+C+D
        % A
        rA   = reshape(fA,ntA,1);
        zA   = reshape(r1A,ntA,1).*reshape(gA,ntA,1);
        
        % B
        rB   = reshape(fB,ntB,1);
        zB   = reshape(gB,ntB,1)+reshape(r1B,ntB,1).*(H1-reshape(gB,ntB,1));
        
        % C
        rC   = reshape(fC,ntC,1);
        zC   = H1+r1C.*(H2-H1);
        
        % Changing points in D
        rD(Linea9Db) = rA(Linea9Ab);
        zD(Linea9Db) = zA(Linea9Ab);
        rD(Linea9Dr) = rA(Linea9Ar);
        zD(Linea9Dr) = zA(Linea9Ar);
        rD(Linea10D) = rB(Linea10B(2:length(Linea10B)));
        zD(Linea10D) = zB(Linea10B(2:length(Linea10B)));
        
        % Compute neighbours
        id           = find(rD >= 0 & rD <= (max(rD(Linea9Db)) + 0.007) & zD >= -0.007 & zD <= (max(zD(Linea10D)) + 0.015));
        [la]         = setneighbours_FEM2(rD, zD, nid, Nneig, WF, la, id');
        
        %correcting matrix
        [dd0rD,dd0zD,dd0rrD,dd0zzD,C] = collocationmatrixMLS(la,rD,zD);
        dd0rzD                        = dd0rD*dd0zD;
        
        % Display message
        fprintf('      Subiteration %d: error = %f, hnominal = %f\n', iter,error,hnominal)
    end
    
    % REPRESENT RESULTS
%     visualization(rA, rB, rC, rD, rE, zA, zB, zC, zD, zE, nrA, nrB,...
%     nrC, nrE, nzA, nzB, nzC, nzE, phiA, phiB, phiC, phiD, phiE, HrA, HrB,...
%     HrC, HrD, HrE, HzA, HzB, HzC, HzD, HzE, MAz0, MAr0, gA, fA, P1, P2, la)
xotoredeable
%cambiamos x0
        
mass  
stop
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
    
    %%  Save solution
    if( ksave==1)
        % Case name
        finame = ['soluciones/apolydt2casezeta' po2com(zeta), 'tau' po2com(tau) 'Ib' po2com(Ib),'nrA' num2str(nrA)     'nzA' num2str(nzA)  'nrB' num2str(nrB) 'nrC' num2str(nrC) 'nrE' num2str(nrE)    'nzE' num2str(nzE)];
        clear a b
        
        % Save file
        hmin(ll)=gA(nzA,1);
        save(finame);     
        
        % Update index
        ksave=0;
        
        plot(tiempo(ll), gA(nzA,1),'k.')
        %plot(fA(nrA,:),gA(nzA,:),'linewidth',0.5)
        pause(0.001)
        
        % Display command line information
        fprintf('Iteration %d of %d completed. Results saved as %s.mat\n-----\n',ll, Nsteps, finame)
    end
end


