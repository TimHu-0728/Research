function visualizationHerrada(rA, rB, rC, rD, rE, zA, zB, zC, zD, zE, nrA, nrB,...
    nrC, nrD, nrE, nzA, nzB, nzC, nzD, nzE, phiA, phiB, phiC, phiD, phiE,psiD, HrA, HrB,...
    HrC, HrD, HrE, HzA, HzB, HzC, HzD, HzE, dphidrD,dphidzD,dphidrrD,...
    dphidzzD,dphidrzD,dpsidrD,dpsidzD,dpsidrrD,dpsidzzD,...
    dpsidrzD,MAz0, MAr0, gA, fA, Forcer,Forcez, pA, etaA, Ib_fe, gaxisA,...
    dd0rD,dd0rrD, dd0rzD,dd0zD, dd0zzD,ZminD,ZmaxD,RminD,RmaxD)
% Represents variables of interest of the magnetic sloshing problem

%% Preprocessing
% Reshaping
ra=reshape(rA,nrA,nzA);
za=reshape(zA,nrA,nzA);
rb=reshape(rB,nrB,nzB);
zb=reshape(zB,nrB,nzB);
rc=reshape(rC,nrC,nzC);
zc=reshape(zC,nrC,nzC);
re=reshape(rE,nrE,nzE);
ze=reshape(zE,nrE,nzE);
rd=rD;
zd=zD;
% Forces
Forcerx=reshape(Forcer,nrA,nzA);
Forcezx=reshape(Forcez,nrA,nzA);

% Viscosity
etaA = reshape(etaA,nrA,nzA);

% Pressure
pA   = reshape(pA,nrA,nzA);

% Reshape potentials
% phiA=reshape(phiA,nrA*nzA,1);
% phiB=reshape(phiB,nrB*nzB,1);
% phiC=reshape(phiC,nrC*nzC,1);
% phiE=reshape(phiE,nrE*nzE,1);
%HzA=reshape(HzA,nrA,nzA);

% Magnetic field norm
HA=sqrt(HrA.^2+HzA.^2);
HB=sqrt(HrB.^2+HzB.^2);
HC=sqrt(HrC.^2+HzC.^2);
HD=sqrt(HrD.^2+HzD.^2);
HE=sqrt(HrE.^2+HzE.^2);
HA=reshape(HA,nrA,nzA);
HB=reshape(HB,nrB,nzB);
HC=reshape(HC,nrC,nzC);
HD=reshape(HD,nrD,nzD);
HE=reshape(HE,nrE,nzE);

% Magnetization grid norm
Mfield=(MAz0.^2+MAr0.^2).^0.5;
Mfield=reshape(Mfield,nrA,nzA);

% Move gA
gA    = gA -gaxisA;

%% Plot mesh
% figure
% hold on
% plot(r,z,'k.')
%
% % External and internal contours
% plot(P1(3:(P1(2)+2)), P1((P1(2)+3):(2*P1(2)+2)),'k-','LineWidth',2)
% plot(P2(3:(P2(2)+2)), P2((P2(2)+3):(2*P2(2)+2)),'r-','LineWidth',2)
%
% % Equilibrium surface
% plot(fA(nrA,:),gA(nzA,:)-gaxisA,'r','linewidth',2)
%
% % Format
% axis equal
% xlabel('r [m]')
% ylabel('z [m]')

%% Plot magnetic field
figure
hold on

% Magnetic field H (A/m)
velocities(za,ra,zb,rb,zc,rc,zd,rd,ze,re,HA,HB,HC,HD,HE)
h = colorbar; ylabel(h,'$H$ [A/m]', 'Interpreter', 'latex','FontSize', 14)

% Equilibrium surface
hold on
plot(fA(nrA,:),gA(nrA,:),'k','linewidth',2)
plot([0,0,1,1]*max(fA(nrA,:)),[gA(nrA,1),0,0,gA(nrA,end)],'k-','linewidth',1)

% Container and coil
plot([RminD,RmaxD,RmaxD,RminD,RminD],[ZminD,ZminD,ZmaxD,ZmaxD,ZminD],'k-','linewidth',1)

% Mesh points
%plot(rD,zD,'k.','MarkerSize',0.1)

% Format
axis equal
xlabel('r [m]')
ylabel('z [m]')

%% Plot magnetization
figure
hold on

% Magnetization (A/m)
contourf(ra,za,Mfield,linspace(1778.1,2185.7,20))
h = colorbar; ylabel(h,'$M$ [A/m]', 'Interpreter', 'latex','FontSize', 14)

% External and internal contours
plot([0,0,1,1]*max(fA(nrA,:)),[gA(nrA,1),0,0,gA(nrA,end)],'k-','linewidth',1)

% Equilibrium surface
plot(fA(nrA,:),gA(nrA,:),'k','linewidth',1)

% Format
axis equal
xlabel('r [m]')
ylabel('z [m]')

%% Forces
figure
hold on

% Forces vector
moduleF = (Forcerx.^2+Forcezx.^2).^0.5;
contourf(ra,za,moduleF, linspace(0,2000,14))
caxis([0,2000])
h = colorbar; ylabel(h,'$f_m$ [N/m$^3$]', 'Interpreter', 'latex','FontSize', 14)
Forcerx(2:2:end,:) = NaN;
Forcerx(:,2:2:end) = NaN;
Forcezx(2:2:end,:) = NaN;
Forcezx(:,2:2:end) = NaN;
Forcerx(ra==0.055 & za==0) = NaN;
Forcezx(ra==0.055 & za==0) = NaN;
quiver(ra,za, Forcerx, Forcezx,3,'r')

% Equilibrium surface and contours
plot(fA(nrA,:),gA(nrA,:),'k','linewidth',1)
plot([0,0,1,1]*max(fA(nrA,:)),[gA(nrA,1),0,0,gA(nrA,end)],'k-','linewidth',1)

% Format
axis equal
xlabel('r [m]')
ylabel('z [m]')

figure
mesh(ra,za,moduleF)
xlabel('r [m]')
ylabel('z [m]')
grid on

%% Viscosity
figure
hold on

% Forces vector
contourf(ra,za,etaA, 100,'LineStyle','none')
h = colorbar; ylabel(h,'$\eta$ [Pa$\cdot$s]', 'Interpreter', 'latex','FontSize', 14)

% Equilibrium surface and contours
plot(fA(nrA,:),gA(nrA,:),'k','linewidth',1)
plot([0,0,1,1]*max(fA(nrA,:)),[gA(nrA,1),0,0,gA(nrA,end)],'k-','linewidth',1)

% Format
axis equal
xlabel('r [m]')
ylabel('z [m]')

%% Pressure
figure, hold on
contourf(ra,za,pA, 40,'LineStyle','-')
h = colorbar; ylabel(h,'$P$ [Pa]', 'Interpreter', 'latex','FontSize', 14)

% Equilibrium surface and contours
plot(fA(nrA,:),gA(nrA,:),'k','linewidth',1)
plot([0,0.055,0.055], [0,0,gA(nrA,end)],'k','linewidth',1)

% Format
axis equal
xlabel('r [m]')
ylabel('z [m]')


%% Theoretical pressure distribution
xiP = 0.1;
mu0 = 4*pi*1e-7;
figure, hold on
contourf(ra,za,HA,0:10000:2e5,'LineStyle','-')
h = colorbar; ylabel(h,'$H$ [A/m]', 'Interpreter', 'latex','FontSize', 14)

% Equilibrium surface and contours
plot(fA(nrA,:),gA(nrA,:),'k','linewidth',1)
plot([0,0,1,1]*max(fA(nrA,:)),[gA(nrA,1),0,0,gA(nrA,end)],'k-','linewidth',1)

% Format
axis equal
xlabel('r [m]')
ylabel('z [m]')

%% Scalar potential phi in D
figure
subplot(2,3,1)
scatter(rD(:),zD(:),30,(phiD(:)),'filled')
title('$log(\phi)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar
subplot(2,3,2)
scatter(rD(:),zD(:),30,(dphidrD),'filled')
title('$log(d\phi/dr)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar
subplot(2,3,3)
scatter(rD(:),zD(:),30,(dphidrrD),'filled')
title('$log(d^2\phi/dr^2)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar
subplot(2,3,4)
scatter(rD(:),zD(:),30,(dphidzD),'filled')
title('$log(d\phi/dz)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar
subplot(2,3,5)
scatter(rD(:),zD(:),30,(dphidzzD),'filled')
title('$log(d^2\phi/dz^2)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar
subplot(2,3,6)
scatter(rD(:),zD(:),30,(dphidrzD),'filled')
title('$log(d^2\phi/drdz)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar

figure
subplot(2,3,1)
scatter(rD(:),zD(:),30,(psiD(:)),'filled')
title('$log(\psi)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar
subplot(2,3,2)
scatter(rD(:),zD(:),30,(dpsidrD),'filled')
title('$log(d\psi/dr)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar
subplot(2,3,3)
scatter(rD(:),zD(:),30,(dpsidrrD),'filled')
title('$log(d^2\psi/dr^2)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar
subplot(2,3,4)
scatter(rD(:),zD(:),30,(dpsidzD),'filled')
title('$log(d\psi/dz)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar
subplot(2,3,5)
scatter(rD(:),zD(:),30,(dpsidzzD),'filled')
title('$log(d^2\psi/dz^2)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar
subplot(2,3,6)
scatter(rD(:),zD(:),30,(dpsidrzD),'filled')
title('$log(d^2\psi/drdz)$')
xlabel('r (m)')
ylabel('z (m)')
axis equal
colorbar


%% Style
set(findobj('Type', 'Text'), 'Interpreter', 'latex');
set(findobj('Type', 'Legend'), 'Interpreter', 'latex');
set(findobj('Type', 'axes'), 'FontName','Times New Roman','FontSize', 14,...
    'XMinorGrid','off','YMinorGrid','off','TickLabelInterpreter', 'latex');
set(findobj('Type', 'ColorBar'), 'TickLabelInterpreter', 'latex');
set(findobj('Type', 'figure'), 'Units', 'centimeters', 'Position', [5 5 15 12]); % [0 0 16 12]
set(findall(findobj('Type', 'axes'), 'Type', 'Text'), 'Interpreter', 'latex');
set(findall(findobj('Type', 'axes'), 'Type', 'Line'), 'LineWidth', 1);

