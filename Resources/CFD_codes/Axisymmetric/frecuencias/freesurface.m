%% Analyze and compare free surface shape

% Parameters
rootdir = 'D:\Documents\Research\PhD\Ferrofluids\Viscous CFD\magentotalpsif2Dpotential3Dmass2final\soluciones\';
filename= @(I) ['DYT-05-Nov-2020_Ib',num2str(I*200),'zeta0tau0nrA81nzA81nrB121nrC101nrE51nzE51.mat'];
I       = 0:20;

% Plot FS
figure
for i = 1:length(I)
    % Load data
    load([rootdir,filename(I(i))],'fA','gA','nrA','nzA','hnominal','gaxisA')
    hdif = -gaxisA(1,1);
    
    % Equilibrium surface
    subplot(1,4,1)
    plot(fA(nrA,:),gA(nzA,:)+hdif,'linewidth',0.5)
    hold on, grid on
    
    % Center
    subplot(1,4,2)
    plot(I(i),gA(nzA,1)+hdif,'x')
    hold on, grid on
    
    % Contour
    subplot(1,4,3)
    plot(I(i),gA(nzA,end)+hdif,'x')
    hold on, grid on
    
    % Contour-Center
    subplot(1,4,4)
    plot(I(i),gA(nzA,end)-gA(nzA,1),'x')
    hold on, grid on
    ylim([0.002,0.020])
end

% Vessel contour
subplot(1,4,1)
plot([0,0.055,0.055], [0,0,0.06],'k','linewidth',1)
axis equal