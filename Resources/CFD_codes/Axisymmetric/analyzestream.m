%% Plots stream function as a GIF

% Path locations
path_jacobian  = './eigen/jacobians2D/';
path_jacobian1 = './subroutinesMatlab/';
path_mesher    = './meshgenerator/';
path_jacmesh   = './meshgenerator/jacobians2D';
path_cmap      = './customcolormap';


% Path operations
restoredefaultpath
addpath(path_jacobian)
addpath(path_jacobian1)
addpath(path_mesher)
addpath(path_jacmesh)
addpath(path_cmap)

% Colormap
J     = customcolormap([0, 0.5, 1], [0, 0, 0.5156; 1,1,1; 0.5,0,0]);

% List files
files = dir('./soluciones/STREAM_v2-DYT-29-Apr-2021_Ib4000_zeta1,2571e-05tau5,5194e-07nrA101nzA101nrB41nrC21nrE61nzE101_T100.*.mat');

% Plot
figure(1);
axis tight manual % this ensures that getframe() returns a consistent size
for kind = 1:length(files)
    filename = files(kind).name;
    load(['./soluciones/',filename])
    time = str2double(filename(end-9:end-4));
    
    % Compute stream
    % % % x0=eval([order ']']);
    pA=zeros(nrA,nzA);
    order='[';
    for j=1:nbl
        bl = list_block{j};
        NVAR = eval(['NV' bl]);
        for i=1:NVAR
            lv = eval(['list_var_' bl '{' num2str(i) '};']);
            order = [order 'reshape(' lv ',nt' bl ',1);'];
        end
    end
    x0=eval([order ']']);
    x0m = x0;
    x0mm = x0m;
    
    xotoredeable
    r1A=repmat(r0A', [1 nzA]);
    ZA=(gA-gaxisA).*r1A;
    
    matrixAstream
    xotoredeable
    psi0=a\b; %MEJORAR EL TIEMPO DE INVERSÓN MENDIANTE MÉTODOS ITERATIVOS...
    psi0=reshape(psi0,nrA,nzA);
    
    plot(fA(nrA,:),gA(nrA,:),'k-','LineWidth',2)
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
    
    ns=24;
    for i=1:ns+1
        %va(i)=(i/15)^2*c1;
        vsa(i)=((i-1)/ns)*c1;
        vsb(i)=((i-1)/ns)*c2;
    end
    cmap = jet(2*(ns+1));
    
    % Record center position
    tsim(kind)  = time;
    zaxis(kind) = gA(nrA,1);
    
    % Plot
    %     contour(fA,ZA,psi3a ,vsa,'-')
    %     hold on
    %     contour(fA,ZA,psi3a,vsb,'-')
    colormap(J)
    contour(fA,ZA,psi3a ,[vsa,vsb])
    title([num2str(time) ' s'])
    colorbar
    xlabel('r [m]')
    ylabel('z [m]')
    ylim([0,0.06])
    xlim([0,0.055])
    caxis([-0.0025,0.0025])
    axis equal
    set(gca,'tickdir','out')
    set(gca,'box','off')
    set(gca,'FontSize',11)
    set(gca,'LineWidth',1.5)
    set(gca,'FontName','SansSerif')
    set(gca,'Layer','top')
    set(gcf,'Color',[1,1,1])
    hold off
    drawnow
    
    % Capture image
    frame = getframe(1);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the GIF File
    if kind == 1
        imwrite(imind,cm,'stream_05.gif','gif', 'Loopcount',inf,'DelayTime',0.25);
    elseif kind ~= length(files)
        imwrite(imind,cm,'stream_05.gif','gif','WriteMode','append','DelayTime',0.25);
    else
        imwrite(imind,cm,'stream_05.gif','gif','WriteMode','append','DelayTime',2);
    end
end

figure
plot(tsim, zaxis)
xlabel('t [s]')
ylabel('z [m]')
