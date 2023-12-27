clc
clear
close all

% [User Input] Position of B field investigated (m) 
cartesian_coord=[2.8,2.4,0.5];              
x_user=cartesian_coord(1);
y_user=cartesian_coord(2);
z_user=cartesian_coord(3);

% [User Input] Coils positions [x y z] (m)
coil_positions=[4 0 0;          % first coil position
                -4 0 0;         % second coil position
                0 4 0;          % third coil position 
                0 -4 0;         % forth coil position
                0 0 4;          % fifth coil position 
                0 0 -4];        % sixth coil position 
       

% [User Input] Coils orientation angles [psi (z-axis), theta (y-axis), phi (x-axis)] (deg), 
Euler_angles=deg2rad([0 -90 0;       % first coil euler angles
                      0 90 0;      % second coil euler angles
                      0 0 90;     % third coil euler angles
                      0 0 -90;    % forth coil euler angles
                      0 0 0; 
                      0 0 180]);
                 
% [User Input] Rectangular Coil Parameters
a=3;                        % Rectangular coil length (m)
b=3;                        % Rectangular coil width (m)
h=1;                        % Total height (m)
N=30;                       % Number of stacks
I=60000;                    % Current (A)
mu=1.25663706212e-6;        % Vacuum magnetic permeability (N/(A^2))


%
% Output the magnetic field information
[Bx_user,By_user,Bz_user]=B_net(x_user,y_user,z_user,a,b,N,I,h,mu,coil_positions,Euler_angles);
fprintf('Rectangular Coil:\nLength a = %6.3f (m)\nWidth b = %6.3f (m)\nHeight h = %6.3f (m)\nNumber of Stack N = %i\nCurrent I = %6.4f (A)\n\nAt position\nr = [%6.3f %6.3f %6.3f ] (m)\n\nMagnetic Field B is\nB = [ %6.4g %6.4g %6.4g ] (Tesla)\n',a,b,h,N,I,x_user,y_user,z_user,Bx_user,By_user,Bz_user)



% Set dimension range of the space
xrange=max(coil_positions(:,1))+b/2+1;
yrange=max(coil_positions(:,2))+a/2+1;
zrange=max(coil_positions(:,3))+2*h+1;

% Set mesh size
[x,y,z]=meshgrid(linspace(-xrange,xrange,20),linspace(-yrange,yrange,20),linspace(-zrange,zrange,15));

% Calculate magnetic field at each mesh point
Bx=zeros(size(x));
By=zeros(size(y));
Bz=zeros(size(z));
for i=1:size(Bx,1)
    for j=1:size(Bx,2)
        for k=1:size(Bx,3)
            [Bx(i,j,k),By(i,j,k),Bz(i,j,k)]=B_net(x(i,j,k),y(i,j,k),z(i,j,k),a,b,N,I,h,mu,coil_positions,Euler_angles);
        end
    end
end

% Use quiver3 to plot B field in 3d space
scal=1.2;      % vector scale in quiver3
q=quiver3(x,y,z,Bx,By,Bz,scal,'LineWidth',0.75,'MaxHeadSize',1,'Color',[0.2 0.6 0.2]);
axis equal
xlabel('X (m)','FontSize',14)
ylabel('Y (m)','FontSize',14)
zlabel('Z (m)','FontSize',14)
title('Magnetic field of Rectangular coil (Tesla)','FontSize',16)
hold on

% plot 3d coil shape
netplot(a,b,h,N,coil_positions,Euler_angles)




%% Functions
% Plotting rectangular coil shape with length a (m), width b (m), number of stack N, and total height h (m)
function plotcoil(a,b,h,N,x_coil,y_coil,z_coil,psi,theta,phi)
slope=h/N/(2*a+2*b);
xc=-b/2*ones(4*N+1,1);
yc=-a/2*ones(4*N+1,1);
zc=zeros(4*N+1,1);
a_el=a*slope;
b_el=b*slope;

for i=1:4:(length(zc)-2)
    xc(i+2)=xc(i+2)*(-1);
    xc(i+3)=xc(i+3)*(-1);
    yc(i+1)=yc(i+1)*(-1);
    yc(i+2)=yc(i+2)*(-1);
end

for j=2:length(zc)
    if mod(j,2)==0
        zc(j)=zc(j)+a_el+zc(j-1);
    elseif mod(j,2)==1
        zc(j)=zc(j)+b_el+zc(j-1);
    end
end
zc=zc-h/2;

r_new=zeros(4*N+1,3);
for i=1:size(r_new,1)
    r_new(i,:)=quatrotate(eul2quat([psi theta phi]),[xc(i),yc(i),zc(i)])+[x_coil,y_coil,z_coil];
end

plot3(r_new(:,1),r_new(:,2),r_new(:,3),'LineWidth',1.5,'Color',[0 0 0])
hold on
plot3(x_coil,y_coil,z_coil,'+r')

end



% Plot all the coils in one graph
function netplot(a,b,h,N,coil_positions,Euler_angles)
    for i=1:size(coil_positions,1)
        plotcoil(a,b,h,N,coil_positions(i,1),coil_positions(i,2),coil_positions(i,3),Euler_angles(i,1),Euler_angles(i,2),Euler_angles(i,3));
        hold on
    end
    hold off
end



% Function calculate components of Magnetic Field Bx, By, Bx with input parameters:
% Position r = [x,y,z] (m),
% Coil length a (m), width b (m), total height h (m), number of stack N, current I (A), mu (N*A^(-2))
% Coil orientation angles, psi,theta,phi (rad)
%
function [Bx,By,Bz]=getB(x,y,z,a,b,N,I,h,mu,psi,theta,phi) 
%Define sin(alpha) function for each segment
sina12=@(x,y,zp,a,b) (a/2-y)./sqrt((x-b/2)^2+(-zp).^2+(a/2-y)^2);
sina11=@(x,y,zp,a,b) -(a/2+y)./sqrt((x-b/2)^2+(-zp).^2+(a/2+y)^2);
sina22=@(x,y,zp,a,b) (b/2-x)./sqrt((-y-a/2)^2+(-zp).^2+(b/2-x)^2);
sina21=@(x,y,zp,a,b) -(b/2+x)./sqrt((-y-a/2)^2+(-zp).^2+(b/2+x)^2);
sina32=@(x,y,zp,a,b) (a/2+y)./sqrt((-x-b/2)^2+(-zp).^2+(a/2+y)^2);
sina31=@(x,y,zp,a,b) -(a/2-y)./sqrt((-x-b/2)^2+(-zp).^2+(a/2-y)^2);
sina42=@(x,y,zp,a,b) (b/2+x)./sqrt((y-a/2)^2+(-zp).^2+(b/2+x)^2);
sina41=@(x,y,zp,a,b) -(b/2-x)./sqrt((y-a/2)^2+(-zp).^2+(b/2-x)^2);

%Define B field for each segment
B1={@(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((x-b/2)^2+(-zp).^2).*(sina12(x,y,zp,a,b)-sina11(x,y,zp,a,b)).*(-sin(atan2(-zp,x-b/2)));
    @(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((x-b/2)^2+(-zp).^2).*(sina12(x,y,zp,a,b)-sina11(x,y,zp,a,b)).*0;
    @(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((x-b/2)^2+(-zp).^2).*(sina12(x,y,zp,a,b)-sina11(x,y,zp,a,b)).*(-cos(atan2(-zp,x-b/2)))};

B2={@(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((-y-a/2)^2+(-zp).^2).*(sina22(x,y,zp,a,b)-sina21(x,y,zp,a,b)).*0;
    @(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((-y-a/2)^2+(-zp).^2).*(sina22(x,y,zp,a,b)-sina21(x,y,zp,a,b)).*sin(atan2(-zp,-y-a/2));
    @(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((-y-a/2)^2+(-zp).^2).*(sina22(x,y,zp,a,b)-sina21(x,y,zp,a,b)).*-cos(atan2(-zp,-y-a/2))};

B3={@(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((-x-b/2)^2+(-zp).^2).*(sina32(x,y,zp,a,b)-sina31(x,y,zp,a,b)).*sin(atan2(-zp,-x-b/2));
    @(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((-x-b/2)^2+(-zp).^2).*(sina32(x,y,zp,a,b)-sina31(x,y,zp,a,b)).*0;
    @(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((-x-b/2)^2+(-zp).^2).*(sina32(x,y,zp,a,b)-sina31(x,y,zp,a,b)).*(-cos(atan2(-zp,-x-b/2)))};

B4={@(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((y-a/2)^2+(-zp).^2).*(sina42(x,y,zp,a,b)-sina41(x,y,zp,a,b)).*0;
    @(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((y-a/2)^2+(-zp).^2).*(sina42(x,y,zp,a,b)-sina41(x,y,zp,a,b)).*(-sin(atan2(-zp,y-a/2)));
    @(x,y,zp,a,b,mu,I) mu*I/4/pi./sqrt((y-a/2)^2+(-zp).^2).*(sina42(x,y,zp,a,b)-sina41(x,y,zp,a,b)).*(-cos(atan2(-zp,y-a/2)))};

% Define B field for one stack of rectangular coil
B_rect={@(x,y,zp,a,b,mu,I) B1{1}(x,y,zp,a,b,mu,I)+B2{1}(x,y,zp,a,b,mu,I)+B3{1}(x,y,zp,a,b,mu,I)+B4{1}(x,y,zp,a,b,mu,I);
    @(x,y,zp,a,b,mu,I) B1{2}(x,y,zp,a,b,mu,I)+B2{2}(x,y,zp,a,b,mu,I)+B3{2}(x,y,zp,a,b,mu,I)+B4{2}(x,y,zp,a,b,mu,I);
    @(x,y,zp,a,b,mu,I) B1{3}(x,y,zp,a,b,mu,I)+B2{3}(x,y,zp,a,b,mu,I)+B3{3}(x,y,zp,a,b,mu,I)+B4{3}(x,y,zp,a,b,mu,I)};

% Define B field for N stacks of rectangular coil 
elvt=h/N;
Bx=@(x,y,zp,a,b,mu,I) 0;
By=@(x,y,zp,a,b,mu,I) 0;
Bz=@(x,y,zp,a,b,mu,I) 0;
for i=1:N
    Bx=@(x,y,zp,a,b,mu,I) Bx(x,y,zp,a,b,mu,I)+B_rect{1}(x,y,zp+(-1)^i*floor(i/2)*elvt,a,b,mu,I);
    By=@(x,y,zp,a,b,mu,I) By(x,y,zp,a,b,mu,I)+B_rect{2}(x,y,zp+(-1)^i*floor(i/2)*elvt,a,b,mu,I);
    Bz=@(x,y,zp,a,b,mu,I) Bz(x,y,zp,a,b,mu,I)+B_rect{3}(x,y,zp+(-1)^i*floor(i/2)*elvt,a,b,mu,I);
end

% Calculate the position vector before the rotation r_old (m)
r_old= quatrotate(eul2quat(-[psi theta phi]),[x,y,z]);

% Evaluate the magnetic field components Bx, By, Bz at r_old
Bx_old=Bx(r_old(1),r_old(2),r_old(3),a,b,mu,I);
By_old=By(r_old(1),r_old(2),r_old(3),a,b,mu,I);
Bz_old=Bz(r_old(1),r_old(2),r_old(3),a,b,mu,I);

% Rotate the magnetic field B evaluated at r_old to the new position
B=quatrotate(eul2quat([psi theta phi]),[Bx_old By_old Bz_old]);
Bx=-B(1);
By=-B(2);
Bz=-B(3);
end



% Function calculates the net magnetic field B_net with coils centered at different positions.
function [Bx_net,By_net,Bz_net]=B_net(x,y,z,a,b,N,I,h,mu,coil_positions,Euler_angles)
    Bx=zeros(size(coil_positions,1),1);
    By=zeros(size(coil_positions,1),1);
    Bz=zeros(size(coil_positions,1),1);
% Calculate the magnetic field components Bx, By, Bz of each coil with shifted coil origin
    for i=1:size(coil_positions,1)
        [Bx(i),By(i),Bz(i)]=getB(x-coil_positions(i,1),y-coil_positions(i,2),z-coil_positions(i,3),a,b,N,I,h,mu,Euler_angles(i,1),Euler_angles(i,2),Euler_angles(i,3));
    end
% Merge together to form the net magnetic field
    Bx_net=sum(Bx);
    By_net=sum(By);
    Bz_net=sum(Bz);
end
