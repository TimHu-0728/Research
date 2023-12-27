clc
clear
close all

% [User Input] Position of B field investigated (m) 
cartesian_coord=[0,0,0];              
x_user=cartesian_coord(1);
y_user=cartesian_coord(2);
z_user=cartesian_coord(3);

% [User Input] Coils positions [x y z] (m)
coil_positions=[0 0 0];          % first coil position;  
       

% [User Input] Coils orientation angles [psi (z-axis), theta (y-axis), phi (x-axis)] (deg), 
Euler_angles=deg2rad([0 45 0]);
                 
% [User Input] Rectangular Coil Parameters
a=3.175;                        % Rectangular coil length (m)
b=3.175;                        % Rectangular coil width (m)
h=3.175;                        % Total height (m)
N=50;                       % Number of stacks
I=1177.7e3*h/N;                    % Current (A)
mu=4*pi*1e-7;        % Vacuum magnetic permeability (N/(A^2))

%
% Output the magnetic field information
% [Bx_user,By_user,Bz_user]=B_net(x_user,y_user,z_user,a,b,N,I,h,mu,coil_positions,Euler_angles);
% fprintf('Rectangular Coil:\nLength a = %6.3f (m)\nWidth b = %6.3f (m)\nHeight h = %6.3f (m)\nNumber of Stack N = %i\nCurrent I = %6.4f (A)\n\nAt position\nr = [%6.4f %6.4f %6.4f ] (m)\n\nMagnetic Field B is\nB = [%6.4g %6.4g %6.4g ] (Tesla)\n',a,b,h,N,I,x_user,y_user,z_user,Bx_user,By_user,Bz_user)
% 

range      = 6;                          % [mm]
gridsz     = 10;
x_range    = linspace(-range,range,gridsz);
y_range    = linspace(-range,range,gridsz);
z_range    = linspace(-range,range,gridsz);
[x,y,z]    = meshgrid(x_range,y_range,z_range);


% % Set dimension range of the space
% xrange=max(coil_positions(:,1))+b/2+0.05;
% yrange=max(coil_positions(:,2))+a/2+0.05;
% zrange=max(coil_positions(:,3))+2*h+0.025;
% 
% % Set mesh size
% [x,y,z]=meshgrid(linspace(-xrange,xrange,20),linspace(-yrange,yrange,20),linspace(-zrange,zrange,15));

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


% Function calculates the B field around a rectangular loop
function B=B_rect(x,y,zp,a,b,mu,I)
    % Define sin(alpha) function for each segment
    sina12=(a/2-y)/sqrt((x-b/2)^2+(-zp)^2+(a/2-y)^2);
    sina11=-(a/2+y)/sqrt((x-b/2)^2+(-zp)^2+(a/2+y)^2);
    sina22=(b/2-x)/sqrt((-y-a/2)^2+(-zp)^2+(b/2-x)^2);
    sina21=-(b/2+x)/sqrt((-y-a/2)^2+(-zp)^2+(b/2+x)^2);
    sina32=(a/2+y)/sqrt((-x-b/2)^2+(-zp)^2+(a/2+y)^2);
    sina31=-(a/2-y)/sqrt((-x-b/2)^2+(-zp)^2+(a/2-y)^2);
    sina42=(b/2+x)/sqrt((y-a/2)^2+(-zp)^2+(b/2+x)^2);
    sina41=-(b/2-x)/sqrt((y-a/2)^2+(-zp)^2+(b/2-x)^2);
    % Define B field for each segment
    B1=[mu*I/4/pi/sqrt((x-b/2)^2+(-zp)^2)*(sina12-sina11)*(-sin(atan2(-zp,x-b/2)));
    mu*I/4/pi/sqrt((x-b/2)^2+(-zp)^2)*(sina12-sina11)*0;
    mu*I/4/pi/sqrt((x-b/2)^2+(-zp)^2)*(sina12-sina11)*(-cos(atan2(-zp,x-b/2)))];
    B2=[mu*I/4/pi/sqrt((-y-a/2)^2+(-zp)^2)*(sina22-sina21)*0;
    mu*I/4/pi/sqrt((-y-a/2)^2+(-zp)^2)*(sina22-sina21)*sin(atan2(-zp,-y-a/2));
    mu*I/4/pi/sqrt((-y-a/2)^2+(-zp)^2)*(sina22-sina21)*-cos(atan2(-zp,-y-a/2))];
    B3=[mu*I/4/pi/sqrt((-x-b/2)^2+(-zp)^2)*(sina32-sina31)*sin(atan2(-zp,-x-b/2));
    mu*I/4/pi/sqrt((-x-b/2)^2+(-zp)^2)*(sina32-sina31)*0;
    mu*I/4/pi/sqrt((-x-b/2)^2+(-zp)^2)*(sina32-sina31)*(-cos(atan2(-zp,-x-b/2)))];
    B4=[mu*I/4/pi/sqrt((y-a/2)^2+(-zp)^2)*(sina42-sina41)*0;
    mu*I/4/pi/sqrt((y-a/2)^2+(-zp)^2)*(sina42-sina41)*(-sin(atan2(-zp,y-a/2)));
    mu*I/4/pi/sqrt((y-a/2)^2+(-zp)^2)*(sina42-sina41)*(-cos(atan2(-zp,y-a/2)))];
    % Define B field for one stack of rectangular coil
    B=B1+B2+B3+B4;
end


% Function calculate B field around N stack rectangular coil 
function B=B_coil(x,y,zp,a,b,mu,I,N,h)
elvt=h/N;    
B=zeros(3,1);
    for i=1:N
        B=B+B_rect(x,y,zp+(-1)^i*floor(i/2)*elvt,a,b,mu,I);
    end
end





% Function calculate components of Magnetic Field Bx, By, Bx with input parameters:
% Position r = [x,y,z] (m),
% Coil length a (m), width b (m), total height h (m), number of stack N, current I (A), mu (N*A^(-2))
% Coil orientation angles, psi,theta,phi (rad)
%
function B=getB(x,y,z,a,b,N,I,h,mu,psi,theta,phi) 
% Calculate the position vector before the rotation r_old (m)
r_old= quatrotate(eul2quat(-[psi theta phi]),[x,y,z]);

% Evaluate the magnetic field components Bx, By, Bz at r_old
B_old=B_coil(r_old(1),r_old(2),r_old(3),a,b,mu,I,N,h);

% Rotate the magnetic field B evaluated at r_old to the new position
B=-quatrotate(eul2quat([psi theta phi]),B_old')';
end





% Function calculates the net magnetic field B_net with coils centered at different positions.
function [Bx_net,By_net,Bz_net]=B_net(x,y,z,a,b,N,I,h,mu,coil_positions,Euler_angles)
% Computations
    B=zeros(size(coil_positions,1),3);
% Calculate the magnetic field components Bx, By, Bz of each coil with shifted coil origin
    for i=1:size(coil_positions,1)
        B(i,:)=getB(x-coil_positions(i,1),y-coil_positions(i,2),z-coil_positions(i,3),a,b,N,I,h,mu,Euler_angles(i,1),Euler_angles(i,2),Euler_angles(i,3))';
    end
% Merge together to form the net magnetic field
    Bnet=sum(B,1);
    Bx_net=Bnet(1);
    By_net=Bnet(2);
    Bz_net=Bnet(3);
end

