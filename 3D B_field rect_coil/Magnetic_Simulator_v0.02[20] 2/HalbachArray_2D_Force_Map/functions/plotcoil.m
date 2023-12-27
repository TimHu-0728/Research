function plotcoil(a,b,h,N,x_coil,y_coil,z_coil,psi,theta,phi,color)
% Plotting rectangular coil shape with length a (m), width b (m), number of stack N, and total height h (m)
%
% plotcoil(a,b,h,N,x_coil,y_coil,z_coil,psi,theta,phi)
%
% INPUTS:
%
% OUTPUTS:
%

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

% Plot
plot3(r_new(:,1),r_new(:,2),r_new(:,3),'LineWidth',1.5,'Color',color)
hold on
plot3(x_coil,y_coil,z_coil,'+r')
axis equal
end