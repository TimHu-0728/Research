%% Setting equations of Block A
%%%% xi coordinate
z1=sym('z1','real'); % z
%%%% eta coordinate
r1=sym('r1','real'); %inner and outer r  (block 1 and block 2)
%theta coordinate;
t1=sym('t1','real');
%time
time1=sym('time1','real')


NVAR = NVC;
NDER = NDC;
dimension = 2;

for k=1:NVAR  %variable
    for i=1:NDER  %derivatives
        %%
        % We create an array of symbolic variables, _yfA_
        yf1(k,i)=sym([ 'yf1','v',num2str(k),'d',num2str(i)],'real');
    end
    tco = list_var_E {k};
    [bas nto] = size(tco);
    list_variables {k} = tco(1:nto-1);
end
list_derivatives = list_der_E;

%%
% We assign now variables names to the created array. In the case of block
% 1 the name is created combining the name of the variables involved and
% name of the derivatives. Some examples,
%
% $$ \frac{\partial w}{\partial \xi} = dwdz0 $$
%
% or
%
% $$ \frac{\partial^2 p}{\partial \xi \partial \eta} = dpdzr0 $$
%

for k=1:NVAR
    eval([list_variables{k} list_derivatives{1} ' = yf1(k,1);']);
    for i=2:NDER
        eval(['d' list_variables{k}  list_derivatives{i} ' = yf1(k,i);']);
    end
end


%2D
for k=1:NVAR
    up=list_variables{k}
    name=['d',up,'0t=0;'];
    eval(name);
    name=['d',up,'0tt=0;'];
    eval(name);
    name=['d',up,'0zt=0;'];
    eval(name);
    name=['d',up,'0rt=0;'];
    eval(name);
    % df0tt=0;
    % dg0t=0;
    % dg0tt=0;
    % df0zt=0;
    % dg0zt=0;
end




%% Mapping 1 for liquid bridge
r=f;
z=zout-r1.*zout;

t=t1;
time=time1;



df0r_res=df0r;
df0rr_res=df0rr;
df0rz_res=df0rz;
df0rt_res=df0rt;
g=0*f;
dg0z=0*f;
dg0zz=0*f;
dg0zt=0*f;
dg0t=0*f;
dg0tt=0*f;
dg0r=0;
dg0rr=0;
dg0rz=0;
dg0rt=0;

df0r=0;
df0rr=0;
df0rz=0;
df0rt=0;

J(1,1)=diff(r,r1)+diff(r,f)*df0r+diff(r,g)*dg0r; %dr/dz1
J(1,2)=diff(r,z1)+diff(r,f)*df0z+diff(r,g)*dg0z;%dr/dr1
J(1,3)=diff(r,t1)+diff(r,f)*df0t+diff(r,g)*dg0t;%dr/dt1

J(2,1)=diff(z,r1)+diff(z,f)*df0r+diff(z,g)*dg0r; %dz/dz1
J(2,2)=diff(z,z1)+diff(z,f)*df0z+diff(z,g)*dg0z;%dz/dr1
J(2,3)=diff(z,t1)+diff(z,f)*df0t+diff(z,g)*dg0t;%dz/dt1

J(3,2)=diff(t,r1)+diff(t,f)*df0r+diff(t,g)*dg0r; %dt/dz1
J(3,1)=diff(t,z1)+diff(t,f)*df0z+diff(t,g)*dg0z;%dt/dr1
J(3,3)=diff(t,t1)+diff(t,f)*df0t+diff(t,g)*dg0t;%dt/dt1

%inverse Jacobian
J1n=inv(J);
A11=J1n(1,1)   %dr1/dr
A12=J1n(2,1)   %dz1/dr
A13=J1n(3,1)   %dt1/dr

A21=J1n(1,2)   %dr1/dz
A22=J1n(2,2)   %dz1/dz
A23=J1n(3,2)   %dt1/dz

A31=J1n(1,3)   %dr1/dt
A32=J1n(2,3)   %dz1/dt
A33=J1n(3,3)   %dt1/dt



A11r1=diff(A11,r1)+diff(A11,f)*df0r+diff(A11,g)*dg0r+diff(A11,df0r)*df0rr+diff(A11,dg0r)*dg0rr+diff(A11,df0z)*df0rz+diff(A11,dg0z)*dg0rz+diff(A11,df0t)*df0rt+diff(A11,dg0t)*dg0rt;
A11z1=diff(A11,z1)+diff(A11,f)*df0z+diff(A11,g)*dg0z+diff(A11,df0r)*df0rz+diff(A11,dg0r)*dg0rz+diff(A11,df0z)*df0zz+diff(A11,dg0z)*dg0zz+diff(A11,df0t)*df0zt+diff(A11,dg0t)*dg0zt;
A11t1=diff(A11,t1)+diff(A11,f)*df0t+diff(A11,g)*dg0t+diff(A11,df0r)*df0rt+diff(A11,dg0r)*dg0rt+diff(A11,df0z)*df0zt+diff(A11,dg0z)*dg0zt+diff(A11,df0t)*df0tt+diff(A11,dg0t)*dg0tt;

A12r1=diff(A12,r1)+diff(A12,f)*df0r+diff(A12,g)*dg0r+diff(A12,df0r)*df0rr+diff(A12,dg0r)*dg0rr+diff(A12,df0z)*df0rz+diff(A12,dg0z)*dg0rz+diff(A12,df0t)*df0rt+diff(A12,dg0t)*dg0rt;
A12z1=diff(A12,z1)+diff(A12,f)*df0z+diff(A12,g)*dg0z+diff(A12,df0r)*df0rz+diff(A12,dg0r)*dg0rz+diff(A12,df0z)*df0zz+diff(A12,dg0z)*dg0zz+diff(A12,df0t)*df0zt+diff(A12,dg0t)*dg0zt;
A12t1=diff(A12,t1)+diff(A12,f)*df0t+diff(A12,g)*dg0t+diff(A12,df0r)*df0rt+diff(A12,dg0r)*dg0rt+diff(A12,df0z)*df0zt+diff(A12,dg0z)*dg0zt+diff(A12,df0t)*df0tt+diff(A12,dg0t)*dg0tt;

A13r1=diff(A13,r1)+diff(A13,f)*df0r+diff(A13,g)*dg0r+diff(A13,df0r)*df0rr+diff(A13,dg0r)*dg0rr+diff(A13,df0z)*df0rz+diff(A13,dg0z)*dg0rz+diff(A13,df0t)*df0rt+diff(A13,dg0t)*dg0rt;
A13z1=diff(A13,z1)+diff(A13,f)*df0z+diff(A13,g)*dg0z+diff(A13,df0r)*df0rz+diff(A13,dg0r)*dg0rz+diff(A13,df0z)*df0zz+diff(A13,dg0z)*dg0zz+diff(A13,df0t)*df0zt+diff(A13,dg0t)*dg0zt;
A13t1=diff(A13,t1)+diff(A13,f)*df0t+diff(A13,g)*dg0t+diff(A13,df0r)*df0rt+diff(A13,dg0r)*dg0rt+diff(A13,df0z)*df0zt+diff(A13,dg0z)*dg0zt+diff(A13,df0t)*df0tt+diff(A13,dg0t)*dg0tt;

A11r=A11r1.*A11 + A11z1.*A12 + A11t1.*A13;
A12r=A12r1.*A11 + A12z1.*A12 + A12t1.*A13;
A13r=A13r1.*A11 + A13z1.*A12 + A13t1.*A13;


A11z=A11r1.*A21 + A11z1.*A22 + A11t1.*A23;
A12z=A12r1.*A21 + A12z1.*A22 + A12t1.*A23;
A13z=A13r1.*A21 + A13z1.*A22 + A13t1.*A23;



for k=1:NVAR
    up=list_variables{k};
    name=['d',up,'rr1=A11.*d',up,'0rr + A12.*d',up,'0rz + A13.*d',up,'0rt;'];
    eval(name)
    name=['d',up,'rz1=A11.*d',up,'0rz + A12.*d',up,'0zz + A13.*d',up,'0zt;'];
    eval(name)
    name=['d',up,'rt1=A11.*d',up,'0rt + A12.*d',up,'0zt + A13.*d',up,'0tt;'];
    eval(name)
    name=['d',up,'rr=A11r.*d',up,'0r + A11.*d',up,'rr1 + A12r.*d',up,'0z + A12.*d',up,'rz1 + A13r.*d',up,'0t + A13.*d',up,'rt1;'];
    eval(name)
end
% durr1=A11.*du0rr + A12.*du0rz + A13.*du0rt;
% durz1=A11.*du0rz + A12.*du0zz + A13.*du0zt;
% durt1=A11.*du0rt + A12.*du0zt + A13.*du0tt;
%durr=A11r.*du0r + A11.*durr1 + A12r.*du0z + A12.*durz1 + A13r.*du0t + A13.*durt1;


A21r1=diff(A21,r1)+diff(A21,f)*df0r+diff(A21,g)*dg0r+diff(A21,df0r)*df0rr+diff(A21,dg0r)*dg0rr+diff(A21,df0z)*df0rz+diff(A21,dg0z)*dg0rz+diff(A21,df0t)*df0rt+diff(A21,dg0t)*dg0rt;
A21z1=diff(A21,z1)+diff(A21,f)*df0z+diff(A21,g)*dg0z+diff(A21,df0r)*df0rz+diff(A21,dg0r)*dg0rz+diff(A21,df0z)*df0zz+diff(A21,dg0z)*dg0zz+diff(A21,df0t)*df0zt+diff(A21,dg0t)*dg0zt;
A21t1=diff(A21,t1)+diff(A21,f)*df0t+diff(A21,g)*dg0t+diff(A21,df0r)*df0rt+diff(A21,dg0r)*dg0rt+diff(A21,df0z)*df0zt+diff(A21,dg0z)*dg0zt+diff(A21,df0t)*df0tt+diff(A21,dg0t)*dg0tt;

A22r1=diff(A22,r1)+diff(A22,f)*df0r+diff(A22,g)*dg0r+diff(A22,df0r)*df0rr+diff(A22,dg0r)*dg0rr+diff(A22,df0z)*df0rz+diff(A22,dg0z)*dg0rz+diff(A22,df0t)*df0rt+diff(A22,dg0t)*dg0rt;
A22z1=diff(A22,z1)+diff(A22,f)*df0z+diff(A22,g)*dg0z+diff(A22,df0r)*df0rz+diff(A22,dg0r)*dg0rz+diff(A22,df0z)*df0zz+diff(A22,dg0z)*dg0zz+diff(A22,df0t)*df0zt+diff(A22,dg0t)*dg0zt;
A22t1=diff(A22,t1)+diff(A22,f)*df0t+diff(A22,g)*dg0t+diff(A22,df0r)*df0rt+diff(A22,dg0r)*dg0rt+diff(A22,df0z)*df0zt+diff(A22,dg0z)*dg0zt+diff(A22,df0t)*df0tt+diff(A22,dg0t)*dg0tt;

A23r1=diff(A23,r1)+diff(A23,f)*df0r+diff(A23,g)*dg0r+diff(A23,df0r)*df0rr+diff(A23,dg0r)*dg0rr+diff(A23,df0z)*df0rz+diff(A23,dg0z)*dg0rz+diff(A23,df0t)*df0rt+diff(A23,dg0t)*dg0rt;
A23z1=diff(A23,z1)+diff(A23,f)*df0z+diff(A23,g)*dg0z+diff(A23,df0r)*df0rz+diff(A23,dg0r)*dg0rz+diff(A23,df0z)*df0zz+diff(A23,dg0z)*dg0zz+diff(A23,df0t)*df0zt+diff(A23,dg0t)*dg0zt;
A23t1=diff(A23,t1)+diff(A23,f)*df0t+diff(A23,g)*dg0t+diff(A23,df0r)*df0rt+diff(A23,dg0r)*dg0rt+diff(A23,df0z)*df0zt+diff(A23,dg0z)*dg0zt+diff(A23,df0t)*df0tt+diff(A23,dg0t)*dg0tt;

A21z=A21r1.*A21 + A21z1.*A22 + A21t1.*A23;
A22z=A22r1.*A21 + A22z1.*A22 + A22t1.*A23;
A23z=A23r1.*A21 + A23z1.*A22 + A23t1.*A23;




for k=1:NVAR
    up=list_variables{k};
    name=['d',up,'zr1=A21.*d',up,'0rr + A22.*d',up,'0rz + A23.*d',up,'0rt;'];
    eval(name)
    name=['d',up,'zz1=A21.*d',up,'0rz + A22.*d',up,'0zz + A23.*d',up,'0zt;'];
    eval(name)
    name=['d',up,'zt1=A21.*d',up,'0rt + A22.*d',up,'0zt + A23.*d',up,'0tt;'];
    eval(name)
    name=['d',up,'zz=A21z.*d',up,'0r + A21.*d',up,'zr1 + A22z.*d',up,'0z + A22.*d',up,'zz1 + A23z.*d',up,'0t + A23.*d',up,'zt1'];
    eval(name)
end

% duzr1=A21.*du0rr + A22.*du0rz + A23.*du0rt;
% duzz1=A21.*du0rz + A22.*du0zz + A23.*du0zt;
% duzt1=A21.*du0rt + A22.*du0zt + A23.*du0tt;
% duzz=A21z.*du0r + A21.*duzr1 + A22z.*du0z + A22.*duzz1 + A23z.*du0t + A23.*duzt1;


A31r1=diff(A31,r1)+diff(A31,f)*df0r+diff(A31,g)*dg0r+diff(A31,df0r)*df0rr+diff(A31,dg0r)*dg0rr+diff(A31,df0z)*df0rz+diff(A31,dg0z)*dg0rz+diff(A31,df0t)*df0rt+diff(A31,dg0t)*dg0rt;
A31z1=diff(A31,z1)+diff(A31,f)*df0z+diff(A31,g)*dg0z+diff(A31,df0r)*df0rz+diff(A31,dg0r)*dg0rz+diff(A31,df0z)*df0zz+diff(A31,dg0z)*dg0zz+diff(A31,df0t)*df0zt+diff(A31,dg0t)*dg0zt;
A31t1=diff(A31,t1)+diff(A31,f)*df0t+diff(A31,g)*dg0t+diff(A31,df0r)*df0rt+diff(A31,dg0r)*dg0rt+diff(A31,df0z)*df0zt+diff(A31,dg0z)*dg0zt+diff(A31,df0t)*df0tt+diff(A31,dg0t)*dg0tt;

A32r1=diff(A32,r1)+diff(A32,f)*df0r+diff(A32,g)*dg0r+diff(A32,df0r)*df0rr+diff(A32,dg0r)*dg0rr+diff(A32,df0z)*df0rz+diff(A32,dg0z)*dg0rz+diff(A32,df0t)*df0rt+diff(A32,dg0t)*dg0rt;
A32z1=diff(A32,z1)+diff(A32,f)*df0z+diff(A32,g)*dg0z+diff(A32,df0r)*df0rz+diff(A32,dg0r)*dg0rz+diff(A32,df0z)*df0zz+diff(A32,dg0z)*dg0zz+diff(A32,df0t)*df0zt+diff(A32,dg0t)*dg0zt;
A32t1=diff(A32,t1)+diff(A32,f)*df0t+diff(A32,g)*dg0t+diff(A32,df0r)*df0rt+diff(A32,dg0r)*dg0rt+diff(A32,df0z)*df0zt+diff(A32,dg0z)*dg0zt+diff(A32,df0t)*df0tt+diff(A32,dg0t)*dg0tt;

A33r1=diff(A33,r1)+diff(A33,f)*df0r+diff(A33,g)*dg0r+diff(A33,df0r)*df0rr+diff(A33,dg0r)*dg0rr+diff(A33,df0z)*df0rz+diff(A33,dg0z)*dg0rz+diff(A33,df0t)*df0rt+diff(A33,dg0t)*dg0rt;
A33z1=diff(A33,z1)+diff(A33,f)*df0z+diff(A33,g)*dg0z+diff(A33,df0r)*df0rz+diff(A33,dg0r)*dg0rz+diff(A33,df0z)*df0zz+diff(A33,dg0z)*dg0zz+diff(A33,df0t)*df0zt+diff(A33,dg0t)*dg0zt;
A33t1=diff(A33,t1)+diff(A33,f)*df0t+diff(A33,g)*dg0t+diff(A33,df0r)*df0rt+diff(A33,dg0r)*dg0rt+diff(A33,df0z)*df0zt+diff(A33,dg0z)*dg0zt+diff(A33,df0t)*df0tt+diff(A33,dg0t)*dg0tt;

A31t=A31r1.*A31 + A31z1.*A32 + A31t1.*A33;
A32t=A32r1.*A31 + A32z1.*A32 + A32t1.*A33;
A33t=A33r1.*A31 + A33z1.*A32 + A33t1.*A33;

for k=1:NVAR
    up=list_variables{k};
    name=['d',up,'tr1=A31.*d',up,'0rr + A32.*d',up,'0rz + A33.*d',up,'0rt;']
    eval(name)
    name=['d',up,'tz1=A31.*d',up,'0rz + A32.*d',up,'0zz + A33.*d',up,'0zt;']
    eval(name)
    name=['d',up,'tt1=A31.*d',up,'0rt + A32.*d',up,'0zt + A33.*d',up,'0tt;']
    eval(name)
    name=['d',up,'tt=A31t.*d',up,'0r + A31.*d',up,'tr1 + A32t.*d',up,'0z + A32.*d',up,'tz1 + A33t.*d',up,'0t + A33.*d',up,'tt1;']
    eval(name)
    name=['d',up,'r=A11.*d',up,'0r+A12.*d',up,'0z+A13.*d',up,'0t; ']
    eval(name)
    name=['d',up,'z=A21.*d',up,'0r+A22.*d',up,'0z+A23.*d',up,'0t;']
    eval(name)
    name=['d',up,'t=A31.*d',up,'0r+A32.*d',up,'0z+A33.*d',up,'0t;']
    eval(name)
end


% dutr1=A31.*du0rr + A32.*du0rz + A33.*du0rt;
% dutz1=A31.*du0rz + A32.*du0zz + A33.*du0zt;
% dutt1=A31.*du0rt + A32.*du0zt + A33.*du0tt;
% dutt=A31t.*du0r + A31.*dutr1 + A32t.*du0z + A32.*dutz1 + A33t.*du0t + A33.*dutt1;
% dur=A11.*du0r+A12.*du0z+A13.*du0t;
% duz=A21.*du0r+A22.*du0z+A23.*du0t;
% dut=A31.*du0r+A32.*du0z+A33.*du0t;



%2D
%duz=A21.*du0r
%dur=A11.*du0r+A12.*du0z
%dphirz0=A11.*A21*dphi0rr+A12.*A21*dphi0rz+A11z*dphi0r+A12z*dphi0z;

for k=1:NVAR
    up=list_variables{k};
    name=['d',up,'rz=A11.*A21*d',up,'0rr+A12.*A21*d',up,'0rz+A11z*d',up,'0r+A12z*d',up,'0z;']
    eval(name)
end
%mesh
r=f;
z=zout+r1.*(-zout);


%computing scalar C given the curve M(H)
mu0=4*pi*1e-7;

Hr=(-dphiz/r)-dpsir;
dHrdr=(-dphirz/r+dphiz/r^2)-dpsirr;
dHrdz=-dphizz/r-dpsirz;
Hz   = dphir/r-dpsiz;
dHzdr= (dphirr/r-dphir/r^2)-dpsirz;
dHzdz= dphirz/r-dpsizz;

%computing vector B
Br=mu0*Hr;
Bz=mu0*Hz;
dBzdz=mu0*dHzdz;
dBrdr=mu0*dHrdr;

Hra= 0*Hr;
Hza= dphirr - dpsiz;

%% potential
syms phi_ref
FEE(1)=dBzdz+dBrdr+Br/r;
FEE(3)=phi-phi_ref;dHrdz-dHzdr;
%Equation for f
FEE(2)=df0rr_res;

%BClinea12++++
%axis (m=0)

FEElinea2(1)=phi; %v;
FEElinea2(2)=f;%df0r;
FEElinea2(3)=dpsir;

z1f    = z - mz1;
z2f    = z - mz2;
r1f    = sqrt(r.^2 + z1f.^2);
r2f    = sqrt(r.^2 + z2f.^2);

%potential
phim=-1/4/pi * ((m1*z1f)./r1f.^3 + (m2*z2f)./r2f.^3);

% Hr (2 is coil -c- or magnet -m-)
Hrm_c    = m1/4/pi * ((3*z1f.*r)./r1f.^5);
Hrm_m    = m1/4/pi * ((3*z1f.*r)./r1f.^5)        + m2/4/pi * ((3*z2f.*r)./r2f.^5);

% Hz (2 is coil -c- or magnet -m-)
Hzm_c    = m1/4/pi * ((2*z1f.^2 - r.^2)./r1f.^5);
Hzm_m    = m1/4/pi * ((2*z1f.^2 - r.^2)./r1f.^5) + m2/4/pi * ((2*z2f.^2 - r.^2)./r2f.^5);

% Phi_dipole (2 is coil -c- or magnet -m-)
% phid_c   = m2/4/pi * r^2/r2f.^3;
phid_m   = 0;

% Exact phi caused by a circular loop (18/04/2021)
syms phid_c

% Psi dipole (2 is coil -c- or magnet -m-)
psi_c    = m1/4/pi * z1f/r1f.^3;
psi_m    = m1/4/pi * z1f/r1f.^3 + m2/4/pi * z2f/r2f.^3;

% bottom
FEElinea1_c(1)=psi-psi_c; % -phi is in DE  (variable 1)
FEElinea1_c(2)=df0r_res;
FEElinea1_c(3)=phi-phid_c; % -phi is in DE  (variable 1)

FEElinea1_m(1)=psi-psi_m; % -phi is in DE  (variable 1)
FEElinea1_m(2)=df0r_res;
FEElinea1_m(3)=phi-phid_m; % -phi is in DE  (variable 1)

%conexion ADE
FVertexEAD(1)=phi;   %DDAE
FVertexEAD(2)=f-R;
FVertexEAD(3)=psi;   %DDAE

FVertexDEAE(1)=dphir;    %DEA
FVertexDEAE(2)=dpsir;
FVertexDEAE(3)=r;
FVertexDEAE(4)=z;
FVertexDEAE(5)=0;

%conexion with A
FEElinea9(1)=-Hr;%
FEElinea9(2)=-f;%
FEElinea9(3)=-Bz; %-Br is in Block DD (line 9) (variable 1)

FAElinea9(4)=-phi;%
FAElinea9(9)=-psi; %-Br is in Block DD (line 9) (variable 1)
FAElinea9(11)=0;

%right wall conexion with D
FEElinea19(1)=phi; % -phi is in DE  (variable 1)
FEElinea19(2)=f-R;
FEElinea19(3)=psi; % -psi is in DE  (variable 1)

FDElinea19(1)=dphir;%
FDElinea19(2)=dpsir; %-Br is in Block DD (line 9) (variable 1)
FDElinea19(3)=r;
FDElinea19(4)=z;
FDElinea19(5)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=reshape(yf1',NVAR *NDER,1); %variables

%x=[yf1(5,1);yf1(6,1)]


%jacobians
for k=1:NVAR
    
    dFEE(k,:)=jacobian(FEE(k),x);
    dFEElinea1_c(k,:)=jacobian(FEElinea1_c(k),x);
    dFEElinea1_m(k,:)=jacobian(FEElinea1_m(k),x);
    dFEElinea2(k,:)=jacobian(FEElinea2(k),x);
    dFEElinea9(k,:)=jacobian(FEElinea9(k),x);
    dFEElinea19(k,:)=jacobian(FEElinea19(k),x);
    dFVertexEAD(k,:)=jacobian(FVertexEAD(k),x);
end

matlabFunction(dphir,dphiz,dphirr,dphizz,dphirz,dpsir,dpsiz,dpsirr,dpsizz,dpsirz,...
    'file',[path_jacobian 'potentialsEE.m'],'vars',{z1,r1,x,pa});
matlabFunction(dHrdz-dHzdr,'file',  [path_jacobian 'curlH_E.m'],'vars',{z1,r1,x,pa});
matlabFunction(Hra,Hza,'file',  [path_jacobian 'MagneticfieldEEa.m'],'vars',{z1,r1,x,pa});

matlabFunction(Hr,Hz,'file',  [path_jacobian 'MagneticfieldEE.m'],'vars',{z1,r1,x,pa});
matlabFunction(FEE,dFEE,'file',  [path_jacobian 'equationFEE.m'],'vars',{z1,r1,x,pa,phi_ref});
matlabFunction(FEElinea1_c,dFEElinea1_c,'file',[path_jacobian 'equationFEElinea1_c.m'],'vars',{z1,r1,x,pa,phid_c});
matlabFunction(FEElinea1_m,dFEElinea1_m,'file',[path_jacobian 'equationFEElinea1_m.m'],'vars',{z1,r1,x,pa});
matlabFunction(FEElinea2,dFEElinea2,'file',[path_jacobian 'equationFEElinea2.m'],'vars',{z1,r1,x,pa});
matlabFunction(FEElinea9,dFEElinea9,'file',[path_jacobian 'equationFEElinea9.m'],'vars',{z1,r1,x,pa});
matlabFunction(FEElinea19,dFEElinea19,'file',[path_jacobian 'equationFEElinea19.m'],'vars',{z1,r1,x,pa});
matlabFunction(FVertexEAD,dFVertexEAD,'file',[path_jacobian 'equationVertexEAD.m'],'vars',{z1,r1,x,pa});

for k=1:NVA
    dFAElinea9(k,:)=jacobian(FAElinea9(k),x);
end
matlabFunction(FAElinea9,dFAElinea9,'file',[path_jacobian 'equationFAElinea9.m'],'vars',{z1,r1,x,pa});
%

for k=1:NVD
    dFDElinea19(k,:)=jacobian(FDElinea19(k),x);
    dFVertexDEAE(k,:)=jacobian(FVertexDEAE(k),x);
end
matlabFunction(FDElinea19,dFDElinea19,'file',[path_jacobian 'equationFDElinea19.m'],'vars',{z1,r1,x,pa});
matlabFunction(FVertexDEAE,dFVertexDEAE,'file',[path_jacobian 'equationFVertexDEAE.m'],'vars',{z1,r1,x,pa});