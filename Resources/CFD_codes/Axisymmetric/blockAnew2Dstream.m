%% Setting equations of Block A
%%%% xi coordinate
z1=sym('z1','real'); % z
%%%% eta coordinate
r1=sym('r1','real'); %inner and outer r  (block 1 and block 2)
%theta coordinate;
t1=sym('t1','real');
%time
time1=sym('time1','real')


NVAR = NVA;
NDER = NDA;
dimension = 2;

for k=1:NVAR  %variable
    for i=1:NDER  %derivatives
        %%
        % We create an array of symbolic variables, _yfA_
        yf1(k,i)=sym([ 'yf1','v',num2str(k),'d',num2str(i)],'real');
    end
    tco = list_var_A {k};
    [bas nto] = size(tco);
    list_variables {k} = tco(1:nto-1);
end
list_derivatives = list_der_A;

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





%% Mapping 1 for sphere-case
% f10=sin(acos(1-g10));
% z=g.*z1;
% G1=sin(acos(1-g10*z1))-f10*z1;
% r=f.*z1+r1.*G1;
% t=t1;
% time=time1;

%% Mapping 1 for liquid bridge
z=r1.*(g-gaxis);   %%% z=g(xi,theta,time)
r=f; %%% r=eta*f(xi,theta,time)
t=t1;
time=time1;

dg0r_res=dg0r;
dg0rr_res=dg0rr;
dg0rz_res=dg0rz;
dg0rt_res=dg0rt;

df0r_res=df0r;
df0rr_res=df0rr;
df0rz_res=df0rz;
df0rt_res=df0rt;

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

 %temporal derivatives
 %mapping temporal  
 syms dz1t dr1t dt1t;
                  %%%% Sphere - case
%  eqn1 = (dg0r.*dr1t+dg0t.*dt1t+dgtime).*z1 + g.*dz1t==0;
%  eqn2 = (df0r.*dr1t+df0t.*dt1t+dftime).*z1 + f.*dz1t + G1.*dr1t+(diff(G1,z1).*dz1t+diff(G1,t1).*dt1t+diff(G1,r1).*dr1t+diff(G1,time1)).*r1==0;
%  eqn3 = dt1t==0;
 
                  %%%% liquid bridge - case
 eqn1 = df0z.*dz1t+df0t.*dt1t+df0time==0;
 eqn2 = (dg0z.*dz1t+dg0t.*dt1t+dg0time).*r1 + g.*dr1t==0;
 eqn3 = dt1t==0;
 [Aeqn,Beqn]=equationsToMatrix([eqn1,eqn2,eqn3],[dr1t,dz1t,dt1t]);
 Xeqn=linsolve(Aeqn,Beqn)
 
 dr1t=Xeqn(1);
 dz1t=Xeqn(2);
 dt1t=Xeqn(3);
 
 ut1=du0time+dr1t.*du0r+dz1t.*du0z+dt1t.*du0t;
 wt1=dw0time+dr1t.*dw0r+dz1t.*dw0z+dt1t.*dw0t;
 

  
  





%mesh
r=f;
z=r1.*(g-gaxis);












% %Equation for f
% 
% 
% 


% Forcer=subs(Forcer,Bz0,Bz);
% Forcer=subs(Forcer,Br0,Br);
% 
% Forcez=subs(Forcez,Bz0,Bz);
% Forcez=subs(Forcez,Br0,Br);

%% axial momentum
eta=duz-dwr;
u0=-dpz/r;
w0=+dpr/r;
  % eta=1
%% axial momentum
FAA(1)=(dpzz+dprr+dpr/r)-2*dpr./r+eta.*r;  %BC++++
%axis 
FAAlinea3(1)=p; %dwr;
%right wall conexion with vexte9
FAAlinea9r(1)=p;
FVertexABD(1)=p;
%wall
FAAlinea9b(1)=p;  
FVertexAD(1)=p;  
% Normal Balance
FAAlinea11(1)=-((u0-df0time).*dg0z-(w0-dg0time).*df0z);



 
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=reshape(yf1',NVAR *NDER,1); %variables

%x=[yf1(5,1);yf1(6,1)]


%jacobians
for k=1:1
    dFAA(k,:)=jacobian(FAA(k),x);
    dFAAlinea3(k,:)=jacobian(FAAlinea3(k),x);
    dFAAlinea11(k,:)=jacobian(FAAlinea11(k),x);
    dFAAlinea9b(k,:)=jacobian(FAAlinea9b(k),x);
    dFAAlinea9r(k,:)=jacobian(FAAlinea9r(k),x);
     dFVertexAD(k,:)=jacobian(FVertexAD(k),x);
     dFVertexABD(k,:)=jacobian(FVertexABD(k),x);
end


 matlabFunction(FVertexAD,dFVertexAD,'file',  [path_jacobian 'equationFVertexADstream.m'],'vars',{z1,r1,x,pa});
 matlabFunction(FVertexABD,dFVertexABD,'file',  [path_jacobian 'equationFVertexABDstream.m'],'vars',{z1,r1,x,pa});

 
 






  matlabFunction(FAA,dFAA,'file',  [path_jacobian 'equationFAAstream.m'],'vars',{z1,r1,x,pa});
% % 
%   % 
  matlabFunction(FAAlinea3,dFAAlinea3,'file',[path_jacobian 'equationFAAlinea3stream.m'],'vars',{z1,r1,x,pa});
 matlabFunction(FAAlinea11,dFAAlinea11,'file',[path_jacobian 'equationFAAlinea11stream.m'],'vars',{z1,r1,x,pa});

matlabFunction(FAAlinea9b,dFAAlinea9b,'file',[path_jacobian 'equationFAAlinea9bstream.m'],'vars',{z1,r1,x,pa});
matlabFunction(FAAlinea9r,dFAAlinea9r,'file',[path_jacobian 'equationFAAlinea9rstream.m'],'vars',{z1,r1,x,pa});


