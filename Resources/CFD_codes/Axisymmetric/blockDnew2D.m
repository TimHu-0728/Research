%% Setting equations of Block B
%

%% Creation of blocks
% The block 1 corresponds to the inner liquid fluid. At the main it has
% been defined that the number of main variables are NVAR which are derived
% in the equations pertaining to the fluid labelled as block 1 in
% NDER different ways. Note that in NDER is included the "zeroth order derivative" case.
% In this way a total of _NVAR_ x  _NDER_ variables are defined. Note that
% the derivative will be performed respect to the *mapped* variables, $\xi$ and
% $\eta$
%

NVAR = NVD;
NDER = NDD;
dimension = 2;


for k=1:NVAR  %variable
    for i=1:NDER  %derivatives
        %%
        % We create an array of symbolic variables, _yfB_
        yfD(k,i)=sym([ 'yfD','v',num2str(k),'d',num2str(i)],'real');
    end
    tco = list_var_D {k};
    [bas nto] = size(tco);
    list_variables {k} = tco(1:nto-1);
end
list_derivatives =  list_der_D;

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
    eval([list_variables{k} list_derivatives{1} ' = yfD(k,1);']);
for i=2:NDER    
eval(['d' list_variables{k} 'd' list_derivatives{i} ' = yfD(k,i);']);
end
end

%% 
% The index 0 indicates that the symbolic variables created corresponds to 
% mapped coordinates: $\xi$ (z-direction), $\eta$ (radial-direction) and 
% $\tau$ (time).

z0=sym('z0','real'); % xi=z0  
r0=sym('r0','real'); % eta=r0=r/f(z,t)
t0=sym('t0','real'); % tau = t0 = t

list_mapped = {'r0' 'z0'};

%%
% Auxiliaries variables of the mapping transformation depend in general of
% the mapped vartiables $\xi$, $\eta$ and $\tau$
%
%syms F(z0, r0, t0) 
%%
% In the present case only depends on $\xi$ and $\tau$
%

syms F0(z0,r0,t0) G0(z0,r0,t0)
r=F0;
z=G0;   
t=t0;

%%
% The geometrical jacobian is given by
%
% $$ J= \left(\begin{array}{cc} {\partial r}/{\partial \eta} & {\partial r}/{\partial \xi} \\ 
% {\partial z}/{\partial \eta} & {\partial z}/{\partial \xi}
% \end{array}\right) $$
%

J= jacobian([r z],[r0 z0]);

%%
% The inverse Jacobian will be needed,
%
% $$ J^{-1}= \left(\begin{array}{cc} {\partial \eta}/{\partial r} & {\partial \eta}/{\partial z} \\ 
% {\partial \xi}/{\partial r} & {\partial \xi}/{\partial z}
% \end{array}\right) $$
%

Jinv=inv(J);
Jinv=Jinv(z0,r0,t0);

%%
% It is important to remove the functional character of $J^{-1}$ by
% evaluating the tensor function, 
%

% Jinv=Jinv(z0, r0, t0); 

%%
% In our case since _F_ only depends on _z0_ and _t0_

% Jinv = Jinv(z0,t0)

%%
% We create also an array for the second order derivative,
%
% $$ J2 = \left(\begin{array}{cc} {\partial^2 \eta}/{\partial r^2} & {\partial^2 \eta}/{\partial z^2} \\ 
% {\partial^2 \xi}/{\partial r^2} & {\partial^2 \xi}/{\partial z^2}
% \end{array}\right) $$
%
%
% Since each term of the inverse Jacobian are functions of the mapped
% variables. See, for example, the 11 term,
% 
% $$ J^{-1}_{11} (\xi, \eta) = \frac{\partial \eta}{\partial r}   $$
% 
% The second order derivative can be therefore evaluated with the chain
% rule,
%
% $$ J2_{11} = \frac{\partial^2 \eta}{\partial r^2} = \frac{\partial J^{-1}_{11}}{\partial r}  
% = \frac{\partial J^{-1}_{11}}{\partial \eta} \frac{\partial \eta}{\partial r}  +
% \frac{\partial J^{-1}_{11}}{\partial \xi} \frac{\partial \xi}{\partial r} = 
% \frac{\partial J^{-1}_{11}}{\partial \eta} J^{-1}_{11}  +
% \frac{\partial J^{-1}_{11}}{\partial \xi}J^{-1}_{21} $$
%

for i=1:dimension
    for j=1:dimension     
        name =['J2(i,j) ='];
        for k=1:dimension
            name = [ name ' + Jinv(' num2str(k) ',j)* diff( Jinv(i,j) ,'  list_mapped{k} ')' ];
        end
        eval([name ';']);
    end
end

 name1 =['J22r =diff( Jinv(2,2),r0)*Jinv(1,1)+diff( Jinv(2,2),z0)*Jinv(2,1)'];
 eval([name1 ';']);
 name2 =['J12r =diff( Jinv(1,2),r0)*Jinv(1,1)+diff( Jinv(1,2),z0)*Jinv(2,1)'];
   eval([name2 ';']);

%%
% We have all the elements necessary for the geometrical mapping. We need
% to consider their consequences in the time transformation. In the 2D case 
% we are considering, the mapping is,
%
% $$ z = z(\xi, \eta, \tau) \quad  r = r(\xi, \eta, \tau) $$
% 
% Deriving the above expression respect to the (physical) time t we obtain,
%
% $$ \frac{d z}{d t} =  \frac{\partial z}{\partial \tau}  \frac{\partial
% \tau}{\partial t} + \frac{\partial z}{\partial \xi}  \frac{\partial
% \xi}{\partial t} + \frac{\partial z}{\partial \eta}  \frac{\partial
% \eta}{\partial t} $$
%
% Note that $z$ and $t$ are independent variables and in this case $\tau = t$. 
% Therefore, the above equation can be written as 
%
% $$  \frac{\partial z}{\partial \tau} + \frac{\partial z}{\partial \xi}  \frac{\partial
% \xi}{\partial t} + \frac{\partial z}{\partial \eta}  \frac{\partial
% \eta}{\partial t}  = 0$$
%
% Similarly, for the radial coordinate $r$,
%
% $$  \frac{\partial r}{\partial \tau} + \frac{\partial r}{\partial \xi}  \frac{\partial
% \xi}{\partial t} + \frac{\partial r}{\partial \eta}  \frac{\partial
% \eta}{\partial t}  = 0$$
%
% Note that the above equations form a system which allows to compute the
% unknowns
%
% $$ \frac{\partial \xi}{\partial t}, \quad \frac{\partial \eta}{\partial t}$
%
% Named symbolically _dr0dt_ and  _dz0dt_
%

 syms dr0dt  dz0dt ;
%          
  eqn1 = diff(r,t0) +diff(r,z0)*dz0dt +diff(r,r0)*dr0dt==0;
  eqn2 = diff(z,t0) +diff(z,z0)*dz0dt +diff(z,r0)*dr0dt==0;
%  
[Aeqn,Beqn]=equationsToMatrix([eqn1,eqn2],[dr0dt,dz0dt]);

 Xeqn=linsolve(Aeqn,Beqn);
 
%%
% Now it is necessary to perform a substitution to recover the orginal
% variables. It is important to keep the order of the substitutions; from
% order 0 to second order in this case, 
%

list_1 = {'F' 'dFdr0' 'dFdz0' 'dFdrr0'  'dFdzz0' 'dFdrz0' 'dFdt0' 'G' 'dGdr0' 'dGdz0' 'dGdrr0'  'dGdzz0' 'dGdrz0' 'dGdt0' };
F1=F0(z0, r0, t0);
dF1dr0=diff(F0(z0, r0, t0),r0);
dF1dz0=diff(F0(z0, r0, t0),z0);
dF1drr0=diff(F0(z0, r0, t0), r0, r0);
dF1dzz0=diff(F0(z0, r0, t0), z0, z0);
dF1dt0=diff(F0(z0, r0, t0), t0);
dF1drz0=diff(F0(z0, r0, t0), r0, z0);
G1=G0(z0, r0, t0);
dG1dr0=diff(G0(z0, r0, t0),r0);
dG1dz0=diff(G0(z0, r0, t0),z0);
dG1drr0=diff(G0(z0, r0, t0), r0, r0);
dG1dzz0=diff(G0(z0, r0, t0), z0, z0);
dG1dt0=diff(G0(z0, r0, t0), t0);
dG1drz0=diff(G0(z0, r0, t0), r0, z0);
list_2 = {'F1' 'dF1dr0' 'dF1dz0' 'dF1drr0'  'dF1dzz0' 'dF1drz0' 'dF1dt0' 'G1' 'dG1dr0' 'dG1dz0' 'dG1drr0'  'dG1dzz0' 'dG1drz0' 'dG1dt0' };


%list_2 = {'etar(z0, t0)', 'diff(etar(z0, t0), z0)', 'diff(etar(z0, t0), t0)', 'diff(etar(z0, t0), z0, z0)', 'etaz(z0, t0)' 'diff(etaz(z0, t0), z0)', 'diff(etaz(z0, t0), t0)', 'diff(etaz(z0, t0), z0, z0)'};

% list_1 = {'F' 'dFdr0' 'dFdz0' 'dFdrr0'  'dFdzz0' 'dFdrz0'}
% list_2 = {'F1' 'dF1dr0' 'dF1dz0' 'dF1drr0'  'dF1dzz0' 'dF1drz0'}


a = size(list_1);
Jinvrz=Jinv;
J2rz=J2;
J22rrz=J22r;
J12rrz=J12r;
for j=a(2):-1:1
    eval(['Jinv=subs(Jinv,' list_2{j} ',' list_1{j} ')']);
    eval(['J2=subs(J2,' list_2{j} ',' list_1{j} ')']);
    eval(['Xeqn=subs(Xeqn,' list_2{j} ',' list_1{j} ')']);
     eval(['Jinvrz=subs(Jinvrz,' list_2{j} ',' list_1{j} ')']);
    eval(['J2rz=subs(J2rz,' list_2{j} ',' list_1{j} ')']);
     eval(['J22rrz=subs(J22rrz,' list_2{j} ',' list_1{j} ')']);
     eval(['J12rrz=subs(J12rrz,' list_2{j} ',' list_1{j} ')']);
end

r=subs(r,F0,F);
z=subs(z,G0,G);

r=r(0,0,0);
z=z(0,0,0);

%%
% Now the correspondence between physical derivatives and mapped 
% derivatives has to be performed. For example for the first order 
% derivative of $r$ the relationship would be
%
% $$ \frac{\partial}{ \partial r} = \frac{\partial}{ \partial \xi} 
% \frac{\partial \xi}{ \partial r} + \frac{\partial}{ \partial \eta} 
% \frac{\partial \eta}{ \partial r} $$  
%
% and for the second order
% 
% $$ \frac{\partial^2}{ \partial r^2} = \frac{\partial^2}{ \partial^2 \xi} 
% \left( \frac{\partial \xi}{ \partial r} \right)^2 +
% \frac{\partial}{ \partial \xi} 
% \frac{\partial^2 \xi}{ \partial r^2} + \frac{\partial^2}{ \partial^2 \eta} 
% \left( \frac{\partial \eta}{ \partial r} \right)^2 +
% \frac{\partial}{ \partial \eta} 
% \frac{\partial^2 \eta}{ \partial r^2} 
% + 2 \frac{\partial^2}{ \partial \eta \partial \xi} \left( \frac{\partial \xi}{ \partial r} \right) \left( \frac{\partial \eta}{ \partial r} \right)  $$  
%
% First we set the name of the time derivative calculated above,
%

dr0dt=Xeqn(1);
dz0dt=Xeqn(2); 
eps=0;
%%
% *Pay attention to the order of the list below!* 
%
lder = {'r' 'z'};  %derivadas centradas
for k=1:NVAR
    var = list_variables{k}; 
    %%
    % zeroth order.
    eval([var ' =yfD(k,1);']);
    %%
    % time derivative    
    dvart = strcat('d', var,'dt = d', var, 'dt0');
    for j=1:dimension 
        dvart = [  dvart '+ d' var 'd' lder{j} '0*d' lder{j} '0dt'];
        %%
        % first order and second order spatial derivative
        dv1 = strcat('d', var, 'd',  lder{j},'=');
        dv2 = strcat('d', var, 'd', lder{j}, lder{j},'=');
        for i=1:dimension
            dv1 = [dv1 ' + d' var 'd' lder{i} '0*Jinv(' num2str(i) ',' num2str(j) ')' ];
            dv2 = [dv2 ' + d' var 'd' lder{i} '0*J2('   num2str(i) ',' num2str(j) ')' ];
            dv2 = [dv2 ' + d' var 'd' lder{i} lder{i} '0*Jinv(' num2str(i) ',' num2str(j) ')^2' ];
        end
        dv2 = [dv2 '+2*d' var 'drz0*Jinv(1,' num2str(j) ')*Jinv(2,' num2str(j) ');' ];
        eval([dvart ';'])
        eval([dv1 ';'])
        eval([dv2 ';'])
    end
    
      dv3=strcat('d', var, 'drz=');
    dv3=[dv3 '+ (d' var 'drr0*Jinvrz(1,1)+ d' var 'drz0*Jinvrz(2,1))* Jinvrz(1,2)'   '+ (d' var 'drz0*Jinvrz(1,1)+ d' var 'dzz0*Jinvrz(2,1))* Jinvrz(2,2)'] ;
    dv3=[dv3 '+ d' var 'dr0*J12rrz'  '+ d' var 'dz0*J22rrz'];
  eval([dv3 ';']) 
end




%computing scalar C given the curve M(H)
mu0=4*pi*1e-7;






Hr=-dphidz/r-dpsidr;
Hz= dphidr/r-dpsidz;
dHrdr=-dphidrz/r+dphidz/r^2-dpsidrr;
dHrdz=-dphidzz/r-dpsidrz;

dHzdr= dphidrr/r-dphidr/r^2-dpsidrz;
dHzdz= dphidrz/r-dpsidzz;

Br=mu0*Hr;
Bz=mu0*(Hz+Mz);
dBrdr=mu0*dHrdr;
dBrdz=mu0*dHrdz;
dBzdr=mu0*(dHzdr+dMzdr);
dBzdz=mu0*(dHzdz+dMzdz);

Hra= 0*Hr;
Hza= dphidrr - dpsidz;


 g11=dGdz0^2+dFdz0^2;
g22=dGdr0^2+dFdr0^2;
g12=dGdr0*dGdz0+dFdr0*dFdz0;
J=dGdr0*dFdz0-dGdz0*dFdr0;
  
  
  eps1=0;
  D1=eps1*((dFdz0^2+dGdz0^2)/(dFdr0^2+dGdr0^2))^0.5+(1-eps1);
  dD1dr0=diff(D1,dGdz0)*dGdrz0+diff(D1,dGdr0)*dGdrr0+diff(D1,dFdz0)*dFdrz0+diff(D1,dFdr0)*dFdrr0;
  dD1dz0=diff(D1,dGdz0)*dGdzz0+diff(D1,dGdr0)*dGdrz0+diff(D1,dFdz0)*dFdzz0+diff(D1,dFdr0)*dFdrz0;
  Q1=-(dD1dr0*dFdz0-dD1dz0*dFdr0)*J/D1;


%% potential
FDD(1)=dBzdz+dBrdr+Br/r;
FDD(2)=dHrdz-dHzdr;
FDD(3)=F-Fcoil;g22*dFdzz0+g11*dFdrr0-2*g12*dFdrz0;  
FDD(4)=g22*dGdzz0+g11*dGdrr0-2*g12*dGdrz0-Q1;
FDD(5)= Mz;



Jf = Ib/Sb; %current
FDDCoil(1)=dBzdz+dBrdr+Br/r;
FDDCoil(2)=dHrdz-dHzdr-Jf;
FDDCoil(3)=F-Fcoil;  
FDDCoil(4)=G-Gcoil;
FDDCoil(5)= Mz-Mmagnet;

FDDCoilContour(1)=dBzdz+dBrdr+Br/r;
FDDCoilContour(2)=dHrdz-dHzdr-Jf/2;
FDDCoilContour(3)=F-Fcoil;  
FDDCoilContour(4)=G-Gcoil;
FDDCoilContour(5)= Mz-Mmagnet/2;

FDDCoil1(1)=dBzdz+dBrdr+Br/r;
FDDCoil1(2)=dHrdz-dHzdr;
FDDCoil1(3)=F-Fcoil;  
FDDCoil1(4)=G-Gcoil;
FDDCoil1(5)= Mz;

FDDRightCoil(1)=dBzdz+dBrdr+Br/r;
FDDRightCoil(2)=dHrdz-dHzdr;
FDDRightCoil(3)=F-Fcoil;  
FDDRightCoil(4)=G-Gcoil;
FDDRightCoil(5)= Mz;

 %H=gradient phim
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


%conexion wih E+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FDDlinea19(1)=-dphidr;  % conecting with EDD 
FDDlinea19(2)=-dpsidr;
FDDlinea19(3)=-F;     % conecting with EDD -F
FDDlinea19(4)=-G;     % conecting with EDD -G
FDDlinea19(5)=Mz;

%E
FEDlinea19(1)=-phi;  % conecting with EDD 
FEDlinea19(2)=0; 
FEDlinea19(3)=-psi;
    % conecting with EDD -F


%Vertex DEB  psiD=psiE dpsiDr=dpsiEeR%phiA=psiD+++++++++++++++++++++++++
FVertexDEA(1)=-dphidr;  % conecting with E
FVertexDEA(2)=-dpsidr;
FVertexDEA(3)=-F;    
FVertexDEA(4)=-G;
FVertexDEA(5)=Mz;
%E
FVertexEADD(1)=-phi;  
FVertexEADD(3)=-psi;
FVertexEADD(2)=0;   
%A
FVertexADED(4)=-phi;  
FVertexADED(9)=-psi;
FVertexADED(11)=0;    
%Conexions with A+++++++++++++++++++++++++++++++++++++++++++++++++++
FDDlinea20(1)=-Br; %v;
FDDlinea20(2)=-Hz; %v;
FDDlinea20(3)=-F; %v;
FDDlinea20(4)=-G; %v;
FDDlinea20(5)=Mz;
%A
FADlinea20(4)=-phi;
FADlinea20(9)=-psi;
FADlinea20(11)=0;
%Vertex DAB  psiD=psiE dpsiDr=dpsiEeR   phiA=psiD+++++++++++++++++++++++
FVertexDAB(1)=-dphidr;  % conecting with B
FVertexDAB(2)=-dpsidr;
FVertexDAB(3)=-F;    
FVertexDAB(4)=-G;
FVertexDAB(5)=Mz;
%B
FVertexBDAD(1)=-phi;  % conecting with FDEEVertexDEA
FVertexBDAD(4)=-psi;

%A
FVertexABDD(4)=-phi;  % conecting with FDEEVertexDEA
FVertexABDD(9)=-psi;
FVertexABDD(11)=0;
%Conexions with  B++++++++++++++++++++++++++++++++++++++++++++++++++
FDDlinea10(1)=-dphidr; %v;
FDDlinea10(2)=-dpsidr; %v;
FDDlinea10(3)=-F; %v;
FDDlinea10(4)=-G; %v;
FDDlinea10(5)=Mz;
%B
FBDlinea10(1)=-phi;
FBDlinea10(4)=-psi;


%VertexDBC  psiD=psiC dpsiDr=dpsiCr phiB=psiD+++++++++++++++++++++++++
FVertexDBC(1)=-phi; %v; in C+++++++++++++++++++++++++++++++++
FVertexDBC(2)=-psi; %v;
FVertexDBC(3)=-F; %v;
FVertexDBC(4)=-G; %v;
FVertexDBC(5)=Mz;
%B
FVertexBCDD(1)=-phi;
FVertexBCDD(4)=-psi;
%C
FVertexCDBD(1)=-dphidr;
FVertexCDBD(3)=-dpsidr;

%Conexions with  C++++++++++++++++++++++++++++++++++++++++++++++++++
FDDlinea13(1)=-dphidr; %v;
FDDlinea13(2)=-dpsidr; %v;
FDDlinea13(3)=-F; %v;
FDDlinea13(4)=-G; %v;
FDDlinea13(5)=Mz;
%C
FCDlinea13(1)=-phi;
FCDlinea13(3)=-psi;

%superior
FDDlinea7_c(1)=psi-psi_c;
FDDlinea7_c(2)=phi-phid_c;
FDDlinea7_c(3)=F-Fcoil;dFdr0;
FDDlinea7_c(4)=G-H2;
FDDlinea7_c(5)=Mz;

FDDlinea7_m(1)=psi-psi_m;
FDDlinea7_m(2)=phi-phid_m;
FDDlinea7_m(3)=F-Fcoil;dFdr0;
FDDlinea7_m(4)=G-H2;
FDDlinea7_m(5)=Mz;

%lateral
FDDlinea8_c(1)=psi-psi_c;%dpsi0z+Hzm_m;
FDDlinea8_c(2)=phi-phid_c;
FDDlinea8_c(3)=F-Rout;
FDDlinea8_c(4)=dGdz0;
FDDlinea8_c(5)=Mz;

FDDlinea8_m(1)=psi-psi_m;%dpsi0z+Hzm_m;
FDDlinea8_m(2)=phi-phid_m;
FDDlinea8_m(3)=F-Rout;
FDDlinea8_m(4)=dGdz0;
FDDlinea8_m(5)=Mz;

FDDlinea8_center_c(1) = psi-psi_c;
FDDlinea8_center_c(2) = phi-phid_c;
FDDlinea8_center_c(3) = F-Rout;
FDDlinea8_center_c(4) = dGdz0;
FDDlinea8_center_c(5)=Mz;

FDDlinea8_center_m(1) = psi-psi_m;
FDDlinea8_center_m(2) = phi-phid_m;
FDDlinea8_center_m(3) = F-Rout;
FDDlinea8_center_m(4) = dGdz0;
FDDlinea8_center_m(5)=Mz;

%pared inferior
FDDlinea21_c(1)=psi-psi_c;
FDDlinea21_c(2)=phi-phid_c;
FDDlinea21_c(3)=F-Fcoil;%  dFdr0;
FDDlinea21_c(4)=G-zout;
FDDlinea21_c(5)=Mz;

FDDlinea21_m(1)=psi-psi_m;
FDDlinea21_m(2)=phi-phid_m;
FDDlinea21_m(3)=F-Fcoil;%  dFdr0;
FDDlinea21_m(4)=G-zout;
FDDlinea21_m(5)=Mz;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=reshape(yfD',NVAR *NDER,1); %variables

%x=[yf1(5,1);yf1(6,1)]


%jacobians
for k=1:NVAR
    
    dFDD(k,:)=jacobian(FDD(k),x);
    dFDDCoil(k,:)=jacobian(FDDCoil(k),x);
    dFDDCoil1(k,:)=jacobian(FDDCoil1(k),x);
    dFDDCoilContour(k,:)=jacobian(FDDCoilContour(k),x);
    dFDDRightCoil(k,:)=jacobian(FDDRightCoil(k),x);
    dFDDlinea21_c(k,:)=jacobian(FDDlinea21_c(k),x);
    dFDDlinea21_m(k,:)=jacobian(FDDlinea21_m(k),x);
     dFDDlinea7_c(k,:)=jacobian(FDDlinea7_c(k),x);
    dFDDlinea7_m(k,:)=jacobian(FDDlinea7_m(k),x);
    dFDDlinea8_c(k,:)=jacobian(FDDlinea8_c(k),x);
    dFDDlinea8_m(k,:)=jacobian(FDDlinea8_m(k),x);
    dFDDlinea8_center_c(k,:) = jacobian(FDDlinea8_center_c(k),x);
    dFDDlinea8_center_m(k,:) = jacobian(FDDlinea8_center_m(k),x);
   
    %block A
    dFDDlinea20(k,:)=jacobian(FDDlinea20(k),x);
    %block B
    dFDDlinea10(k,:)=jacobian(FDDlinea10(k),x);
     %block C
    dFDDlinea13(k,:)=jacobian(FDDlinea13(k),x);
    %block E
    dFDDlinea19(k,:)=jacobian(FDDlinea19(k),x);
    %vertices
    dFVertexDAB(k,:)=jacobian(FVertexDAB(k),x);
    dFVertexDBC(k,:)=jacobian(FVertexDBC(k),x);
    dFVertexDEA(k,:)=jacobian(FVertexDEA(k),x);
end


matlabFunction(dHrdz-dHzdr,'file',  [path_jacobian 'curlH_D.m'],'vars',{z0,r0,x,pa});
matlabFunction(Hr,Hz,'file',  [path_jacobian 'MagneticfieldDD.m'],'vars',{z0,r0,x,pa});
matlabFunction(Hra,Hza,'file',  [path_jacobian 'MagneticfieldDDa.m'],'vars',{z0,r0,x,pa});
matlabFunction(dphidr,dphidz,dphidrr,dphidzz,dphidrz,dpsidr,dpsidz,dpsidrr,dpsidzz,dpsidrz,...
    'file',[path_jacobian 'potentialsDD.m'],'vars',{z0,r0,x,pa});

matlabFunction(FDD,dFDD,'file',  [path_jacobian 'equationFDD.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDCoil,dFDDCoil,'file',  [path_jacobian 'equationFDDCoil.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDCoil1,dFDDCoil1,'file',  [path_jacobian 'equationFDDCoil1.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDCoilContour,dFDDCoilContour,'file',  [path_jacobian 'equationFDDCoilContour.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDRightCoil,dFDDRightCoil,'file',  [path_jacobian 'equationFDDRightCoil.m'],'vars',{z0,r0,x,pa});

matlabFunction(FDDlinea21_c,dFDDlinea21_c,'file',  [path_jacobian 'equationFDDlinea21_c.m'],'vars',{z0,r0,x,pa,phid_c});
matlabFunction(FDDlinea21_m,dFDDlinea21_m,'file',  [path_jacobian 'equationFDDlinea21_m.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDlinea19,dFDDlinea19,'file',  [path_jacobian 'equationFDDlinea19.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDlinea7_c,dFDDlinea7_c,'file',  [path_jacobian 'equationFDDlinea7_c.m'],'vars',{z0,r0,x,pa,phid_c});
matlabFunction(FDDlinea7_m,dFDDlinea7_m,'file',  [path_jacobian 'equationFDDlinea7_m.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDlinea8_c,dFDDlinea8_c,'file',  [path_jacobian 'equationFDDlinea8_c.m'],'vars',{z0,r0,x,pa,phid_c});
matlabFunction(FDDlinea8_m,dFDDlinea8_m,'file',  [path_jacobian 'equationFDDlinea8_m.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDlinea8_center_c,dFDDlinea8_center_c,'file',  [path_jacobian 'equationFDDlinea8_center_c.m'],'vars',{z0,r0,x,pa,phid_c});
matlabFunction(FDDlinea8_center_m,dFDDlinea8_center_m,'file',  [path_jacobian 'equationFDDlinea8_center_m.m'],'vars',{z0,r0,x,pa});


%A
matlabFunction(FDDlinea20,dFDDlinea20,'file',  [path_jacobian 'equationFDDlinea20.m'],'vars',{z0,r0,x,pa});
%B
matlabFunction(FDDlinea10,dFDDlinea10,'file',  [path_jacobian 'equationFDDlinea10.m'],'vars',{z0,r0,x,pa});
%C
matlabFunction(FDDlinea13,dFDDlinea13,'file',  [path_jacobian 'equationFDDlinea13.m'],'vars',{z0,r0,x,pa});
%E
matlabFunction(FDDlinea19,dFDDlinea19,'file',  [path_jacobian 'equationFDDlinea19.m'],'vars',{z0,r0,x,pa});
%vertices

matlabFunction(FVertexDAB,dFVertexDAB,'file',[path_jacobian 'equationFVertexDAB.m'],'vars',{z0,r0,x,pa});
matlabFunction(FVertexDBC,dFVertexDBC,'file',[path_jacobian 'equationFVertexDBC.m'],'vars',{z0,r0,x,pa});
matlabFunction(FVertexDEA,dFVertexDEA,'file',[path_jacobian 'equationFVertexDEA.m'],'vars',{z0,r0,x,pa});


%conexions with B
 for k=1:NVB
 dFBDlinea10(k,:)=jacobian(FBDlinea10(k),x);
 dFVertexBDAD(k,:)=jacobian(FVertexBDAD(k),x);
 dFVertexBCDD(k,:)=jacobian(FVertexBCDD(k),x);
 end
 matlabFunction(FBDlinea10,dFBDlinea10,'file',[path_jacobian 'equationFBDlinea10.m'],'vars',{z0,r0,x,pa});
 matlabFunction(FVertexBDAD,dFVertexBDAD,'file',[path_jacobian 'equationFVertexBDAD.m'],'vars',{z0,r0,x,pa});
 matlabFunction(FVertexBCDD,dFVertexBCDD,'file',[path_jacobian 'equationFVertexBCDD.m'],'vars',{z0,r0,x,pa});
%conexions with B
for k=1:NVC
 dFCDlinea13(k,:)=jacobian(FCDlinea13(k),x);
 dFVertexCDBD(k,:)=jacobian(FVertexCDBD(k),x);
 end
 matlabFunction(FCDlinea13,dFCDlinea13,'file',[path_jacobian 'equationFCDlinea13.m'],'vars',{z0,r0,x,pa});
 matlabFunction(FVertexCDBD,dFVertexCDBD,'file',[path_jacobian 'equationFVertexCDBD.m'],'vars',{z0,r0,x,pa});

%conexions with A
for k=1:NVA
dFADlinea20(k,:)=jacobian(FADlinea20(k),x);
dFVertexADED(k,:)=jacobian(FVertexADED(k),x);
dFVertexABDD(k,:)=jacobian(FVertexABDD(k),x);

end
matlabFunction(FADlinea20,dFADlinea20,'file',[path_jacobian 'equationFADlinea20.m'],'vars',{z0,r0,x,pa});
matlabFunction(FVertexADED,dFVertexADED,'file',[path_jacobian 'equationFVertexADED.m'],'vars',{z0,r0,x,pa});
matlabFunction(FVertexABDD,dFVertexABDD,'file',[path_jacobian 'equationFVertexABDD.m'],'vars',{z0,r0,x,pa});


%conexions with E
for k=1:NVE
dFEDlinea19(k,:)=jacobian(FEDlinea19(k),x);
dFVertexEADD(k,:)=jacobian(FVertexEADD(k),x);
end
matlabFunction(FEDlinea19,dFEDlinea19,'file',[path_jacobian 'equationFEDlinea19.m'],'vars',{z0,r0,x,pa});
matlabFunction(FVertexEADD,dFVertexEADD,'file',[path_jacobian 'equationFVertexEADD.m'],'vars',{z0,r0,x,pa});



