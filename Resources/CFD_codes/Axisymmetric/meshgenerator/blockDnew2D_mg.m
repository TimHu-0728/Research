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
FDD(1)=F-Fcoil;g22*dFdzz0+g11*dFdrr0-2*g12*dFdrz0-Q1;
FDD(2)=g22*dGdzz0+g11*dGdrr0-2*g12*dGdrz0;

FDDa(1)=dFdzz0+dFdrr0;
FDDa(2)=dGdzz0+dGdrr0;

FDDCoil(1) = F-Fcoil;
FDDCoil(2) = G-Gcoil;
FDDCoil1   = FDDCoil;

FDDRightCoil(1) = F-Fcoil;
FDDRightCoil(2) = G-Gcoil;

%conexion wih E+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
FDDlinea19(1)=F-Fcoil;     % conecting with EDD -F
FDDlinea19(2)=G-Gcoil;     % conecting with EDD -G

% conecting with EDD -F
%Vertex DEB  psiD=psiE dpsiDr=dpsiEeR%phiA=psiD+++++++++++++++++++++++++
FVertexDEA(1)=F-Fcoil;
FVertexDEA(2)=G-Gcoil;

%Conexions with A+++++++++++++++++++++++++++++++++++++++++++++++++++
FDDlinea20(1)=F-Fcoil; %v;
FDDlinea20(2)=G-Gcoil; %v;

%Vertex DAB  psiD=psiE dpsiDr=dpsiEeR   phiA=psiD+++++++++++++++++++++++
FVertexDAB(1)=F-Fcoil;
FVertexDAB(2)=G-Gcoil;

%Conexions with  B++++++++++++++++++++++++++++++++++++++++++++++++++
FDDlinea10(1)=F-Fcoil; %v;
FDDlinea10(2)=G-Gcoil; %v;

%VertexDBC  psiD=psiC dpsiDr=dpsiCr phiB=psiD+++++++++++++++++++++++++
FVertexDBC(1)=F-Fcoil; %v;
FVertexDBC(2)=G-Gcoil;%v;


%Conexions with  C++++++++++++++++++++++++++++++++++++++++++++++++++
FDDlinea13(1)=F-Fcoil; %v;
FDDlinea13(2)=G-Gcoil; %v;

%superior
FDDlinea7(1)=F-Fcoil;dFdr0;
FDDlinea7(2)=G-H2;

%lateral
FDDlinea8(1)=F-Rout;
FDDlinea8(2)=G-Gcoil;

%pared inferior
FDDlinea21(1)=F-Fcoil;dFdr0;
FDDlinea21(2)=G-zout;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=reshape(yfD',NVAR *NDER,1); %variables

%x=[yf1(5,1);yf1(6,1)]

%jacobians
for k=1:NVAR
    
    dFDD(k,:)=jacobian(FDD(k),x);
    dFDDa(k,:)=jacobian(FDDa(k),x);
    dFDDCoil(k,:)=jacobian(FDDCoil(k),x);
    dFDDCoil1(k,:)=jacobian(FDDCoil1(k),x);
    dFDDRightCoil(k,:)=jacobian(FDDRightCoil(k),x);
    dFDDlinea21(k,:)=jacobian(FDDlinea21(k),x);
    
    dFDDlinea7(k,:)=jacobian(FDDlinea7(k),x);
    dFDDlinea8(k,:)=jacobian(FDDlinea8(k),x);
    
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







matlabFunction(FDD,dFDD,'file',  [path_jacobian 'equationFDD_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDa,dFDDa,'file',  [path_jacobian 'equationFDDa_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDCoil,dFDDCoil,'file',  [path_jacobian 'equationFDDCoil_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDCoil1,dFDDCoil1,'file',  [path_jacobian 'equationFDDCoil1_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDRightCoil,dFDDRightCoil,'file',  [path_jacobian 'equationFDDRightCoil_mg.m'],'vars',{z0,r0,x,pa});

matlabFunction(FDDlinea20,dFDDlinea20,'file',  [path_jacobian 'equationFDDlinea20_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDlinea10,dFDDlinea10,'file',  [path_jacobian 'equationFDDlinea10_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDlinea13,dFDDlinea13,'file',  [path_jacobian 'equationFDDlinea13_mg.m'],'vars',{z0,r0,x,pa});

matlabFunction(FDDlinea21,dFDDlinea21,'file',  [path_jacobian 'equationFDDlinea21_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDlinea19,dFDDlinea19,'file',  [path_jacobian 'equationFDDlinea19_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDlinea7,dFDDlinea7,'file',  [path_jacobian 'equationFDDlinea7_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FDDlinea8,dFDDlinea8,'file',  [path_jacobian 'equationFDDlinea8_mg.m'],'vars',{z0,r0,x,pa});


matlabFunction(FVertexDAB,dFVertexDAB,'file',[path_jacobian 'equationFVertexDAB_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FVertexDBC,dFVertexDBC,'file',[path_jacobian 'equationFVertexDBC_mg.m'],'vars',{z0,r0,x,pa});
matlabFunction(FVertexDEA,dFVertexDEA,'file',[path_jacobian 'equationFVertexDEA_mg.m'],'vars',{z0,r0,x,pa});