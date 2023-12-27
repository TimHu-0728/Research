% Path locations
path_jacobian  = './eigen/jacobians2Dnew/';
path_jacobian1 = './subroutinesMatlab/';
path_jacobian1 = 'C:\Users\miguel\Dropbox/subroutinesMatlab/';
path_meshD     = './meshless_D';

% Path operations
restoredefaultpath
addpath(path_meshD)
addpath(path_jacobian)
addpath(path_jacobian1)
name='soluciones/21-Jul-2020_Ib2000zeta0tau0nrA81nzA81nrB91nrC81nrE51nzE51'

load(name)

n=1
switch n
    case 1
        hold on
plot (tiempo,hmin)
case 2
    a=3
visualization(rA, rB, rC, rD, rE, zA, zB, zC, zD, zE, nrA, nrB,...
    nrC, nrE, nzA, nzB, nzC, nzE, phiA, phiB, phiC, phiD, phiE, HrA, HrB,...
    HrC, HrD, HrE, HzA, HzB, HzC, HzD, HzE, MAz0, MAr0, gA, fA, Forcer,...
    Forcez, pA, etaA, P1, P2, la, Ib_fe)
    case 3
      rA=reshape(rA,nrA,nzA);  
      zA=reshape(zA,nrA,nzA);
      contourf(rA,zA,uA)
end