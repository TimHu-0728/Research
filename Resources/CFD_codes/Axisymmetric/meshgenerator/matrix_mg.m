%% Time derivative 
 dt2=dt1+dt; 
 bm= -(dt2/dt)/(dt2-dt); 
 
 
  bmm= (dt/dt2)/(dt2-dt); 
 bp=-((dt2/dt)^2-1)/((dt2/dt)*(dt-dt2)); 
 
 
 xt = bp*x0 + bm*x0m + bmm*x0mm;
 
 % Initializing the variables 
 
nbi = 0;
yfD= zeros(NVD,NDD,ntD);

inic = nbi+0*ntD+1; 
ifin = nbi+1*ntD;
FD= x0(inic:ifin);
yfD(1,1,:)=FD;
yfD(1,2,: ) = ddr0D*FD;
yfD(1,3,: ) = ddz0D*FD;
yfD(1,4,: ) = ddrr0D*FD;
yfD(1,5,: ) = ddzz0D*FD;
yfD(1,6,: ) = ddrz0D*FD;
yfD(1,7,: ) = xt(inic:ifin);
 
 
inic = nbi+1*ntD+1; 
ifin = nbi+2*ntD;
GD= x0(inic:ifin);
yfD(2,1,:)=GD;
yfD(2,2,: ) = ddr0D*GD;
yfD(2,3,: ) = ddz0D*GD;
yfD(2,4,: ) = ddrr0D*GD;
yfD(2,5,: ) = ddzz0D*GD;
yfD(2,6,: ) = ddrz0D*GD;
yfD(2,7,: ) = xt(inic:ifin);
 
 
nbi = nbi + ntD*NVD; 
 

%
% Initializing the jacobian submatrix 
%
 
dd = 0*speye(ntD*NVD,ntD*NVD); 

% Block D
 
FDD= zeros(NVD,ntD);
DFDD= zeros(NVD,NDD*NVD,ntD);

 % bulk 
for i=2:nzD-1 
for j=2:nrD-1 
l=sub2ind([nrD,nzD],j,i);
x = reshape(yfD(:,:,l)',NDD*NVD,1);
[FDD(:,l),DFDD(:,:,l)]=equationFDD(z0D(i), r0D(j),x,pa);
end

%
% *top (j=NR)* and *bottom(j=1)* 
%
 
l=sub2ind([nrD,nzD],1,i); % bottom
x = reshape(yfD(:,:,l)',NDD*NVD,1);
[FDD(:,l),DFDD(:,:,l)]=equationFDDb(z0D(i), r0D(1),x,pa);

l=sub2ind([nrD,nzD],nrD,i); %top 
 x = reshape(yfD(:,:,l)',NDD*NVD,1);
[FDD(:,l),DFDD(:,:,l)]=equationFDDt(z0D(i), r0D(nrD),x,pa);
end

%
% *Left (i=1)* and *right(i=NZ)* 
%
 
for j=2:nrD-1 
l = sub2ind([nrD,nzD],j,1); 
x = reshape(yfD(:,:,l)',NDD*NVD,1);
[FDD(:,l),DFDD(:,:,l)]=equationFDDl(z0D(1), r0D(j),x,pa);
 
l = sub2ind([nrD,nzD],j,nzD); 
x = reshape(yfD(:,:,l)',NDD*NVD,1);
[FDD(:,l),DFDD(:,:,l)]=equationFDDr(z0D(nzD), r0D(j),x,pa);
end

%corners (Must match with those in connections!) 


 % corner (z) i=1 and (r) j=1 

l = sub2ind([nrD,nzD],1,1); 
x = reshape(yfD(:,:,l)',NDD*NVD,1);
[FDD(:,l),DFDD(:,:,l)]=equationFDDl(z0D(1), r0D(1),x,pa);

 % corner (z) i=nrD and (r) j=1 

l = sub2ind([nrD,nzD],nrD,1); 
x = reshape(yfD(:,:,l)',NDD*NVD,1);
[FDD(:,l),DFDD(:,:,l)]=equationFDDl(z0D(1), r0D(nrD),x,pa);

 % corner (z) i=1 and (r) j=nzD 

l = sub2ind([nrD,nzD],1,nzD); 
x = reshape(yfD(:,:,l)',NDD*NVD,1);
[FDD(:,l),DFDD(:,:,l)]=equationFDDr(z0D(nzD), r0D(1),x,pa);

 % corner (z) i=nrD and (r) j=nzD 

l = sub2ind([nrD,nzD],nrD,nzD); 
x = reshape(yfD(:,:,l)',NDD*NVD,1);
[FDD(:,l),DFDD(:,:,l)]=equationFDDr(z0D(nzD), r0D(nrD),x,pa);

%
% Mounting the jacobians 
%
 
xd = 0*speye(ntD*NVD,1); 

for jj=1:NVD 
rowi = (jj-1)*ntD + 1;
rowf = jj*ntD;
for kk=1:NVD 
km=(kk-1)*NDD+1;
kp=kk*NDD;
D=squeeze(DFDD(jj,km:kp,:));;
B=spdiags(D(1,:)',0,ntD, ntD) ...
+ spdiags(D(2,:)',0, ntD, ntD)* ddr0D...
+ spdiags(D(3,:)',0, ntD, ntD)* ddz0D...
+ spdiags(D(4,:)',0, ntD, ntD)* ddrr0D...
+ spdiags(D(5,:)',0, ntD, ntD)* ddzz0D...
+ spdiags(D(6,:)',0, ntD, ntD)* ddrz0D...
+ spdiags(D(7,:)',0, ntD, ntD)* bp;
coli = (kk-1)*ntD + 1;
colf = kk*ntD;
dd(rowi:rowf,coli:colf) = sparse(B);
end
xd(rowi:rowf) = - FDD(jj,:)';
end

%
% Connections 
%
 
% 
% Mounting of vector b and off diagonal jacobian
% 


a = [ dd ;  ];
b = [ xd; ];
