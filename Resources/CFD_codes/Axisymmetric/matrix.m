%% Time derivative 
 dt2=dt1+dt; 
 bm= -(dt2/dt)/(dt2-dt); 
 
 
  bmm= (dt/dt2)/(dt2-dt); 
 bp=-((dt2/dt)^2-1)/((dt2/dt)*(dt-dt2)); 
 
 
 xt = bp*x0 + bm*x0m + bmm*x0mm;
 
 % Initializing the variables 
 
nbi = 0;
yfA= zeros(NVA,NDA,ntA);

inic = nbi+0*ntA+1; 
ifin = nbi+1*ntA;
wA= x0(inic:ifin);
yfA(1,1,:)=wA;
yfA(1,2,: ) = dd0rA*wA;
yfA(1,3,: ) = dd0zA*wA;
yfA(1,4,: ) = dd0rrA*wA;
yfA(1,5,: ) = dd0zzA*wA;
yfA(1,6,: ) = dd0rzA*wA;
yfA(1,7,: ) = xt(inic:ifin);
 
 
inic = nbi+1*ntA+1; 
ifin = nbi+2*ntA;
uA= x0(inic:ifin);
yfA(2,1,:)=uA;
yfA(2,2,: ) = dd0rA*uA;
yfA(2,3,: ) = dd0zA*uA;
yfA(2,4,: ) = dd0rrA*uA;
yfA(2,5,: ) = dd0zzA*uA;
yfA(2,6,: ) = dd0rzA*uA;
yfA(2,7,: ) = xt(inic:ifin);
 
 
inic = nbi+2*ntA+1; 
ifin = nbi+3*ntA;
pA= x0(inic:ifin);
yfA(3,1,:)=pA;
yfA(3,2,: ) = dd0rA*pA;
yfA(3,3,: ) = dd0zA*pA;
yfA(3,4,: ) = dd0rrA*pA;
yfA(3,5,: ) = dd0zzA*pA;
yfA(3,6,: ) = dd0rzA*pA;
yfA(3,7,: ) = xt(inic:ifin);
 
 
inic = nbi+3*ntA+1; 
ifin = nbi+4*ntA;
phiA= x0(inic:ifin);
yfA(4,1,:)=phiA;
yfA(4,2,: ) = dd0rA*phiA;
yfA(4,3,: ) = dd0zA*phiA;
yfA(4,4,: ) = dd0rrA*phiA;
yfA(4,5,: ) = dd0zzA*phiA;
yfA(4,6,: ) = dd0rzA*phiA;
yfA(4,7,: ) = xt(inic:ifin);
 
 
inic = nbi+4*ntA+1; 
ifin = nbi+5*ntA;
fA= x0(inic:ifin);
yfA(5,1,:)=fA;
yfA(5,2,: ) = dd0rA*fA;
yfA(5,3,: ) = dd0zA*fA;
yfA(5,4,: ) = dd0rrA*fA;
yfA(5,5,: ) = dd0zzA*fA;
yfA(5,6,: ) = dd0rzA*fA;
yfA(5,7,: ) = xt(inic:ifin);
 
 
inic = nbi+5*ntA+1; 
ifin = nbi+6*ntA;
gA= x0(inic:ifin);
yfA(6,1,:)=gA;
yfA(6,2,: ) = dd0rA*gA;
yfA(6,3,: ) = dd0zA*gA;
yfA(6,4,: ) = dd0rrA*gA;
yfA(6,5,: ) = dd0zzA*gA;
yfA(6,6,: ) = dd0rzA*gA;
yfA(6,7,: ) = xt(inic:ifin);
 
 
inic = nbi+6*ntA+1; 
ifin = nbi+7*ntA;
CA= x0(inic:ifin);
yfA(7,1,:)=CA;
yfA(7,2,: ) = dd0rA*CA;
yfA(7,3,: ) = dd0zA*CA;
yfA(7,4,: ) = dd0rrA*CA;
yfA(7,5,: ) = dd0zzA*CA;
yfA(7,6,: ) = dd0rzA*CA;
yfA(7,7,: ) = xt(inic:ifin);
 
 
inic = nbi+7*ntA+1; 
ifin = nbi+8*ntA;
etaA= x0(inic:ifin);
yfA(8,1,:)=etaA;
yfA(8,2,: ) = dd0rA*etaA;
yfA(8,3,: ) = dd0zA*etaA;
yfA(8,4,: ) = dd0rrA*etaA;
yfA(8,5,: ) = dd0zzA*etaA;
yfA(8,6,: ) = dd0rzA*etaA;
yfA(8,7,: ) = xt(inic:ifin);
 
 
inic = nbi+8*ntA+1; 
ifin = nbi+9*ntA;
psiA= x0(inic:ifin);
yfA(9,1,:)=psiA;
yfA(9,2,: ) = dd0rA*psiA;
yfA(9,3,: ) = dd0zA*psiA;
yfA(9,4,: ) = dd0rrA*psiA;
yfA(9,5,: ) = dd0zzA*psiA;
yfA(9,6,: ) = dd0rzA*psiA;
yfA(9,7,: ) = xt(inic:ifin);
 
 
inic = nbi+9*ntA+1; 
ifin = nbi+10*ntA;
VA= x0(inic:ifin);
yfA(10,1,:)=VA;
yfA(10,2,: ) = dd0rA*VA;
yfA(10,3,: ) = dd0zA*VA;
yfA(10,4,: ) = dd0rrA*VA;
yfA(10,5,: ) = dd0zzA*VA;
yfA(10,6,: ) = dd0rzA*VA;
yfA(10,7,: ) = xt(inic:ifin);
 
 
inic = nbi+10*ntA+1; 
ifin = nbi+11*ntA;
gaxisA= x0(inic:ifin);
yfA(11,1,:)=gaxisA;
yfA(11,2,: ) = dd0rA*gaxisA;
yfA(11,3,: ) = dd0zA*gaxisA;
yfA(11,4,: ) = dd0rrA*gaxisA;
yfA(11,5,: ) = dd0zzA*gaxisA;
yfA(11,6,: ) = dd0rzA*gaxisA;
yfA(11,7,: ) = xt(inic:ifin);
 
 
nbi = nbi + ntA*NVA; 
 
yfB= zeros(NVB,NDB,ntB);

inic = nbi+0*ntB+1; 
ifin = nbi+1*ntB;
phiB= x0(inic:ifin);
yfB(1,1,:)=phiB;
yfB(1,2,: ) = dd0rB*phiB;
yfB(1,3,: ) = dd0zB*phiB;
yfB(1,4,: ) = dd0rrB*phiB;
yfB(1,5,: ) = dd0zzB*phiB;
yfB(1,6,: ) = dd0rzB*phiB;
yfB(1,7,: ) = xt(inic:ifin);
 
 
inic = nbi+1*ntB+1; 
ifin = nbi+2*ntB;
fB= x0(inic:ifin);
yfB(2,1,:)=fB;
yfB(2,2,: ) = dd0rB*fB;
yfB(2,3,: ) = dd0zB*fB;
yfB(2,4,: ) = dd0rrB*fB;
yfB(2,5,: ) = dd0zzB*fB;
yfB(2,6,: ) = dd0rzB*fB;
yfB(2,7,: ) = xt(inic:ifin);
 
 
inic = nbi+2*ntB+1; 
ifin = nbi+3*ntB;
gB= x0(inic:ifin);
yfB(3,1,:)=gB;
yfB(3,2,: ) = dd0rB*gB;
yfB(3,3,: ) = dd0zB*gB;
yfB(3,4,: ) = dd0rrB*gB;
yfB(3,5,: ) = dd0zzB*gB;
yfB(3,6,: ) = dd0rzB*gB;
yfB(3,7,: ) = xt(inic:ifin);
 
 
inic = nbi+3*ntB+1; 
ifin = nbi+4*ntB;
psiB= x0(inic:ifin);
yfB(4,1,:)=psiB;
yfB(4,2,: ) = dd0rB*psiB;
yfB(4,3,: ) = dd0zB*psiB;
yfB(4,4,: ) = dd0rrB*psiB;
yfB(4,5,: ) = dd0zzB*psiB;
yfB(4,6,: ) = dd0rzB*psiB;
yfB(4,7,: ) = xt(inic:ifin);
 
 
nbi = nbi + ntB*NVB; 
 
yfC= zeros(NVC,NDC,ntC);

inic = nbi+0*ntC+1; 
ifin = nbi+1*ntC;
phiC= x0(inic:ifin);
yfC(1,1,:)=phiC;
yfC(1,2,: ) = dd0rC*phiC;
yfC(1,3,: ) = dd0zC*phiC;
yfC(1,4,: ) = dd0rrC*phiC;
yfC(1,5,: ) = dd0zzC*phiC;
yfC(1,6,: ) = dd0rzC*phiC;
yfC(1,7,: ) = xt(inic:ifin);
 
 
inic = nbi+1*ntC+1; 
ifin = nbi+2*ntC;
fC= x0(inic:ifin);
yfC(2,1,:)=fC;
yfC(2,2,: ) = dd0rC*fC;
yfC(2,3,: ) = dd0zC*fC;
yfC(2,4,: ) = dd0rrC*fC;
yfC(2,5,: ) = dd0zzC*fC;
yfC(2,6,: ) = dd0rzC*fC;
yfC(2,7,: ) = xt(inic:ifin);
 
 
inic = nbi+2*ntC+1; 
ifin = nbi+3*ntC;
psiC= x0(inic:ifin);
yfC(3,1,:)=psiC;
yfC(3,2,: ) = dd0rC*psiC;
yfC(3,3,: ) = dd0zC*psiC;
yfC(3,4,: ) = dd0rrC*psiC;
yfC(3,5,: ) = dd0zzC*psiC;
yfC(3,6,: ) = dd0rzC*psiC;
yfC(3,7,: ) = xt(inic:ifin);
 
 
nbi = nbi + ntC*NVC; 
 
yfD= zeros(NVD,NDD,ntD);

inic = nbi+0*ntD+1; 
ifin = nbi+1*ntD;
phiD= x0(inic:ifin);
yfD(1,1,:)=phiD;
yfD(1,2,: ) = dd0rD*phiD;
yfD(1,3,: ) = dd0zD*phiD;
yfD(1,4,: ) = dd0rrD*phiD;
yfD(1,5,: ) = dd0zzD*phiD;
yfD(1,6,: ) = dd0rzD*phiD;
yfD(1,7,: ) = xt(inic:ifin);
 
 
inic = nbi+1*ntD+1; 
ifin = nbi+2*ntD;
psiD= x0(inic:ifin);
yfD(2,1,:)=psiD;
yfD(2,2,: ) = dd0rD*psiD;
yfD(2,3,: ) = dd0zD*psiD;
yfD(2,4,: ) = dd0rrD*psiD;
yfD(2,5,: ) = dd0zzD*psiD;
yfD(2,6,: ) = dd0rzD*psiD;
yfD(2,7,: ) = xt(inic:ifin);
 
 
inic = nbi+2*ntD+1; 
ifin = nbi+3*ntD;
FD= x0(inic:ifin);
yfD(3,1,:)=FD;
yfD(3,2,: ) = dd0rD*FD;
yfD(3,3,: ) = dd0zD*FD;
yfD(3,4,: ) = dd0rrD*FD;
yfD(3,5,: ) = dd0zzD*FD;
yfD(3,6,: ) = dd0rzD*FD;
yfD(3,7,: ) = xt(inic:ifin);
 
 
inic = nbi+3*ntD+1; 
ifin = nbi+4*ntD;
GD= x0(inic:ifin);
yfD(4,1,:)=GD;
yfD(4,2,: ) = dd0rD*GD;
yfD(4,3,: ) = dd0zD*GD;
yfD(4,4,: ) = dd0rrD*GD;
yfD(4,5,: ) = dd0zzD*GD;
yfD(4,6,: ) = dd0rzD*GD;
yfD(4,7,: ) = xt(inic:ifin);
 
 
nbi = nbi + ntD*NVD; 
 
yfE= zeros(NVE,NDE,ntE);

inic = nbi+0*ntE+1; 
ifin = nbi+1*ntE;
phiE= x0(inic:ifin);
yfE(1,1,:)=phiE;
yfE(1,2,: ) = dd0rE*phiE;
yfE(1,3,: ) = dd0zE*phiE;
yfE(1,4,: ) = dd0rrE*phiE;
yfE(1,5,: ) = dd0zzE*phiE;
yfE(1,6,: ) = dd0rzE*phiE;
yfE(1,7,: ) = xt(inic:ifin);
 
 
inic = nbi+1*ntE+1; 
ifin = nbi+2*ntE;
fE= x0(inic:ifin);
yfE(2,1,:)=fE;
yfE(2,2,: ) = dd0rE*fE;
yfE(2,3,: ) = dd0zE*fE;
yfE(2,4,: ) = dd0rrE*fE;
yfE(2,5,: ) = dd0zzE*fE;
yfE(2,6,: ) = dd0rzE*fE;
yfE(2,7,: ) = xt(inic:ifin);
 
 
inic = nbi+2*ntE+1; 
ifin = nbi+3*ntE;
psiE= x0(inic:ifin);
yfE(3,1,:)=psiE;
yfE(3,2,: ) = dd0rE*psiE;
yfE(3,3,: ) = dd0zE*psiE;
yfE(3,4,: ) = dd0rrE*psiE;
yfE(3,5,: ) = dd0zzE*psiE;
yfE(3,6,: ) = dd0rzE*psiE;
yfE(3,7,: ) = xt(inic:ifin);
 
 
nbi = nbi + ntE*NVE; 
 

%
% Initializing the jacobian submatrix 
%
 
aa = 0*speye(ntA*NVA,ntA*NVA); 
ab = 0*speye(ntA*NVA,ntB*NVB); 
ac = 0*speye(ntA*NVA,ntC*NVC); 
ad = 0*speye(ntA*NVA,ntD*NVD); 
ae = 0*speye(ntA*NVA,ntE*NVE); 
ba = 0*speye(ntB*NVB,ntA*NVA); 
bb = 0*speye(ntB*NVB,ntB*NVB); 
bc = 0*speye(ntB*NVB,ntC*NVC); 
bd = 0*speye(ntB*NVB,ntD*NVD); 
be = 0*speye(ntB*NVB,ntE*NVE); 
ca = 0*speye(ntC*NVC,ntA*NVA); 
cb = 0*speye(ntC*NVC,ntB*NVB); 
cc = 0*speye(ntC*NVC,ntC*NVC); 
cd = 0*speye(ntC*NVC,ntD*NVD); 
ce = 0*speye(ntC*NVC,ntE*NVE); 
da = 0*speye(ntD*NVD,ntA*NVA); 
db = 0*speye(ntD*NVD,ntB*NVB); 
dc = 0*speye(ntD*NVD,ntC*NVC); 
dd = 0*speye(ntD*NVD,ntD*NVD); 
de = 0*speye(ntD*NVD,ntE*NVE); 
ea = 0*speye(ntE*NVE,ntA*NVA); 
eb = 0*speye(ntE*NVE,ntB*NVB); 
ec = 0*speye(ntE*NVE,ntC*NVC); 
ed = 0*speye(ntE*NVE,ntD*NVD); 
ee = 0*speye(ntE*NVE,ntE*NVE); 

% Block A
 
FAA= zeros(NVA,ntA);
DFAA= zeros(NVA,NDA*NVA,ntA);

 % bulk 
for i=2:nzA-1 
for j=2:nrA-1 
l=sub2ind([nrA,nzA],j,i);
x = reshape(yfA(:,:,l)',NDA*NVA,1);
[FAA(:,l),DFAA(:,:,l)]=equationFAA(z0A(i), r0A(j),x,pa);
end

%
% *top (j=NR)* and *bottom(j=1)* 
%
 
l=sub2ind([nrA,nzA],1,i); % bottom
x = reshape(yfA(:,:,l)',NDA*NVA,1);
[FAA(:,l),DFAA(:,:,l)]=equationFAAb(z0A(i), r0A(1),x,pa);

l=sub2ind([nrA,nzA],nrA,i); %top 
 x = reshape(yfA(:,:,l)',NDA*NVA,1);
[FAA(:,l),DFAA(:,:,l)]=equationFAAt(z0A(i), r0A(nrA),x,pa);
end

%
% *Left (i=1)* and *right(i=NZ)* 
%
 
for j=2:nrA-1 
l = sub2ind([nrA,nzA],j,1); 
x = reshape(yfA(:,:,l)',NDA*NVA,1);
[FAA(:,l),DFAA(:,:,l)]=equationFAAl(z0A(1), r0A(j),x,pa);
 
l = sub2ind([nrA,nzA],j,nzA); 
x = reshape(yfA(:,:,l)',NDA*NVA,1);
[FAA(:,l),DFAA(:,:,l)]=equationFAAr(z0A(nzA), r0A(j),x,pa);
end

%corners (Must match with those in connections!) 


 % corner (z) i=1 and (r) j=1 

l = sub2ind([nrA,nzA],1,1); 
x = reshape(yfA(:,:,l)',NDA*NVA,1);
[FAA(:,l),DFAA(:,:,l)]=equationFAAl(z0A(1), r0A(1),x,pa);

 % corner (z) i=nrA and (r) j=1 

l = sub2ind([nrA,nzA],nrA,1); 
x = reshape(yfA(:,:,l)',NDA*NVA,1);
[FAA(:,l),DFAA(:,:,l)]=equationFAAl(z0A(1), r0A(nrA),x,pa);

 % corner (z) i=1 and (r) j=nzA 

l = sub2ind([nrA,nzA],1,nzA); 
x = reshape(yfA(:,:,l)',NDA*NVA,1);
[FAA(:,l),DFAA(:,:,l)]=equationFAAr(z0A(nzA), r0A(1),x,pa);

 % corner (z) i=nrA and (r) j=nzA 

l = sub2ind([nrA,nzA],nrA,nzA); 
x = reshape(yfA(:,:,l)',NDA*NVA,1);
[FAA(:,l),DFAA(:,:,l)]=equationFAAr(z0A(nzA), r0A(nrA),x,pa);

%
% Mounting the jacobians 
%
 
xa = 0*speye(ntA*NVA,1); 

for jj=1:NVA 
rowi = (jj-1)*ntA + 1;
rowf = jj*ntA;
for kk=1:NVA 
km=(kk-1)*NDA+1;
kp=kk*NDA;
D=squeeze(DFAA(jj,km:kp,:));;
B=spdiags(D(1,:)',0,ntA, ntA) ...
+ spdiags(D(2,:)',0, ntA, ntA)* dd0rA...
+ spdiags(D(3,:)',0, ntA, ntA)* dd0zA...
+ spdiags(D(4,:)',0, ntA, ntA)* dd0rrA...
+ spdiags(D(5,:)',0, ntA, ntA)* dd0zzA...
+ spdiags(D(6,:)',0, ntA, ntA)* dd0rzA...
+ spdiags(D(7,:)',0, ntA, ntA)* bp;
coli = (kk-1)*ntA + 1;
colf = kk*ntA;
aa(rowi:rowf,coli:colf) = sparse(B);
end
xa(rowi:rowf) = - FAA(jj,:)';
end

%
% Connections 
%
 
% Connection of A with A 
%
 
DFAA= zeros(NVA,NDA*NVA,ntA);

%
% Constructing aa 
%
 
for jj=1:NVA
list = jl +(jj-1)*ntA;
for kk=1:NVA
km=(kk-1)*NDA+1;
kp=kk*NDA;
D=squeeze(DFAA(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntA, ntA) ...
+ spdiags(D(2,:)',0, ntA, ntA)* dd0rA...
+ spdiags(D(3,:)',0, ntA, ntA)* dd0zA...
+ spdiags(D(4,:)',0, ntA, ntA)* dd0rrA...
+ spdiags(D(5,:)',0, ntA, ntA)* dd0zzA...
+ spdiags(D(6,:)',0, ntA, ntA)* dd0rzA...
+ spdiags(D(7,:)',0, ntA, ntA)* bp;
aux = sparse(B);
coli = (kk-1)*ntA + 1;
colf = kk*ntA;
aa(list,coli:colf) = aux(jl_c,:);
end
xa(list) = xa(list) - FAA(jj,il)';
end
 
% Connection of A with A 
%
 
DFAA= zeros(NVA,NDA*NVA,ntA);

%
% Constructing aa 
%
 
for jj=1:NVA
list = jl +(jj-1)*ntA;
for kk=1:NVA
km=(kk-1)*NDA+1;
kp=kk*NDA;
D=squeeze(DFAA(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntA, ntA) ...
+ spdiags(D(2,:)',0, ntA, ntA)* dd0rA...
+ spdiags(D(3,:)',0, ntA, ntA)* dd0zA...
+ spdiags(D(4,:)',0, ntA, ntA)* dd0rrA...
+ spdiags(D(5,:)',0, ntA, ntA)* dd0zzA...
+ spdiags(D(6,:)',0, ntA, ntA)* dd0rzA...
+ spdiags(D(7,:)',0, ntA, ntA)* bp;
aux = sparse(B);
coli = (kk-1)*ntA + 1;
colf = kk*ntA;
aa(list,coli:colf) = aux(jl_c,:);
end
xa(list) = xa(list) - FAA(jj,il)';
end
 

% Block B
 
FBB= zeros(NVB,ntB);
DFBB= zeros(NVB,NDB*NVB,ntB);

 % bulk 
for i=2:nzB-1 
for j=2:nrB-1 
l=sub2ind([nrB,nzB],j,i);
x = reshape(yfB(:,:,l)',NDB*NVB,1);
[FBB(:,l),DFBB(:,:,l)]=equationFBB(z0B(i), r0B(j),x,pa);
end

%
% *top (j=NR)* and *bottom(j=1)* 
%
 
l=sub2ind([nrB,nzB],1,i); % bottom
x = reshape(yfB(:,:,l)',NDB*NVB,1);
[FBB(:,l),DFBB(:,:,l)]=equationFBBb(z0B(i), r0B(1),x,pa);

l=sub2ind([nrB,nzB],nrB,i); %top 
 x = reshape(yfB(:,:,l)',NDB*NVB,1);
[FBB(:,l),DFBB(:,:,l)]=equationFBBt(z0B(i), r0B(nrB),x,pa);
end

%
% *Left (i=1)* and *right(i=NZ)* 
%
 
for j=2:nrB-1 
l = sub2ind([nrB,nzB],j,1); 
x = reshape(yfB(:,:,l)',NDB*NVB,1);
[FBB(:,l),DFBB(:,:,l)]=equationFBBl(z0B(1), r0B(j),x,pa);
 
l = sub2ind([nrB,nzB],j,nzB); 
x = reshape(yfB(:,:,l)',NDB*NVB,1);
[FBB(:,l),DFBB(:,:,l)]=equationFBBr(z0B(nzB), r0B(j),x,pa);
end

%corners (Must match with those in connections!) 


 % corner (z) i=1 and (r) j=1 

l = sub2ind([nrB,nzB],1,1); 
x = reshape(yfB(:,:,l)',NDB*NVB,1);
[FBB(:,l),DFBB(:,:,l)]=equationFBBl(z0B(1), r0B(1),x,pa);

 % corner (z) i=nrB and (r) j=1 

l = sub2ind([nrB,nzB],nrB,1); 
x = reshape(yfB(:,:,l)',NDB*NVB,1);
[FBB(:,l),DFBB(:,:,l)]=equationFBBl(z0B(1), r0B(nrB),x,pa);

 % corner (z) i=1 and (r) j=nzB 

l = sub2ind([nrB,nzB],1,nzB); 
x = reshape(yfB(:,:,l)',NDB*NVB,1);
[FBB(:,l),DFBB(:,:,l)]=equationFBBr(z0B(nzB), r0B(1),x,pa);

 % corner (z) i=nrB and (r) j=nzB 

l = sub2ind([nrB,nzB],nrB,nzB); 
x = reshape(yfB(:,:,l)',NDB*NVB,1);
[FBB(:,l),DFBB(:,:,l)]=equationFBBr(z0B(nzB), r0B(nrB),x,pa);

%
% Mounting the jacobians 
%
 
xb = 0*speye(ntB*NVB,1); 

for jj=1:NVB 
rowi = (jj-1)*ntB + 1;
rowf = jj*ntB;
for kk=1:NVB 
km=(kk-1)*NDB+1;
kp=kk*NDB;
D=squeeze(DFBB(jj,km:kp,:));;
B=spdiags(D(1,:)',0,ntB, ntB) ...
+ spdiags(D(2,:)',0, ntB, ntB)* dd0rB...
+ spdiags(D(3,:)',0, ntB, ntB)* dd0zB...
+ spdiags(D(4,:)',0, ntB, ntB)* dd0rrB...
+ spdiags(D(5,:)',0, ntB, ntB)* dd0zzB...
+ spdiags(D(6,:)',0, ntB, ntB)* dd0rzB...
+ spdiags(D(7,:)',0, ntB, ntB)* bp;
coli = (kk-1)*ntB + 1;
colf = kk*ntB;
bb(rowi:rowf,coli:colf) = sparse(B);
end
xb(rowi:rowf) = - FBB(jj,:)';
end

%
% Connections 
%
 
% Connection of B with B 
%
 
DFBB= zeros(NVB,NDB*NVB,ntB);

%
% Constructing bb 
%
 
for jj=1:NVB
list = jl +(jj-1)*ntB;
for kk=1:NVB
km=(kk-1)*NDB+1;
kp=kk*NDB;
D=squeeze(DFBB(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntB, ntB) ...
+ spdiags(D(2,:)',0, ntB, ntB)* dd0rB...
+ spdiags(D(3,:)',0, ntB, ntB)* dd0zB...
+ spdiags(D(4,:)',0, ntB, ntB)* dd0rrB...
+ spdiags(D(5,:)',0, ntB, ntB)* dd0zzB...
+ spdiags(D(6,:)',0, ntB, ntB)* dd0rzB...
+ spdiags(D(7,:)',0, ntB, ntB)* bp;
aux = sparse(B);
coli = (kk-1)*ntB + 1;
colf = kk*ntB;
bb(list,coli:colf) = aux(jl_c,:);
end
xb(list) = xb(list) - FBB(jj,il)';
end
 
% Connection of B with B 
%
 
DFBB= zeros(NVB,NDB*NVB,ntB);

%
% Constructing bb 
%
 
for jj=1:NVB
list = jl +(jj-1)*ntB;
for kk=1:NVB
km=(kk-1)*NDB+1;
kp=kk*NDB;
D=squeeze(DFBB(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntB, ntB) ...
+ spdiags(D(2,:)',0, ntB, ntB)* dd0rB...
+ spdiags(D(3,:)',0, ntB, ntB)* dd0zB...
+ spdiags(D(4,:)',0, ntB, ntB)* dd0rrB...
+ spdiags(D(5,:)',0, ntB, ntB)* dd0zzB...
+ spdiags(D(6,:)',0, ntB, ntB)* dd0rzB...
+ spdiags(D(7,:)',0, ntB, ntB)* bp;
aux = sparse(B);
coli = (kk-1)*ntB + 1;
colf = kk*ntB;
bb(list,coli:colf) = aux(jl_c,:);
end
xb(list) = xb(list) - FBB(jj,il)';
end
 

% Block C
 
FCC= zeros(NVC,ntC);
DFCC= zeros(NVC,NDC*NVC,ntC);

 % bulk 
for i=2:nzC-1 
for j=2:nrC-1 
l=sub2ind([nrC,nzC],j,i);
x = reshape(yfC(:,:,l)',NDC*NVC,1);
[FCC(:,l),DFCC(:,:,l)]=equationFCC(z0C(i), r0C(j),x,pa);
end

%
% *top (j=NR)* and *bottom(j=1)* 
%
 
l=sub2ind([nrC,nzC],1,i); % bottom
x = reshape(yfC(:,:,l)',NDC*NVC,1);
[FCC(:,l),DFCC(:,:,l)]=equationFCCb(z0C(i), r0C(1),x,pa);

l=sub2ind([nrC,nzC],nrC,i); %top 
 x = reshape(yfC(:,:,l)',NDC*NVC,1);
[FCC(:,l),DFCC(:,:,l)]=equationFCCt(z0C(i), r0C(nrC),x,pa);
end

%
% *Left (i=1)* and *right(i=NZ)* 
%
 
for j=2:nrC-1 
l = sub2ind([nrC,nzC],j,1); 
x = reshape(yfC(:,:,l)',NDC*NVC,1);
[FCC(:,l),DFCC(:,:,l)]=equationFCCl(z0C(1), r0C(j),x,pa);
 
l = sub2ind([nrC,nzC],j,nzC); 
x = reshape(yfC(:,:,l)',NDC*NVC,1);
[FCC(:,l),DFCC(:,:,l)]=equationFCCr(z0C(nzC), r0C(j),x,pa);
end

%corners (Must match with those in connections!) 


 % corner (z) i=1 and (r) j=1 

l = sub2ind([nrC,nzC],1,1); 
x = reshape(yfC(:,:,l)',NDC*NVC,1);
[FCC(:,l),DFCC(:,:,l)]=equationFCCl(z0C(1), r0C(1),x,pa);

 % corner (z) i=nrC and (r) j=1 

l = sub2ind([nrC,nzC],nrC,1); 
x = reshape(yfC(:,:,l)',NDC*NVC,1);
[FCC(:,l),DFCC(:,:,l)]=equationFCCl(z0C(1), r0C(nrC),x,pa);

 % corner (z) i=1 and (r) j=nzC 

l = sub2ind([nrC,nzC],1,nzC); 
x = reshape(yfC(:,:,l)',NDC*NVC,1);
[FCC(:,l),DFCC(:,:,l)]=equationFCCr(z0C(nzC), r0C(1),x,pa);

 % corner (z) i=nrC and (r) j=nzC 

l = sub2ind([nrC,nzC],nrC,nzC); 
x = reshape(yfC(:,:,l)',NDC*NVC,1);
[FCC(:,l),DFCC(:,:,l)]=equationFCCr(z0C(nzC), r0C(nrC),x,pa);

%
% Mounting the jacobians 
%
 
xc = 0*speye(ntC*NVC,1); 

for jj=1:NVC 
rowi = (jj-1)*ntC + 1;
rowf = jj*ntC;
for kk=1:NVC 
km=(kk-1)*NDC+1;
kp=kk*NDC;
D=squeeze(DFCC(jj,km:kp,:));;
B=spdiags(D(1,:)',0,ntC, ntC) ...
+ spdiags(D(2,:)',0, ntC, ntC)* dd0rC...
+ spdiags(D(3,:)',0, ntC, ntC)* dd0zC...
+ spdiags(D(4,:)',0, ntC, ntC)* dd0rrC...
+ spdiags(D(5,:)',0, ntC, ntC)* dd0zzC...
+ spdiags(D(6,:)',0, ntC, ntC)* dd0rzC...
+ spdiags(D(7,:)',0, ntC, ntC)* bp;
coli = (kk-1)*ntC + 1;
colf = kk*ntC;
cc(rowi:rowf,coli:colf) = sparse(B);
end
xc(rowi:rowf) = - FCC(jj,:)';
end

%
% Connections 
%
 
% Connection of C with C 
%
 
DFCC= zeros(NVC,NDC*NVC,ntC);

%
% Constructing cc 
%
 
for jj=1:NVC
list = jl +(jj-1)*ntC;
for kk=1:NVC
km=(kk-1)*NDC+1;
kp=kk*NDC;
D=squeeze(DFCC(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntC, ntC) ...
+ spdiags(D(2,:)',0, ntC, ntC)* dd0rC...
+ spdiags(D(3,:)',0, ntC, ntC)* dd0zC...
+ spdiags(D(4,:)',0, ntC, ntC)* dd0rrC...
+ spdiags(D(5,:)',0, ntC, ntC)* dd0zzC...
+ spdiags(D(6,:)',0, ntC, ntC)* dd0rzC...
+ spdiags(D(7,:)',0, ntC, ntC)* bp;
aux = sparse(B);
coli = (kk-1)*ntC + 1;
colf = kk*ntC;
cc(list,coli:colf) = aux(jl_c,:);
end
xc(list) = xc(list) - FCC(jj,il)';
end
 
% Connection of C with C 
%
 
DFCC= zeros(NVC,NDC*NVC,ntC);

%
% Constructing cc 
%
 
for jj=1:NVC
list = jl +(jj-1)*ntC;
for kk=1:NVC
km=(kk-1)*NDC+1;
kp=kk*NDC;
D=squeeze(DFCC(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntC, ntC) ...
+ spdiags(D(2,:)',0, ntC, ntC)* dd0rC...
+ spdiags(D(3,:)',0, ntC, ntC)* dd0zC...
+ spdiags(D(4,:)',0, ntC, ntC)* dd0rrC...
+ spdiags(D(5,:)',0, ntC, ntC)* dd0zzC...
+ spdiags(D(6,:)',0, ntC, ntC)* dd0rzC...
+ spdiags(D(7,:)',0, ntC, ntC)* bp;
aux = sparse(B);
coli = (kk-1)*ntC + 1;
colf = kk*ntC;
cc(list,coli:colf) = aux(jl_c,:);
end
xc(list) = xc(list) - FCC(jj,il)';
end
 

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
+ spdiags(D(2,:)',0, ntD, ntD)* dd0rD...
+ spdiags(D(3,:)',0, ntD, ntD)* dd0zD...
+ spdiags(D(4,:)',0, ntD, ntD)* dd0rrD...
+ spdiags(D(5,:)',0, ntD, ntD)* dd0zzD...
+ spdiags(D(6,:)',0, ntD, ntD)* dd0rzD...
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
 
% Connection of D with D 
%
 
DFDD= zeros(NVD,NDD*NVD,ntD);

%
% Constructing dd 
%
 
for jj=1:NVD
list = jl +(jj-1)*ntD;
for kk=1:NVD
km=(kk-1)*NDD+1;
kp=kk*NDD;
D=squeeze(DFDD(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntD, ntD) ...
+ spdiags(D(2,:)',0, ntD, ntD)* dd0rD...
+ spdiags(D(3,:)',0, ntD, ntD)* dd0zD...
+ spdiags(D(4,:)',0, ntD, ntD)* dd0rrD...
+ spdiags(D(5,:)',0, ntD, ntD)* dd0zzD...
+ spdiags(D(6,:)',0, ntD, ntD)* dd0rzD...
+ spdiags(D(7,:)',0, ntD, ntD)* bp;
aux = sparse(B);
coli = (kk-1)*ntD + 1;
colf = kk*ntD;
dd(list,coli:colf) = aux(jl_c,:);
end
xd(list) = xd(list) - FDD(jj,il)';
end
 
% Connection of D with D 
%
 
DFDD= zeros(NVD,NDD*NVD,ntD);

%
% Constructing dd 
%
 
for jj=1:NVD
list = jl +(jj-1)*ntD;
for kk=1:NVD
km=(kk-1)*NDD+1;
kp=kk*NDD;
D=squeeze(DFDD(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntD, ntD) ...
+ spdiags(D(2,:)',0, ntD, ntD)* dd0rD...
+ spdiags(D(3,:)',0, ntD, ntD)* dd0zD...
+ spdiags(D(4,:)',0, ntD, ntD)* dd0rrD...
+ spdiags(D(5,:)',0, ntD, ntD)* dd0zzD...
+ spdiags(D(6,:)',0, ntD, ntD)* dd0rzD...
+ spdiags(D(7,:)',0, ntD, ntD)* bp;
aux = sparse(B);
coli = (kk-1)*ntD + 1;
colf = kk*ntD;
dd(list,coli:colf) = aux(jl_c,:);
end
xd(list) = xd(list) - FDD(jj,il)';
end
 
% Connection of D with D 
%
 
DFDD= zeros(NVD,NDD*NVD,ntD);

%
% Constructing dd 
%
 
for jj=1:NVD
list = jl +(jj-1)*ntD;
for kk=1:NVD
km=(kk-1)*NDD+1;
kp=kk*NDD;
D=squeeze(DFDD(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntD, ntD) ...
+ spdiags(D(2,:)',0, ntD, ntD)* dd0rD...
+ spdiags(D(3,:)',0, ntD, ntD)* dd0zD...
+ spdiags(D(4,:)',0, ntD, ntD)* dd0rrD...
+ spdiags(D(5,:)',0, ntD, ntD)* dd0zzD...
+ spdiags(D(6,:)',0, ntD, ntD)* dd0rzD...
+ spdiags(D(7,:)',0, ntD, ntD)* bp;
aux = sparse(B);
coli = (kk-1)*ntD + 1;
colf = kk*ntD;
dd(list,coli:colf) = aux(jl_c,:);
end
xd(list) = xd(list) - FDD(jj,il)';
end
 
% Connection of D with D 
%
 
DFDD= zeros(NVD,NDD*NVD,ntD);

%
% Constructing dd 
%
 
for jj=1:NVD
list = jl +(jj-1)*ntD;
for kk=1:NVD
km=(kk-1)*NDD+1;
kp=kk*NDD;
D=squeeze(DFDD(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntD, ntD) ...
+ spdiags(D(2,:)',0, ntD, ntD)* dd0rD...
+ spdiags(D(3,:)',0, ntD, ntD)* dd0zD...
+ spdiags(D(4,:)',0, ntD, ntD)* dd0rrD...
+ spdiags(D(5,:)',0, ntD, ntD)* dd0zzD...
+ spdiags(D(6,:)',0, ntD, ntD)* dd0rzD...
+ spdiags(D(7,:)',0, ntD, ntD)* bp;
aux = sparse(B);
coli = (kk-1)*ntD + 1;
colf = kk*ntD;
dd(list,coli:colf) = aux(jl_c,:);
end
xd(list) = xd(list) - FDD(jj,il)';
end
 

% Block E
 
FEE= zeros(NVE,ntE);
DFEE= zeros(NVE,NDE*NVE,ntE);

 % bulk 
for i=2:nzE-1 
for j=2:nrE-1 
l=sub2ind([nrE,nzE],j,i);
x = reshape(yfE(:,:,l)',NDE*NVE,1);
[FEE(:,l),DFEE(:,:,l)]=equationFEE(z0E(i), r0E(j),x,pa);
end

%
% *top (j=NR)* and *bottom(j=1)* 
%
 
l=sub2ind([nrE,nzE],1,i); % bottom
x = reshape(yfE(:,:,l)',NDE*NVE,1);
[FEE(:,l),DFEE(:,:,l)]=equationFEEb(z0E(i), r0E(1),x,pa);

l=sub2ind([nrE,nzE],nrE,i); %top 
 x = reshape(yfE(:,:,l)',NDE*NVE,1);
[FEE(:,l),DFEE(:,:,l)]=equationFEEt(z0E(i), r0E(nrE),x,pa);
end

%
% *Left (i=1)* and *right(i=NZ)* 
%
 
for j=2:nrE-1 
l = sub2ind([nrE,nzE],j,1); 
x = reshape(yfE(:,:,l)',NDE*NVE,1);
[FEE(:,l),DFEE(:,:,l)]=equationFEEl(z0E(1), r0E(j),x,pa);
 
l = sub2ind([nrE,nzE],j,nzE); 
x = reshape(yfE(:,:,l)',NDE*NVE,1);
[FEE(:,l),DFEE(:,:,l)]=equationFEEr(z0E(nzE), r0E(j),x,pa);
end

%corners (Must match with those in connections!) 


 % corner (z) i=1 and (r) j=1 

l = sub2ind([nrE,nzE],1,1); 
x = reshape(yfE(:,:,l)',NDE*NVE,1);
[FEE(:,l),DFEE(:,:,l)]=equationFEEl(z0E(1), r0E(1),x,pa);

 % corner (z) i=nrE and (r) j=1 

l = sub2ind([nrE,nzE],nrE,1); 
x = reshape(yfE(:,:,l)',NDE*NVE,1);
[FEE(:,l),DFEE(:,:,l)]=equationFEEl(z0E(1), r0E(nrE),x,pa);

 % corner (z) i=1 and (r) j=nzE 

l = sub2ind([nrE,nzE],1,nzE); 
x = reshape(yfE(:,:,l)',NDE*NVE,1);
[FEE(:,l),DFEE(:,:,l)]=equationFEEr(z0E(nzE), r0E(1),x,pa);

 % corner (z) i=nrE and (r) j=nzE 

l = sub2ind([nrE,nzE],nrE,nzE); 
x = reshape(yfE(:,:,l)',NDE*NVE,1);
[FEE(:,l),DFEE(:,:,l)]=equationFEEr(z0E(nzE), r0E(nrE),x,pa);

%
% Mounting the jacobians 
%
 
xe = 0*speye(ntE*NVE,1); 

for jj=1:NVE 
rowi = (jj-1)*ntE + 1;
rowf = jj*ntE;
for kk=1:NVE 
km=(kk-1)*NDE+1;
kp=kk*NDE;
D=squeeze(DFEE(jj,km:kp,:));;
B=spdiags(D(1,:)',0,ntE, ntE) ...
+ spdiags(D(2,:)',0, ntE, ntE)* dd0rE...
+ spdiags(D(3,:)',0, ntE, ntE)* dd0zE...
+ spdiags(D(4,:)',0, ntE, ntE)* dd0rrE...
+ spdiags(D(5,:)',0, ntE, ntE)* dd0zzE...
+ spdiags(D(6,:)',0, ntE, ntE)* dd0rzE...
+ spdiags(D(7,:)',0, ntE, ntE)* bp;
coli = (kk-1)*ntE + 1;
colf = kk*ntE;
ee(rowi:rowf,coli:colf) = sparse(B);
end
xe(rowi:rowf) = - FEE(jj,:)';
end

%
% Connections 
%
 
% Connection of E with E 
%
 
DFEE= zeros(NVE,NDE*NVE,ntE);

%
% Constructing ee 
%
 
for jj=1:NVE
list = jl +(jj-1)*ntE;
for kk=1:NVE
km=(kk-1)*NDE+1;
kp=kk*NDE;
D=squeeze(DFEE(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntE, ntE) ...
+ spdiags(D(2,:)',0, ntE, ntE)* dd0rE...
+ spdiags(D(3,:)',0, ntE, ntE)* dd0zE...
+ spdiags(D(4,:)',0, ntE, ntE)* dd0rrE...
+ spdiags(D(5,:)',0, ntE, ntE)* dd0zzE...
+ spdiags(D(6,:)',0, ntE, ntE)* dd0rzE...
+ spdiags(D(7,:)',0, ntE, ntE)* bp;
aux = sparse(B);
coli = (kk-1)*ntE + 1;
colf = kk*ntE;
ee(list,coli:colf) = aux(jl_c,:);
end
xe(list) = xe(list) - FEE(jj,il)';
end
 
% Connection of E with E 
%
 
DFEE= zeros(NVE,NDE*NVE,ntE);

%
% Constructing ee 
%
 
for jj=1:NVE
list = jl +(jj-1)*ntE;
for kk=1:NVE
km=(kk-1)*NDE+1;
kp=kk*NDE;
D=squeeze(DFEE(jj,km:kp,:));
B = spdiags(D(1,:)',0,ntE, ntE) ...
+ spdiags(D(2,:)',0, ntE, ntE)* dd0rE...
+ spdiags(D(3,:)',0, ntE, ntE)* dd0zE...
+ spdiags(D(4,:)',0, ntE, ntE)* dd0rrE...
+ spdiags(D(5,:)',0, ntE, ntE)* dd0zzE...
+ spdiags(D(6,:)',0, ntE, ntE)* dd0rzE...
+ spdiags(D(7,:)',0, ntE, ntE)* bp;
aux = sparse(B);
coli = (kk-1)*ntE + 1;
colf = kk*ntE;
ee(list,coli:colf) = aux(jl_c,:);
end
xe(list) = xe(list) - FEE(jj,il)';
end
 
% 
% Mounting of vector b and off diagonal jacobian
% 


a = [ aa ab ac ad ae ;  ba bb bc bd be ;  ca cb cc cd ce ;  da db dc dd de ;  ea eb ec ed ee ;  ];
b = [ xa; xb; xc; xd; xe; ];
