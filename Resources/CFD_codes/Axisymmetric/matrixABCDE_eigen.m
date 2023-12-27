xt = 0*x0;
%% temporal derivatives
bp=-1i;

% Initializing the variables
% Initializing the variables

% Initializing the variables

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


inic = nbi+4*ntD+1;
ifin = nbi+5*ntD;
MzD= x0(inic:ifin);
yfD(5,1,:)=MzD;
yfD(5,2,: ) = dd0rD*MzD;
yfD(5,3,: ) = dd0zD*MzD;
yfD(5,4,: ) = dd0rrD*MzD;
yfD(5,5,: ) = dd0zzD*MzD;
yfD(5,6,: ) = dd0rzD*MzD;
yfD(5,7,: ) = xt(inic:ifin);


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









%%%%%%%%%%%%%%BLOCK A%%%%%%%%%%%%%%A+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %% Creating functions and Jacobians

%%%%%%%%%%%%%%BLOCK A%%%%%%%%%%%%%%A+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %% Creating functions and Jacobians

eval(['FAA  ' '= zeros(NVA,ntA);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,
eval(['DFAA  ' '= zeros(NVA,NVA*NDA,ntA);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,
%conections BA
eval(['FBA  ' '= zeros(NVB,ntA);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,
eval(['DFBA  ' '= zeros(NVB,NVA*NDA,ntA);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,
%conections DA
eval(['FDA  ' '= zeros(NVD,ntA);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,
eval(['DFDA  ' '= zeros(NVD,NVA*NDA,ntA);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,


%conections EA
eval(['FEA  ' '= zeros(NVE,ntA);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,
eval(['DFEA  ' '= zeros(NVE,NVA*NDA,ntA);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,

for l=1:ntA
    
    xa=reshape(yfA(:,:,l)',NVA*NDA,1); %xa=[y1(:,i,j);y2(:,i,j);y3(:,i,j);y4(:,i,j);y5(:,i,j)];
    %BULK
    [FAA(:,l),DFAA(:,:,l)]=equationFAA(z1A(l),r1A(l),xa,pa);
    
    %gettig Forcez and Forcez
    [Forcer(l),Forcez(l)]=MagneticForces(z1A(l),r1A(l),xa,pa);
    
    %gettig M and Mz
    [MA(l),MAz(l),HrA(l),HzA(l)]=MagneticfieldAA(z1A(l),r1A(l),xa,pa);
    %getting Magnetization vector Mr0 and Mz0
    [MAr0(l),MAz0(l)]=MagneticvectorAA(z1A(l),r1A(l),xa,pa);
    
    %BC:  linea 9Ab
    %for k=1:length(Linea9Ab)
    
    if(ndA(l)==9)  %
        [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea9(z1A(l),r1A(l),xa,pa);%
        [FEA(:,l),DFEA(:,:,l)]=equationFEAlinea9(z1A(l),r1A(l),xa,pa);%
    end
    if(ndA(l)==11)
        if(l==Linea11A(length(Linea11A)))
            [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea11a(z1A(l),r1A(l),xa,pa);%
        else
            [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea11(z1A(l),r1A(l),xa,pa);%
        end
        %DFAA(1,38,l)=0;
        [FBA(:,l),DFBA(:,:,l)]=equationFBAlinea11(z1A(l),r1A(l),xa,pa);%
    end
    if (ndA(l)==3)
        [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea3(z1A(l),r1A(l),xa,pa);%
        %[MA(l),MAz(l),HrA(l),HzA(l)]=MagneticfieldAAa(z1A(l),r1A(l),xa,pa);
    end
    if(ndA(l)==20)  %
        [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea20(z1A(l),r1A(l),xa,pa);%
        [FDA(:,l),DFDA(:,:,l)]=equationFDAlinea20(z1A(l),r1A(l),xa,pa);%
    end
    if(ndA(l)==145)
        [FAA(:,l),DFAA(:,:,l)]=equationFVertexADE(z1A(l),r1A(l),xa,pa);%
        %[FDA(:,l),DFDA(:,:,l)]=equationFVertexDEAA(z1A(l),r1A(l),xa,pa);%
        
    end
    
    if(ndA(l)==124)
        [FAA(:,l),DFAA(:,:,l)]=equationFVertexABD(z1A(l),r1A(l),xa,pa);%
        
    end
    
    
end



%Mouting AA
for j=1:NVA %equations
    
    for k=1:NVA  %jacobians
        km=(k-1)*NDA+1;
        kp=k*NDA;
        D=squeeze(DFAA(j,km:kp,:));
        B=spdiags(D(1,:)',0,ntA, ntA) ...
            + spdiags(D(2,:)',0, ntA, ntA)* dd0rA...
            + spdiags(D(3,:)',0, ntA, ntA)* dd0zA...
            + spdiags(D(4,:)',0, ntA, ntA)* dd0rrA...
            + spdiags(D(5,:)',0, ntA, ntA)* dd0zzA...
            + spdiags(D(6,:)',0, ntA, ntA)* dd0rzA;
        B1=+spdiags(D(7,:)',0,ntA,ntA)*bp;
        eval(['aa' num2str(j) num2str(k) '=sparse(B);'])
        eval(['caa' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end


%Mouting BA
for j=1:NVB %equations
    
    for k=1:NVA  %jacobians
        km=(k-1)*NDA+1;
        kp=k*NDA;
        
        D=squeeze(DFBA(j,km:kp,:));
        B=spdiags(D(1,:)',0,ntA, ntA) ...
            + spdiags(D(2,:)',0, ntA, ntA)* dd0rA...
            + spdiags(D(3,:)',0, ntA, ntA)* dd0zA...
            + spdiags(D(4,:)',0, ntA, ntA)* dd0rrA...
            + spdiags(D(5,:)',0, ntA, ntA)* dd0zzA...
            + spdiags(D(6,:)',0, ntA, ntA)* dd0rzA;
        B1= spdiags(D(7,:)',0, ntA, ntA)* bp;
        
        eval(['baf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cbaf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end

%Mouting DA
for j=1:NVD %equations
    for k=1:NVA  %jacobians
        km=(k-1)*NDA+1;
        kp=k*NDA;
        
        D=squeeze(DFDA(j,km:kp,:));
        B=spdiags(D(1,:)',0,ntA, ntA) ...
            + spdiags(D(2,:)',0, ntA, ntA)* dd0rA...
            + spdiags(D(3,:)',0, ntA, ntA)* dd0zA...
            + spdiags(D(4,:)',0, ntA, ntA)* dd0rrA...
            + spdiags(D(5,:)',0, ntA, ntA)* dd0zzA...
            + spdiags(D(6,:)',0, ntA, ntA)* dd0rzA;
        B1=spdiags(D(7,:)',0, ntA, ntA)* bp;
        eval(['daf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cdaf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end



%Mouting EA
for j=1:NVE %equations
    for k=1:NVA  %jacobians
        km=(k-1)*NDA+1;
        kp=k*NDA;
        
        D=squeeze(DFEA(j,km:kp,:));
        B=spdiags(D(1,:)',0,ntA, ntA) ...
            + spdiags(D(2,:)',0, ntA, ntA)* dd0rA...
            + spdiags(D(3,:)',0, ntA, ntA)* dd0zA...
            + spdiags(D(4,:)',0, ntA, ntA)* dd0rrA...
            + spdiags(D(5,:)',0, ntA, ntA)* dd0zzA...
            + spdiags(D(6,:)',0, ntA, ntA)* dd0rzA;
        B1=  spdiags(D(7,:)',0, ntA, ntA)* bp;
        eval(['eaf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['ceaf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end


%% Block 2 (BB)+++++++++++++++++++++++++++++++++++++
% %% Creating functions and Jacobians
eval(['FBB  ' '= zeros(NVB,ntB);']);
eval(['DFBB  ' '= zeros(NVB,NVB*NDB,ntB);']);
%conexions CB
eval(['FCB  ' '= zeros(NVC,ntB);']); %
eval(['DFCB  ' '= zeros(NVC,NVB*NDB,ntB);']);

%conexions AB
eval(['FAB  ' '= zeros(NVA,ntB);']); %
eval(['DFAB  ' '= zeros(NVA,NVB*NDB,ntB);']);


%conexions DB
eval(['FDB  ' '= zeros(NVD,ntB);']); %
eval(['DFDB  ' '= zeros(NVD,NVB*NDB,ntB);']);

HrB=zeros(ntB,1);
HzB=zeros(ntB,1);

for l=1:ntB
    
    xa=reshape(yfB(:,:,l)',NVB*NDB,1); %xa=[y1(:,i,j);y2(:,i,j);y3(:,i,j);y4(:,i,j);y5(:,i,j)];
    %BULK
    [FBB(:,l),DFBB(:,:,l)]=equationFBB(z1B(l),r1B(l),xa,pa);
    %BC:  linea 11B
    %gettig M and Mz
    [HrB(l),HzB(l)]=MagneticfieldBB(z1B(l),r1B(l),xa,pa);
    
    if(ndB(l)==11)  %wall
        [FBB(:,l),DFBB(:,:,l)]=equationFBBlinea11(z1B(l),r1B(l),xa,pa);%
        [FAB(:,l),DFAB(:,:,l)]=equationFABlinea11(z1B(l),r1B(l),xa,pa);%
    end
    
    
    %BC: linea 12
    
    if(ndB(l)==12)  %wall
        [FBB(:,l),DFBB(:,:,l)]=equationFBBlinea12(z1B(l),r1B(l),xa,pa);%
        [FCB(:,l),DFCB(:,:,l)]=equationFCBlinea12(z1B(l),r1B(l),xa,pa);%
    end
    
    
    %BC: linea 4
    if(ndB(l)==4)
        [FBB(:,l),DFBB(:,:,l)]=equationFBBlinea4(z1B(l),r1B(l),xa,pa);%
        %gettig M and Mz
        %[HrB(l),HzB(l)]=MagneticfieldBBa(z1B(l),r1B(l),xa,pa);
    end
    
    
    %BC: linea 4
    
    
    if(ndB(l)==10)
        [FBB(:,l),DFBB(:,:,l)]=equationFBBlinea10(z1B(l),r1B(l),xa,pa);%
        [FDB(:,l),DFDB(:,:,l)]=equationFDBlinea10(z1B(l),r1B(l),xa,pa);%
    end
    if(ndB(l)==241)
        [FBB(:,l),DFBB(:,:,l)]=equationFVertexBDA(z1B(l),r1B(l),xa,pa);%
        [FDB(:,l),DFDB(:,:,l)]=equationFVertexDABB(z1B(l),r1B(l),xa,pa);%
    end
    
    if(ndB(l)==234)
        [FBB(:,l),DFBB(:,:,l)]=equationFVertexBCD(z1B(l),r1B(l),xa,pa);%
    end
    
    
end



%2 Mounting Matrix bb

for j=1:NVB %equations
    for k=1:NVB  %jacobians
        km=(k-1)*NDB+1;
        kp=k*NDB;
        D=squeeze(DFBB(j,km:kp,:));
        B=spdiags(D(1,:)',0,ntB, ntB) ...
            + spdiags(D(2,:)',0, ntB, ntB)* dd0rB...
            + spdiags(D(3,:)',0, ntB, ntB)* dd0zB...
            + spdiags(D(4,:)',0, ntB, ntB)* dd0rrB...
            + spdiags(D(5,:)',0, ntB, ntB)* dd0zzB...
            + spdiags(D(6,:)',0, ntB, ntB)* dd0rzB;
        B1= spdiags(D(7,:)',0, ntB, ntB)* bp;
        
        eval(['bb' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cbb' num2str(j) num2str(k) '=-sparse(B1);'])
        
    end
end

%%3 Mounting Matrix abf and vector xabf
for j=1:NVA %equations
    for k=1:NVB  %jacobians
        km=(k-1)*NDB+1;
        kp=k*NDB;
        D=squeeze(DFAB(j,km:kp,:));
        B=spdiags(D(1,:)',0,ntB, ntB) ...
            + spdiags(D(2,:)',0, ntB, ntB)* dd0rB...
            + spdiags(D(3,:)',0, ntB, ntB)* dd0zB...
            + spdiags(D(4,:)',0, ntB, ntB)* dd0rrB...
            + spdiags(D(5,:)',0, ntB, ntB)* dd0zzB...
            + spdiags(D(6,:)',0, ntB, ntB)* dd0rzB;
        B1= spdiags(D(7,:)',0, ntB, ntB)* bp;
        eval(['abf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cabf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end

%%3 Mounting Matrix cbf and vector xcbf
for j=1:NVC %equations
    for k=1:NVB  %jacobians
        km=(k-1)*NDB+1;
        kp=k*NDB;
        D=squeeze(DFCB(j,km:kp,:));
        B=spdiags(D(1,:)',0,ntB, ntB) ...
            + spdiags(D(2,:)',0, ntB, ntB)* dd0rB...
            + spdiags(D(3,:)',0, ntB, ntB)* dd0zB...
            + spdiags(D(4,:)',0, ntB, ntB)* dd0rrB...
            + spdiags(D(5,:)',0, ntB, ntB)* dd0zzB...
            + spdiags(D(6,:)',0, ntB, ntB)* dd0rzB;
        B1= spdiags(D(7,:)',0, ntB, ntB)* bp;
        eval(['cbf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['ccbf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end


%%3 Mounting Matrix cbf and vector xcbf
for j=1:NVD %equations
    for k=1:NVB  %jacobians
        km=(k-1)*NDB+1;
        kp=k*NDB;
        D=squeeze(DFDB(j,km:kp,:));
        B=spdiags(D(1,:)',0,ntB, ntB) ...
            + spdiags(D(2,:)',0, ntB, ntB)* dd0rB...
            + spdiags(D(3,:)',0, ntB, ntB)* dd0zB...
            + spdiags(D(4,:)',0, ntB, ntB)* dd0rrB...
            + spdiags(D(5,:)',0, ntB, ntB)* dd0zzB...
            + spdiags(D(6,:)',0, ntB, ntB)* dd0rzB;
        B1= spdiags(D(7,:)',0, ntB, ntB)* bp;
        
        eval(['dbf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cdbf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end



%% end Block 2 (BB)+++++++++++++++++++++++++++++++++++++


%% Block 3 (CC)+++++++++++++++++++++++++++++++++++++

% %% Creating functions and Jacobians
eval(['FCC  ' '= zeros(NVC,ntC);']);
eval(['DFCC  ' '= zeros(NVC,NVC*NDC,ntC);']);
%conexions BC
eval(['FBC  ' '= zeros(NVB,ntC);']);
eval(['DFBC  ' '= zeros(NVB,NVC*NDC,ntC);']);
%conexions DC
eval(['FDC  ' '= zeros(NVD,ntC);']);
eval(['DFDC  ' '= zeros(NVD,NVC*NDC,ntC);']);
HrC=zeros(ntC,1);
HzC=zeros(ntC,1);


for l=1:ntC
    
    xa=reshape(yfC(:,:,l)',NVC*NDC,1); %xa=[y1(:,i,j);y2(:,i,j);y3(:,i,j);y4(:,i,j);y5(:,i,j)];
    %CULK
    [FCC(:,l),DFCC(:,:,l)]=equationFCC(z1C(l),r1C(l),xa,pa);
    %gettig M and Mz
    [HrC(l),HzC(l)]=MagneticfieldCC(z1C(l),r1C(l),xa,pa);
    
    %CC:  linea 12
    
    if(ndC(l)==12)
        [FCC(:,l),DFCC(:,:,l)]=equationFCClinea12(z1C(l),r1C(l),xa,pa);%
        [FBC(:,l),DFBC(:,:,l)]=equationFBClinea12(z1C(l),r1C(l),xa,pa);%
    end
    
    if(ndC(l)==6)
        if pa(19) == 0
            phid_c                 = exactphiBCs(rC(l),zC(l),pa(13),pa(25), pa(26), pa(27), pa(28));
            [FCC(:,l),DFCC(:,:,l)] = equationFCClinea6_c(z1C(l),r1C(l),xa,pa,phid_c);%
        else
            [FCC(:,l),DFCC(:,:,l)]=equationFCClinea6_m(z1C(l),r1C(l),xa,pa);%
        end
    end
    
    if(ndC(l)==5)
        [FCC(:,l),DFCC(:,:,l)]=equationFCClinea5(z1C(l),r1C(l),xa,pa);%
        %gettig M and Mz
        %[HrC(l),HzC(l)]=MagneticfieldCCa(z1C(l),r1C(l),xa,pa);
    end
    
    
    %CC: linea 13
    
    
    if(ndC(l)==13)
        [FCC(:,l),DFCC(:,:,l)]=equationFCClinea13(z1C(l),r1C(l),xa,pa);%
        [FDC(:,l),DFDC(:,:,l)]=equationFDClinea13(z1C(l),r1C(l),xa,pa);%
    end
    
    if(ndC(l)==342)
        [FCC(:,l),DFCC(:,:,l)]=equationFVertexCDB(z1C(l),r1C(l),xa,pa);%
        [FDC(:,l),DFDC(:,:,l)]=equationFVertexDBCC(z1C(l),r1C(l),xa,pa);%
    end
    
    
end
%CC
for j=1:NVC %equations
    for k=1:NVC  %jacobians
        
        
        km=(k-1)*NDC+1;
        kp=k*NDC;
        D=squeeze(DFCC(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntC, ntC) ...
            + spdiags(D(2,:)',0, ntC, ntC)* dd0rC...
            + spdiags(D(3,:)',0, ntC, ntC)* dd0zC...
            + spdiags(D(4,:)',0, ntC, ntC)* dd0rrC...
            + spdiags(D(5,:)',0, ntC, ntC)* dd0zzC...
            + spdiags(D(6,:)',0, ntC, ntC)* dd0rzC;
        B1= spdiags(D(7,:)',0, ntC, ntC)* bp;
        eval(['cc' num2str(j) num2str(k) '=sparse(B);'])
        eval(['ccc' num2str(j) num2str(k) '=-sparse(B1);'])
        
    end
end

%BC
for j=1:NVB %equations
    for k=1:NVC  %jacobians
        km=(k-1)*NDC+1;
        kp=k*NDC;
        
        D=squeeze(DFBC(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntC, ntC) ...
            + spdiags(D(2,:)',0, ntC, ntC)* dd0rC...
            + spdiags(D(3,:)',0, ntC, ntC)* dd0zC...
            + spdiags(D(4,:)',0, ntC, ntC)* dd0rrC...
            + spdiags(D(5,:)',0, ntC, ntC)* dd0zzC...
            + spdiags(D(6,:)',0, ntC, ntC)* dd0rzC;
        B1= spdiags(D(7,:)',0, ntC, ntC)* bp;
        eval(['bcf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cbcf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end


%DC
for j=1:NVD %equations
    for k=1:NVC  %jacobians
        km=(k-1)*NDC+1;
        kp=k*NDC;
        
        D=squeeze(DFDC(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntC, ntC) ...
            + spdiags(D(2,:)',0, ntC, ntC)* dd0rC...
            + spdiags(D(3,:)',0, ntC, ntC)* dd0zC...
            + spdiags(D(4,:)',0, ntC, ntC)* dd0rrC...
            + spdiags(D(5,:)',0, ntC, ntC)* dd0zzC...
            + spdiags(D(6,:)',0, ntC, ntC)* dd0rzC;
        B1= spdiags(D(7,:)',0, ntC, ntC)* bp;
        eval(['dcf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cdcf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end




%% Block 4 (DD)+++++++++++++++++++++++++++++++++++++

% %% Dreating functions and Jacobians
eval(['FDD  ' '= zeros(NVD,ntD);']);
eval(['DFDD  ' '= zeros(NVD,NVD*NDD,ntD);']);
%conexions AD
eval(['FAD  ' '= zeros(NVA,ntD);']);
eval(['DFAD  ' '= zeros(NVA,NVD*NDD,ntD);']);
%conexions BD
eval(['FBD  ' '= zeros(NVB,ntD);']);
eval(['DFBD  ' '= zeros(NVB,NVD*NDD,ntD);']);
%conexions CD
eval(['FCD  ' '= zeros(NVC,ntD);']);
eval(['DFCD  ' '= zeros(NVC,NVD*NDD,ntD);']);

%conexions ED
eval(['FED  ' '= zeros(NVE,ntD);']);
eval(['DFED  ' '= zeros(NVE,NVD*NDD,ntD);']);
HrD=zeros(ntD,1);
HzD=zeros(ntD,1);
for l=1:ntD
    xa=reshape(yfD(:,:,l)',NVD*NDD,1); %xa=[y1(:,i,j);y2(:,i,j);y3(:,i,j);y4(:,i,j);y5(:,i,j)];
    
    %BULK
    pa(23)=Fcoil(l);
    pa(24)=Gcoil(l);
    [FDD(:,l),DFDD(:,:,l)]=equationFDD(z1D(l),r1D(l),xa,pa);
    
    %gettig M and Mz
    [HrD(l),HzD(l)]=MagneticfieldDD(z1D(l),r1D(l),xa,pa);
    [dphidrD(l),dphidzD(l),dphidrrD(l),dphidzzD(l),dphidrzD(l),dpsidrD(l),dpsidzD(l),...
        dpsidrrD(l),dpsidzzD(l),dpsidrzD(l)] = potentialsDD(z1D(l),r1D(l),xa,pa);
    
    if(ndD(l)==33)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDCoil(z1D(l),r1D(l),xa,pa);%
    end
    if(ndD(l)==34)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDCoil1(z1D(l),r1D(l),xa,pa);%
    end
    if(ndD(l)==35)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDRightCoil(z1D(l),r1D(l),xa,pa);%
    end
    if(ndD(l)==36)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDCoilContour(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==21)
        if pa(19) ==0 % CHECK IF WE ARE USING A MAGNET OR A COIL
            phid_c                 = exactphiBCs(rD(l),zD(l),pa(13),pa(25), pa(26), pa(27), pa(28));
            [FDD(:,l),DFDD(:,:,l)] = equationFDDlinea21_c(z1D(l),r1D(l),xa,pa,phid_c);%
            if (l==Linea8D(2))
                [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea8_center_c(zD(l),rD(l),xa,pa,phid_c);%
            end
        else
            [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea21_m(z1D(l),r1D(l),xa,pa);%
            
            if (l==Linea8D(2))
                [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea8_center_m(zD(l),rD(l),xa,pa);%
            end
            
        end
    end
    
    
    if(ndD(l)==8)
        if pa(19) ==0  % CHECK IF WE ARE USING A MAGNET OR A COIL
            phid_c                 = exactphiBCs(rD(l),zD(l),pa(13),pa(25), pa(26), pa(27), pa(28));
            [FDD(:,l),DFDD(:,:,l)] = equationFDDlinea8_c(z1D(l),r1D(l),xa,pa,phid_c);%
        else
            [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea8_m(z1D(l),r1D(l),xa,pa);%
        end
        
    end
    
    
    
    if(ndD(l)==7)
        if pa(19) ==0 % CHECK IF WE ARE USING A MAGNET OR A COIL
            phid_c                 = exactphiBCs(rD(l),zD(l),pa(13),pa(25), pa(26), pa(27), pa(28));
            [FDD(:,l),DFDD(:,:,l)] = equationFDDlinea7_c(z1D(l),r1D(l),xa,pa,phid_c);%
        else
            [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea7_m(z1D(l),r1D(l),xa,pa);%
        end
    end
    
    if(ndD(l)==19)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea19(z1D(l),r1D(l),xa,pa);%
        [FED(:,l),DFED(:,:,l)]=equationFEDlinea19(z1D(l),r1D(l),xa,pa);%
        
    end
    
    
    
    if(ndD(l)==20)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea20(z1D(l),r1D(l),xa,pa);%
        [FAD(:,l),DFAD(:,:,l)]=equationFADlinea20(z1D(l),r1D(l),xa,pa);%
        
    end
    
    
    if(ndD(l)==10)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea10(z1D(l),r1D(l),xa,pa);%
        [FBD(:,l),DFBD(:,:,l)]=equationFBDlinea10(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==451)
        [FDD(:,l),DFDD(:,:,l)]=equationFVertexDEA(z1D(l),r1D(l),xa,pa);%
        [FAD(:,l),DFAD(:,:,l)]=equationFVertexADED(z1D(l),r1D(l),xa,pa);%
        [FED(:,l),DFED(:,:,l)]=equationFVertexEADD(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==412)
        [FDD(:,l),DFDD(:,:,l)]=equationFVertexDAB(z1D(l),r1D(l),xa,pa);%
        [FAD(:,l),DFAD(:,:,l)]=equationFVertexABDD(z1D(l),r1D(l),xa,pa);%
        [FBD(:,l),DFBD(:,:,l)]=equationFVertexBDAD(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==423)
        [FDD(:,l),DFDD(:,:,l)]=equationFVertexDBC(z1D(l),r1D(l),xa,pa);%
        [FCD(:,l),DFCD(:,:,l)]=equationFVertexCDBD(z1D(l),r1D(l),xa,pa);%
        [FBD(:,l),DFBD(:,:,l)]=equationFVertexBCDD(z1D(l),r1D(l),xa,pa);%
    end
    
    
    
    if(ndD(l)==13)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea13(z1D(l),r1D(l),xa,pa);%
        [FCD(:,l),DFCD(:,:,l)]=equationFCDlinea13(z1D(l),r1D(l),xa,pa);%
    end
    
    
    
    
    
end

%DD
for j=1:NVD %equations
    for k=1:NVD  %jacobians
        
        
        km=(k-1)*NDD+1;
        kp=k*NDD;
        D=squeeze(DFDD(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntD, ntD) ...
            + spdiags(D(2,:)',0, ntD, ntD)* dd0rD...
            + spdiags(D(3,:)',0, ntD, ntD)* dd0zD...
            +spdiags(D(4,:)',0, ntD, ntD)* dd0rrD...
            + spdiags(D(5,:)',0, ntD, ntD)* dd0zzD...
            + spdiags(D(6,:)',0, ntD, ntD)* dd0rzD;
        B1=+ spdiags(D(7,:)',0, ntD, ntD)* bp;
        eval(['dd' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cdd' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end

%AD
for j=1:NVA %equations
    for k=1:NVD  %jacobians
        km=(k-1)*NDD+1;
        kp=k*NDD;
        
        D=squeeze(DFAD(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntD, ntD) ...
            + spdiags(D(2,:)',0, ntD, ntD)* dd0rD...
            + spdiags(D(3,:)',0, ntD, ntD)* dd0zD...
            +spdiags(D(4,:)',0, ntD, ntD)* dd0rrD...
            + spdiags(D(5,:)',0, ntD, ntD)* dd0zzD...
            + spdiags(D(6,:)',0, ntD, ntD)* dd0rzD;
        B1=+ spdiags(D(7,:)',0, ntD, ntD)* bp;
        eval(['adf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cadf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end
%BD
for j=1:NVB %equations
    for k=1:NVD  %jacobians
        km=(k-1)*NDD+1;
        kp=k*NDD;
        
        D=squeeze(DFBD(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntD, ntD) ...
            + spdiags(D(2,:)',0, ntD, ntD)* dd0rD...
            + spdiags(D(3,:)',0, ntD, ntD)* dd0zD...
            +spdiags(D(4,:)',0, ntD, ntD)* dd0rrD...
            + spdiags(D(5,:)',0, ntD, ntD)* dd0zzD...
            + spdiags(D(6,:)',0, ntD, ntD)* dd0rzD;
        B1=+ spdiags(D(7,:)',0, ntD, ntD)* bp;
        eval(['bdf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cbdf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end

%CD
for j=1:NVC %equations
    for k=1:NVD  %jacobians
        km=(k-1)*NDD+1;
        kp=k*NDD;
        
        D=squeeze(DFCD(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntD, ntD) ...
            + spdiags(D(2,:)',0, ntD, ntD)* dd0rD...
            + spdiags(D(3,:)',0, ntD, ntD)* dd0zD...
            +spdiags(D(4,:)',0, ntD, ntD)* dd0rrD...
            + spdiags(D(5,:)',0, ntD, ntD)* dd0zzD...
            + spdiags(D(6,:)',0, ntD, ntD)* dd0rzD;
        B1=+ spdiags(D(7,:)',0, ntD, ntD)* bp;
        eval(['cdf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['ccdf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end

%ED
for j=1:NVE%equations
    for k=1:NVD  %jacobians
        km=(k-1)*NDD+1;
        kp=k*NDD;
        
        D=squeeze(DFED(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntD, ntD) ...
            + spdiags(D(2,:)',0, ntD, ntD)* dd0rD...
            + spdiags(D(3,:)',0, ntD, ntD)* dd0zD...
            +spdiags(D(4,:)',0, ntD, ntD)* dd0rrD...
            + spdiags(D(5,:)',0, ntD, ntD)* dd0zzD...
            + spdiags(D(6,:)',0, ntD, ntD)* dd0rzD;
        B1=+ spdiags(D(7,:)',0, ntD, ntD)* bp;
        eval(['edf' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cedf' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end



%% Block 5 (EE)+++++++++++++++++++++++++++++++++++++

% %% Dreating functions and Jacobians
eval(['FEE  ' '= zeros(NVE,ntE);']);
eval(['DFEE  ' '= zeros(NVE,NVE*NDE,ntE);']);
%conexions DE
eval(['FDE  ' '= zeros(NVD,ntE);']);
eval(['DFDE  ' '= zeros(NVD,NVE*NDE,ntE);']);

%conexions AE
eval(['FAE  ' '= zeros(NVA,ntE);']);
eval(['DFAE  ' '= zeros(NVA,NVE*NDE,ntE);']);


HrE=zeros(ntE,1);
HzE=zeros(ntE,1);

% Virtual Phi
phi_ref = exactphiBCs(rE,zE,pa(13),pa(25), pa(26), pa(27), pa(28));


for l=1:ntE
    xa=reshape(yfE(:,:,l)',NVE*NDE,1); %xa=[y1(:,i,j);y2(:,i,j);y3(:,i,j);y4(:,i,j);y5(:,i,j)];
    %DULK
    [FEE(:,l),DFEE(:,:,l)]=equationFEE(z1E(l),r1E(l),xa,pa,phi_ref(l));
    %DD:  linea 9b
    [HrE(l),HzE(l)]=MagneticfieldEE(z1E(l),r1E(l),xa,pa);
    
    if (ndE(l)==2) %axis
        [FEE(:,l),DFEE(:,:,l)]=equationFEElinea2(zE(l),rE(l),xa,pa);%
        %[HrE(l),HzE(l)]=MagneticfieldEEa(z1E(l),r1E(l),xa,pa);
    end
    
    if(ndE(l)==1)
        if pa(19) ==0  % CHECK IF WE ARE USING A MAGNET OR A COIL
            phid_c                 = exactphiBCs(rE(l),zE(l),pa(13),pa(25), pa(26), pa(27), pa(28));
            [FEE(:,l),DFEE(:,:,l)] = equationFEElinea1_c(z1E(l),r1E(l),xa,pa,phid_c);%
        else
            [FEE(:,l),DFEE(:,:,l)]=equationFEElinea1_m(z1E(l),r1E(l),xa,pa);%
        end
        
    end
    if(ndE(l)==9)
        [FEE(:,l),DFEE(:,:,l)]=equationFEElinea9(z1E(l),r1E(l),xa,pa);%
        [FAE(:,l),DFAE(:,:,l)]=equationFAElinea9(z1E(l),r1E(l),xa,pa);%
    end
    if(ndE(l)==19)
        [FEE(:,l),DFEE(:,:,l)]=equationFEElinea19(z1E(l),r1E(l),xa,pa);%
        [FDE(:,l),DFDE(:,:,l)]=equationFDElinea19(z1E(l),r1E(l),xa,pa);%
    end
    
    if(ndE(l)==514)
        [FEE(:,l),DFEE(:,:,l)]=equationVertexEAD(z1E(l),r1E(l),xa,pa);%
        [FDE(:,l),DFDE(:,:,l)]=equationFVertexDEAE(z1E(l),r1E(l),xa,pa);%
    end
    
    
end



%EE
for j=1:NVE %equations
    for k=1:NVE  %jacobians
        
        
        km=(k-1)*NDE+1;
        kp=k*NDE;
        D=squeeze(DFEE(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntE, ntE) ...
            + spdiags(D(2,:)',0, ntE, ntE)* dd0rE...
            + spdiags(D(3,:)',0, ntE, ntE)* dd0zE...
            +spdiags(D(4,:)',0, ntE, ntE)* dd0rrE...
            + spdiags(D(5,:)',0, ntE, ntE)* dd0zzE...
            + spdiags(D(6,:)',0, ntE, ntE)* dd0rzE;
        B1=spdiags(D(7,:)',0, ntE, ntE)* bp;
        eval(['ee' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cee' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end

%DE
for j=1:NVD %equations
    for k=1:NVE  %jacobians
        
        
        km=(k-1)*NDE+1;
        kp=k*NDE;
        D=squeeze(DFDE(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntE, ntE) ...
            + spdiags(D(2,:)',0, ntE, ntE)* dd0rE...
            + spdiags(D(3,:)',0, ntE, ntE)* dd0zE...
            +spdiags(D(4,:)',0, ntE, ntE)* dd0rrE...
            + spdiags(D(5,:)',0, ntE, ntE)* dd0zzE...
            + spdiags(D(6,:)',0, ntE, ntE)* dd0rzE;
        B1=spdiags(D(7,:)',0, ntE, ntE)* bp;
        eval(['def' num2str(j) num2str(k) '=sparse(B);'])
        eval(['cdef' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end

%AE
for j=1:NVA %equations
    
    for k=1:NVE  %jacobians
        
        
        km=(k-1)*NDE+1;
        kp=k*NDE;
        D=squeeze(DFAE(j,km:kp,:));
        B = spdiags(D(1,:)',0,ntE, ntE) ...
            + spdiags(D(2,:)',0, ntE, ntE)* dd0rE...
            + spdiags(D(3,:)',0, ntE, ntE)* dd0zE...
            +spdiags(D(4,:)',0, ntE, ntE)* dd0rrE...
            + spdiags(D(5,:)',0, ntE, ntE)* dd0zzE...
            + spdiags(D(6,:)',0, ntE, ntE)* dd0rzE;
        B1=spdiags(D(7,:)',0, ntE, ntE)* bp;
        eval(['aef' num2str(j) num2str(k) '=sparse(B);'])
        eval(['caef' num2str(j) num2str(k) '=-sparse(B1);'])
    end
end








%conexions
%% Creating block  2-3  AB

for i = 1:NVA %equation
    for j=1:NVB %variable
        eval(['ab' num2str(i) num2str(j) ' = 0*speye(ntA,ntB);']);
        eval(['cab' num2str(i) num2str(j) ' = 0*speye(ntA,ntB);']);
    end
end
%% Creating block   AC
for i = 1:NVA %equation
    for j=1:NVC %variable
        eval(['ac' num2str(i) num2str(j) ' = 0*speye(ntA,ntC);']);
        eval(['cac' num2str(i) num2str(j) ' = 0*speye(ntA,ntC);']);
    end
end
%% Creating block   AD
for i = 1:NVA %equation
    for j=1:NVD %variable
        eval(['ad' num2str(i) num2str(j) ' = 0*speye(ntA,ntD);']);
        eval(['cad' num2str(i) num2str(j) ' = 0*speye(ntA,ntD);']);
        
    end
end

%% Creating block   AE
for i = 1:NVA %equation
    for j=1:NVE %variable
        eval(['ae' num2str(i) num2str(j) ' = 0*speye(ntA,ntE);']);
        eval(['cae' num2str(i) num2str(j) ' = 0*speye(ntA,ntE);']);
    end
end
%% Creating block   BA
for i = 1:NVB %equation
    for j=1:NVA %variable
        eval(['ba' num2str(i) num2str(j) ' = 0*speye(ntB,ntA);']);
        eval(['cba' num2str(i) num2str(j) ' = 0*speye(ntB,ntA);']);
    end
end
%% Creating block  2-3  BC
for i = 1:NVB %equation
    for j=1:NVC   %variable
        eval(['bc' num2str(i) num2str(j) ' = 0*speye(ntB,ntC);']);
        eval(['cbc' num2str(i) num2str(j) ' = 0*speye(ntB,ntC);']);
    end
end
%% Creating block  2-3  BD
for i = 1:NVB %equation
    for j=1:NVD   %variable
        eval(['bd' num2str(i) num2str(j) ' = 0*speye(ntB,ntD);']);
        eval(['cbd' num2str(i) num2str(j) ' = 0*speye(ntB,ntD);']);
    end
end
%% Creating block  2-3  BE
for i = 1:NVB %equation
    for j=1:NVE   %variable
        eval(['be' num2str(i) num2str(j) ' = 0*speye(ntB,ntE);']);
        eval(['cbe' num2str(i) num2str(j) ' = 0*speye(ntB,ntE);']);
    end
end

%% Creating block  3-3  CA
for i = 1:NVC %equation
    for j=1:NVA   %variable
        eval(['ca' num2str(i) num2str(j) ' = 0*speye(ntC,ntA);']);
        eval(['cca' num2str(i) num2str(j) ' = 0*speye(ntC,ntA);']);
    end
end

%% Creating block  3-3  CB
for i = 1:NVC %equation
    for j=1:NVB   %variable
        eval(['cb' num2str(i) num2str(j) ' = 0*speye(ntC,ntB);']);
        eval(['ccb' num2str(i) num2str(j) ' = 0*speye(ntC,ntB);']);
    end
end

%% Creating block  3-3  CD
for i = 1:NVC %equation
    for j=1:NVD   %variable
        eval(['cd' num2str(i) num2str(j) ' = 0*speye(ntC,ntD);']);
        eval(['ccd' num2str(i) num2str(j) ' = 0*speye(ntC,ntD);']);
    end
end

%% Creating block  3-3  CE
for i = 1:NVC %equation
    for j=1:NVE   %variable
        eval(['ce' num2str(i) num2str(j) ' = 0*speye(ntC,ntE);']);
        eval(['cce' num2str(i) num2str(j) ' = 0*speye(ntC,ntE);']);
    end
end

%% Creating block  3-3  DA
for i = 1:NVD %equation
    for j=1:NVA   %variable
        eval(['da' num2str(i) num2str(j) ' = 0*speye(ntD,ntA);']);
        eval(['cda' num2str(i) num2str(j) ' = 0*speye(ntD,ntA);']);
    end
end

%% Creating block  3-3  DB
for i = 1:NVD %equation
    for j=1:NVB   %variable
        eval(['db' num2str(i) num2str(j) ' = 0*speye(ntD,ntB);']);
        eval(['cdb' num2str(i) num2str(j) ' = 0*speye(ntD,ntB);']);
    end
end

%% Creating block  3-3  DC
for i = 1:NVD %equation
    for j=1:NVC   %variable
        eval(['dc' num2str(i) num2str(j) ' = 0*speye(ntD,ntC);']);
        eval(['cdc' num2str(i) num2str(j) ' = 0*speye(ntD,ntC);']);
    end
end

%% Creating block  3-3  DE
for i = 1:NVD %equation
    for j=1:NVE   %variable
        eval(['de' num2str(i) num2str(j) ' = 0*speye(ntD,ntE);']);
        eval(['cde' num2str(i) num2str(j) ' = 0*speye(ntD,ntE);']);
    end
end

%% Creating block  3-3  EA
for i = 1:NVE %equation
    for j=1:NVA   %variable
        eval(['ea' num2str(i) num2str(j) ' = 0*speye(ntE,ntA);']);
        eval(['cea' num2str(i) num2str(j) ' = 0*speye(ntE,ntA);']);
    end
end

%% Creating block  3-3  EB
for i = 1:NVE %equation
    for j=1:NVB   %variable
        eval(['eb' num2str(i) num2str(j) ' = 0*speye(ntE,ntB);']);
        eval(['ceb' num2str(i) num2str(j) ' = 0*speye(ntE,ntB);']);
    end
end

%% Creating block  3-3  EC
for i = 1:NVE %equation
    for j=1:NVC   %variable
        eval(['ec' num2str(i) num2str(j) ' = 0*speye(ntE,ntC);']);
        eval(['cec' num2str(i) num2str(j) ' = 0*speye(ntE,ntC);']);
    end
end

%% Creating block  3-3  ED
for i = 1:NVE %equation
    for j=1:NVD   %variable
        eval(['ed' num2str(i) num2str(j) ' = 0*speye(ntE,ntD);']);
        eval(['ced' num2str(i) num2str(j) ' = 0*speye(ntE,ntD);']);
    end
end







%correcting blocks and functions

%conexting?
iconect=1
if( iconect==1)
    
    %AB+BA++++++++++++++
    js1=[Linea11A,VertexABD];
    js2=[Linea11B,VertexBDA];
    for ks=1:NVA
        for js=1:NVB
            eval(['ab' num2str(ks) num2str(js) '(js1,:)' '=' 'abf' num2str(ks) num2str(js) '(js2,:);'])
            eval(['cab' num2str(ks) num2str(js) '(js1,:)' '=' 'cabf' num2str(ks) num2str(js) '(js2,:);'])
        end
    end
    
    
    for ks=1:NVB
        for js=1:NVA
            eval(['ba' num2str(ks) num2str(js) '(js2,:)' '=' 'baf' num2str(ks) num2str(js) '(js1,:);'])
            eval(['cba' num2str(ks) num2str(js) '(js2,:)' '=' 'cbaf' num2str(ks) num2str(js) '(js1,:);'])
        end
    end
    
    %AD+DA++++++++++++++
    js1=[Linea20A,VertexABD,VertexADE];
    js2=[Linea20D,VertexDAB,VertexDEA];
    for ks=1:NVA
        for js=1:NVD
            eval(['ad' num2str(ks) num2str(js) '(js1,:)' '=' 'adf' num2str(ks) num2str(js) '(js2,:);'])
            eval(['cad' num2str(ks) num2str(js) '(js1,:)' '=' 'cadf' num2str(ks) num2str(js) '(js2,:);'])
        end
    end
    
    for ks=1:NVD
        for js=1:NVA
            eval(['da' num2str(ks) num2str(js) '(js2,:)' '=' 'daf' num2str(ks) num2str(js) '(js1,:);'])
            eval(['cda' num2str(ks) num2str(js) '(js2,:)' '=' 'cdaf' num2str(ks) num2str(js) '(js1,:);'])
        end
    end
    
    %AE+EA++++++++++++++
    js1=[Linea9A,VertexADE];
    js2=[Linea9E,VertexEAD];
    for ks=1:NVA
        for js=1:NVE
            eval(['ae' num2str(ks) num2str(js) '(js1,:)' '=' 'aef' num2str(ks) num2str(js) '(js2,:);'])
            eval(['cae' num2str(ks) num2str(js) '(js1,:)' '=' 'caef' num2str(ks) num2str(js) '(js2,:);'])
        end
    end
    
    for ks=1:NVE
        for js=1:NVA
            eval(['ea' num2str(ks) num2str(js) '(js2,:)' '=' 'eaf' num2str(ks) num2str(js) '(js1,:);'])
            eval(['cea' num2str(ks) num2str(js) '(js2,:)' '=' 'ceaf' num2str(ks) num2str(js) '(js1,:);'])
        end
    end
    
    
    
    
    
    %CB+BC++++++++++++++
    js1=[Linea12C,VertexCDB];
    js2=[Linea12B,VertexBCD];
    for ks=1:NVC
        for js=1:NVB
            eval(['cb' num2str(ks) num2str(js) '(js1,:)' '=' 'cbf' num2str(ks) num2str(js) '(js2,:);'])
            eval(['ccb' num2str(ks) num2str(js) '(js1,:)' '=' 'ccbf' num2str(ks) num2str(js) '(js2,:);'])
        end
    end
    
    
    for ks=1:NVB
        for js=1:NVC
            eval(['bc' num2str(ks) num2str(js) '(js2,:)' '=' 'bcf' num2str(ks) num2str(js) '(js1,:);'])
            eval(['cbc' num2str(ks) num2str(js) '(js2,:)' '=' 'cbcf' num2str(ks) num2str(js) '(js1,:);'])
        end
    end
    
    
    
    
    
    
    
    
    
    
    %DB+BD++++++++++++++
    js1=[Linea10D,VertexDAB,VertexDBC];
    js2=[Linea10B,VertexBDA,VertexBCD];
    for ks=1:NVD
        for js=1:NVB
            eval(['db' num2str(ks) num2str(js) '(js1,:)' '=' 'dbf' num2str(ks) num2str(js) '(js2,:);'])
            eval(['cdb' num2str(ks) num2str(js) '(js1,:)' '=' 'cdbf' num2str(ks) num2str(js) '(js2,:);'])
        end
    end
    
    for ks=1:NVB
        for js=1:NVD
            eval(['bd' num2str(ks) num2str(js) '(js2,:)' '=' 'bdf' num2str(ks) num2str(js) '(js1,:);'])
            eval(['cbd' num2str(ks) num2str(js) '(js2,:)' '=' 'cbdf' num2str(ks) num2str(js) '(js1,:);'])
        end
    end
    
    
    
    %DC+CD++++++++++++++
    js2=[Linea13C,VertexCDB];
    js1=[Linea13D,VertexDBC];
    for ks=1:NVD
        for js=1:NVC
            eval(['dc' num2str(ks) num2str(js) '(js1,:)' '=' 'dcf' num2str(ks) num2str(js) '(js2,:);'])
            eval(['cdc' num2str(ks) num2str(js) '(js1,:)' '=' 'cdcf' num2str(ks) num2str(js) '(js2,:);'])
        end
    end
    
    
    
    for ks=1:NVC
        for js=1:NVD
            eval(['cd' num2str(ks) num2str(js) '(js2,:)' '=' 'cdf' num2str(ks) num2str(js) '(js1,:);'])
            eval(['ccd' num2str(ks) num2str(js) '(js2,:)' '=' 'ccdf' num2str(ks) num2str(js) '(js1,:);'])
        end
    end
    
    
    
    
    
    
    
    
    %
    %
    %DE+ED++++++++++++++
    js1=[Linea19D,VertexDEA];
    js2=[Linea19E,VertexEAD];
    for ks=1:NVD
        for js=1:NVE
            eval(['de' num2str(ks) num2str(js) '(js1,:)' '=' 'def' num2str(ks) num2str(js) '(js2,:);'])
            eval(['cde' num2str(ks) num2str(js) '(js1,:)' '=' 'cdef' num2str(ks) num2str(js) '(js2,:);'])
        end
    end
    
    
    for ks=1:NVE
        for js=1:NVD
            eval(['ed' num2str(ks) num2str(js) '(js2,:)' '=' 'edf' num2str(ks) num2str(js) '(js1,:);'])
            %% eval(['ced' num2str(ks) num2str(js) '(js2,:)' '=' 'edf' num2str(ks) num2str(js) '(js1,:);']) MALLLLL
            eval(['ced' num2str(ks) num2str(js) '(js2,:)' '=' 'cedf' num2str(ks) num2str(js) '(js1,:);']) %BIEN
            
        end
    end
    
    
    
    
    
end








%BLOCK A
caa=[];
aa=[];
for i=1:NVA
    ax=[];
    ay=[];
    for j=1:NVA
        aa0=eval(cat(2,'aa',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'caa',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    aa=[aa;ax];
    caa=[caa;ay];
end

%BLOCK B
cbb=[];
bb=[];
for i=1:NVB
    ax=[];
    ay=[];
    for j=1:NVB
        aa0=eval(cat(2,'bb',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cbb',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    bb=[bb;ax];
    cbb=[cbb;ay];
end




%BLOCK C;
ccc=[];
cc=[];
for i=1:NVC
    ax=[];
    ay=[];
    for j=1:NVC
        aa0=eval(cat(2,'cc',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'ccc',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    cc=[cc;ax];
    ccc=[ccc;ay];
end

%BLOCK D;
cdd=[];
dd=[];
for i=1:NVD
    ax=[];
    ay=[];
    for j=1:NVD
        aa0=eval(cat(2,'dd',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cdd',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    
    dd=[dd;ax];
    cdd=[cdd;ay];
end


%BLOCK E;
cee=[];
ee=[];
for i=1:NVE
    ax=[];
    ay=[];
    for j=1:NVE
        aa0=eval(cat(2,'ee',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cee',num2str(i),num2str(j)));
        ay=[ay,aa1];
        
    end
    
    ee=[ee;ax];
    cee=[cee;ay];
end



%Conexions
%BLOCK AB;
ab=[];
cab=[];
for i=1:NVA
    ax=[];
    ay=[];
    for j=1:NVB
        aa0=eval(cat(2,'ab',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cab',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    ab=[ab;ax];
    cab=[cab;ay];
end

%BLOCK AC;
ac=[];
cac=[];
for i=1:NVA
    ax=[];
    ay=[];
    for j=1:NVC
        aa0=eval(cat(2,'ac',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cac',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    ac=[ac;ax];
    cac=[cac;ay];
end
%BLOCK AD;
ad=[];
cad=[];
for i=1:NVA
    ax=[];
    ay=[];
    for j=1:NVD
        aa0=eval(cat(2,'ad',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cad',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    ad=[ad;ax];
    cad=[cad;ay];
end
%BLOCK AE;
ae=[];
cae=[];
for i=1:NVA
    ax=[];
    ay=[];
    for j=1:NVE
        aa0=eval(cat(2,'ae',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cae',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    ae=[ae;ax];
    cae=[cae;ay];
end






%BLOCK BA;
ba=[];
cba=[];
for i=1:NVB
    ax=[];
    ay=[];
    for j=1:NVA
        aa0=eval(cat(2,'ba',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cba',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    ba=[ba;ax];
    cba=[cba;ay];
end

%BLOCK BC;
bc=[];
cbc=[];
for i=1:NVB
    ax=[];
    ay=[];
    for j=1:NVC
        aa0=eval(cat(2,'bc',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cbc',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    bc=[bc;ax];
    cbc=[cbc;ay];
end

%BLOCK BD;
bd=[];
cbd=[];
for i=1:NVB
    ax=[];
    ay=[];
    for j=1:NVD
        aa0=eval(cat(2,'bd',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cbd',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    bd=[bd;ax];
    cbd=[cbd;ay];
end

%BLOCK BE;
be=[];
cbe=[];
for i=1:NVB
    ax=[];
    ay=[];
    for j=1:NVE
        aa0=eval(cat(2,'be',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cbe',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    be=[be;ax];
    cbe=[cbe;ay];
end



%BLOCK CA;
ca=[];
cca=[];
for i=1:NVC
    ax=[];
    ay=[];
    for j=1:NVA
        aa0=eval(cat(2,'ca',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cca',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    ca=[ca;ax];
    cca=[cca;ay];
end



%BLOCK CB;
cb=[];
ccb=[];
for i=1:NVC
    ax=[];
    ay=[];
    for j=1:NVB
        aa0=eval(cat(2,'cb',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'ccb',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    cb=[cb;ax];
    ccb=[ccb;ay];
end


%BLOCK CD;
cd=[];
ccd=[];
for i=1:NVC
    ax=[];
    ay=[];
    for j=1:NVD
        aa0=eval(cat(2,'cd',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'ccd',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    cd=[cd;ax];
    ccd=[ccd;ay];
end
%BLOCK CE;
ce=[];
cce=[];
for i=1:NVC
    ax=[];
    ay=[];
    for j=1:NVE
        aa0=eval(cat(2,'ce',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cce',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    ce=[ce;ax];
    cce=[cce;ay];
end






%BLOCK DA;
da=[];
cda=[];
for i=1:NVD
    ax=[];
    ay=[];
    for j=1:NVA
        aa0=eval(cat(2,'da',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cda',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    da=[da;ax];
    cda=[cda;ay];
end



%BLOCK DB;
db=[];
cdb=[];
for i=1:NVD
    ax=[];
    ay=[];
    for j=1:NVB
        aa0=eval(cat(2,'db',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cdb',num2str(i),num2str(j)));
        ay=[ay,aa1];
        
    end
    db=[db;ax];
    cdb=[cdb;ay];
end


%BLOCK DC;
dc=[];
cdc=[];
for i=1:NVD
    ax=[];
    ay=[];
    for j=1:NVC
        aa0=eval(cat(2,'dc',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cdc',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    dc=[dc;ax];
    cdc=[cdc;ay];
end

%BLOCK DE;
de=[];
cde=[];
for i=1:NVD
    ax=[];
    ay=[];
    for j=1:NVE
        aa0=eval(cat(2,'de',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cde',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    de=[de;ax];
    cde=[cde;ay];
end


%BLOCK EA;
ea=[];
cea=[];
for i=1:NVE
    ax=[];
    ay=[];
    for j=1:NVA
        aa0=eval(cat(2,'ea',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cea',num2str(i),num2str(j)));
        ay=[ay,aa1];
        
    end
    ea=[ea;ax];
    cea=[cea;ay];
end



%BLOCK EB;
eb=[];
ceb=[];
for i=1:NVE
    ax=[];
    ay=[];
    for j=1:NVB
        aa0=eval(cat(2,'eb',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'ceb',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    eb=[eb;ax];
    ceb=[ceb;ay];
end


%BLOCK EC;
ec=[];
cec=[];
for i=1:NVE
    ax=[];
    ay=[];
    for j=1:NVC
        aa0=eval(cat(2,'ec',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'cec',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    ec=[ec;ax];
    cec=[cec;ay];
end

%BLOCK ED;
ed=[];
ced=[];
for i=1:NVE
    ax=[];
    ay=[];
    for j=1:NVD
        aa0=eval(cat(2,'ed',num2str(i),num2str(j)));
        ax=[ax,aa0];
        aa1=eval(cat(2,'ced',num2str(i),num2str(j)));
        ay=[ay,aa1];
    end
    ed=[ed;ax];
    ced=[ced;ay];
    
end



%MATRX A and vector B
a=[aa,ab,ac,ad,ae;ba,bb,bc,bd,be;ca,cb,cc,cd,ce;da,db,dc,dd,de;ea,eb,ec,ed,ee];
b=[caa,cab,cac,cad,cae;cba,cbb,cbc,cbd,cbe;cca,ccb,ccc,ccd,cce;cda,cdb,cdc,cdd,cde;cea,ceb,cec,ced,cee];

% a=aa;
% b=caa;
% %
% a=[aa,ab,ac;ba,bb,bc;ca,cb,cc];
% b=([xa;xb;xc]);

% a=[aa,ab,ac;ba,bb];
% b=([xa;xb]);
%  a=[aa,ab,ac,ad;ba,bb,bc,bd;ca,cb,cc,cd;da,db,dc,dd];
%  b=([xa;xb;xc;xd]);

%  a=[aa];
%  b=[caa];

