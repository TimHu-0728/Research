% Initializing the variables

nbi = 0;
yfD= zeros(NVD,NDD,ntD);
xt=0*x0;
inic = nbi+0*ntD+1;
ifin = nbi+1*ntD;
FD= x0(inic:ifin);
yfD(1,1,:)=FD;
yfD(1,2,: ) = dd0rD*FD;
yfD(1,3,: ) = dd0zD*FD;
yfD(1,4,: ) = dd0rrD*FD;
yfD(1,5,: ) = dd0zzD*FD;
yfD(1,6,: ) = dd0rzD*FD;
yfD(1,7,: ) = xt(inic:ifin);

inic = nbi+1*ntD+1;
ifin = nbi+2*ntD;
GD= x0(inic:ifin);
yfD(2,1,:)=GD;
yfD(2,2,: ) = dd0rD*GD;
yfD(2,3,: ) = dd0zD*GD;
yfD(2,4,: ) = dd0rrD*GD;
yfD(2,5,: ) = dd0zzD*GD;
yfD(2,6,: ) = dd0rzD*GD;
yfD(2,7,: ) = xt(inic:ifin);

nbi = nbi + ntD*NVD;

r1D=repmat(r0D', [1 nzD]);
z1D=repmat(z0D, [nrD 1]);

%% Block 4 (DD)+++++++++++++++++++++++++++++++++++++

% %% Dreating functions and Jacobians
eval(['FDD  ' '= zeros(NVD,ntD);']);
eval(['DFDD  ' '= zeros(NVD,NVD*NDD,ntD);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,

for l=1:ntD
    pa(4)=Fcoil(l);
    pa(5)=Gcoil(l);
    
    [is,js]=ind2sub([nrD,nzD],l);
    
    if( js<nz1)
        pa(6)=0;
    else
        pa(6)=alphaf;
    end
    
    xa=reshape(yfD(:,:,l)',NVD*NDD,1); %xa=[y1(:,i,j);y2(:,i,j);y3(:,i,j);y4(:,i,j);y5(:,i,j)];
    
    %BULK
    [FDD(:,l),DFDD(:,:,l)]=equationFDDa_mg(z1D(l),r1D(l),xa,pa);
    
    if( liter==1)
    else
        [FDD(:,l),DFDD(:,:,l)]=equationFDD_mg(z1D(l),r1D(l),xa,pa);
    end
    
    if(ndD(l)==33 || ndD(l)==36)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDCoil_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==34)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDCoil1_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==35)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDRightCoil_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==21)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea21_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==8)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea8_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==7)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea7_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==19)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea19_mg(z1D(l),r1D(l),xa,pa);%
        
    end
    
    if(ndD(l)==20)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea20_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==10)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea10_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==451)
        [FDD(:,l),DFDD(:,:,l)]=equationFVertexDEA_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==412)
        [FDD(:,l),DFDD(:,:,l)]=equationFVertexDAB_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==423)
        [FDD(:,l),DFDD(:,:,l)]=equationFVertexDBC_mg(z1D(l),r1D(l),xa,pa);%
    end
    
    if(ndD(l)==13)
        [FDD(:,l),DFDD(:,:,l)]=equationFDDlinea13_mg(z1D(l),r1D(l),xa,pa);%
    end    
end

%DD
for j=1:NVD %equations
    C=FDD(j,:);
    eval(['xd' num2str(j) '= -C;']) %radial
    
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
        eval(['dd' num2str(j) num2str(k) '=sparse(B);'])
    end
end

%BLOCK D;
xd=[];
dd=[];
for i=1:NVD
    ax=[];
    for j=1:NVD
        aa0=eval(cat(2,'dd',num2str(i),num2str(j)));
        ax=[ax,aa0];
    end
    bb0=eval(cat(2,'xd',num2str(i)));
    dd=[dd;ax];
    xd=[xd;bb0'];
end

a=[dd];
b=([xd]);