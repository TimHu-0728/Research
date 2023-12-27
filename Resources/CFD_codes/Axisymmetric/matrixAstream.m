

    

bp=0;

xt = bp*x0; 






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




%%%%%%%%%%%%%%BLOCK A%%%%%%%%%%%%%%A+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% %% Creating functions and Jacobians

eval(['FAA  ' '= zeros(1,ntA);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,
eval(['DFAA  ' '= zeros(1,NVA*NDA,ntA);']); %1 w-Momentum,2 u-Momentum,3-continuity, 4...7 tauss,
%evaluating numerical Jacobians


for l=1:ntA
    
    xa=reshape(yfA(:,:,l)',NVA*NDA,1); %xa=[y1(:,i,j);y2(:,i,j);y3(:,i,j);y4(:,i,j);y5(:,i,j)];
   [FAA(:,l),DFAA(:,:,l)]=equationFAAstream(z1A(l),r1A(l),xa,pa);
    
    if(ndA(l)==9)  %
      [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea3stream(z1A(l),r1A(l),xa,pa);%   
        
    end
    if(ndA(l)==11)
        if(l==Linea11A(length(Linea11A)))
            [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea3stream(z1A(l),r1A(l),xa,pa);%
        else
            [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea3stream(z1A(l),r1A(l),xa,pa);%
        end
       
    end
    if (ndA(l)==3)
        [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea3stream(z1A(l),r1A(l),xa,pa);%
        %[MA(l),MAz(l),HrA(l),HzA(l)]=MagneticfieldAAa(z1A(l),r1A(l),xa,pa);
    end
    if(ndA(l)==20)  %
        [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea3stream(z1A(l),r1A(l),xa,pa);%
        %[FDA(:,l),DFDA(:,:,l)]=equationFDAlinea20(z1A(l),r1A(l),xa,pa);%
    end
    if(ndA(l)==145)
        [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea3stream(z1A(l),r1A(l),xa,pa);%
        %[FDA(:,l),DFDA(:,:,l)]=equationFVertexDEAA(z1A(l),r1A(l),xa,pa);%
        
      end
    
      if(ndA(l)==124)
        [FAA(:,l),DFAA(:,:,l)]=equationFAAlinea3stream(z1A(l),r1A(l),xa,pa);%
       
    end
    
    
end


%Mouting AA
for j=1:1 %equations
    C=FAA(j,:);
    eval(['xa' num2str(j) '= -C;']) %radial
    
    
    for k=1:NVA  %jacobians
        km=(k-1)*NDA+1;
        kp=k*NDA;
        D=squeeze(DFAA(j,km:kp,:));
        B=spdiags(D(1,:)',0,ntA, ntA) ...
            + spdiags(D(2,:)',0, ntA, ntA)* dd0rA...
            + spdiags(D(3,:)',0, ntA, ntA)* dd0zA...
            + spdiags(D(4,:)',0, ntA, ntA)* dd0rrA...
            + spdiags(D(5,:)',0, ntA, ntA)* dd0zzA...
            + spdiags(D(6,:)',0, ntA, ntA)* dd0rzA...
            + spdiags(D(7,:)',0, ntA, ntA)* bp;
        eval(['aa' num2str(j) num2str(k) '=sparse(B);'])
    end
end


a1=[aa13];
b1=[xa1'];












% 
% %BLOCK A
xa=[];
aa=[];
for i=1:1
    ax=[];
    for j=3:3
        aa0=eval(cat(2,'aa',num2str(i),num2str(j)));
        ax=[ax,aa0];
    end
    bb0=eval(cat(2,'xa',num2str(i)));
    aa=[aa;ax];
    xa=[xa;bb0'];
end










%MATRX A and vector B
a=[aa];
b=([xa]);

c=a-a1

%
% a=[aa,ab,ac;ba,bb,bc;ca,cb,cc];
% b=([xa;xb;xc]);

% a=[aa,ab,ac;ba,bb];
% b=([xa;xb]);
%  a=[aa,ab,ac,ad;ba,bb,bc,bd;ca,cb,cc,cd;da,db,dc,dd];
%  b=([xa;xb;xc;xd]);

%  a=[aa];
%  b=[xa];

