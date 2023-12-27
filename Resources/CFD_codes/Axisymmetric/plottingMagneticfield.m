 rA=reshape(rA,nrA*nzA,1) ;
 zA=reshape(zA,nrA*nzA,1) ;
r=[rA;rB;rC;rD;rE];
     z=[zA;zB;zC;zD;zE];
     
     
     phiA=reshape(phiA,nrA*nzA,1);
     phiB=reshape(phiB,nrB*nzB,1);
     phiC=reshape(phiC,nrC*nzC,1);
     phiE=reshape(phiE,nrE*nzE,1);

     psiA=reshape(psiA,nrA*nzA,1);
     psiB=reshape(psiB,nrB*nzB,1);
     psiC=reshape(psiC,nrC*nzC,1);
     psiE=reshape(psiE,nrE*nzE,1);
        %HzA=reshape(HzA,nrA,nzA);
 
     
     HA=(HrA.^2+HzA.^2).^0.5;
     HB=(HrB.^2+HzB.^2).^0.5;
     HC=(HrC.^2+HzC.^2).^0.5;
     HD=(HrD.^2+HzD.^2).^0.5;
     HE=(HrE.^2+HzE.^2).^0.5;
%      
%      HA=phiA;
%       HB=phiB;
%           HC=phiC;
%              HD=phiD;
%          HE=phiE;
%      
     
     H=[HA;HB;HC;HD;HE];
    %H=[psiA;psiB;psiC;psiD;psiE];
     H=[phiA;phiB;phiC;phiD;phiE];
    z1=[min(z):(max(z)-min(z))/999:max(z)];
    r1=[0:max(r)/999:max(r)];
      [R,Z] = meshgrid(r1,z1);  
      H1 = griddata(r,z,H,R,Z);
       contourf(R,Z,H1)
       hold on
       plot(fA(nrA,:),gA(nrA,:)-gaxisA(nrA,:),'r-')
       stop
%  
         tri=delaunay(r,z);
         tricontour(tri,r,z,H,50)

        scatter(r,z,[],H,'filled');
%        hold on
%             plot(zb(1,:),rb(1,:),'x')
%             hold on
%             plot(zb(1,:),rb(nrB,:),'o')
    
stop