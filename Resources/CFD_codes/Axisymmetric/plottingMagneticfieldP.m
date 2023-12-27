     r=[rA;rB;rC;rD;rE];
     z=[zA;zB;zC;zD;zE];
     
     
     phiAP=reshape(phiAP,nrA*nzA,1);
     phiBP=reshape(phiBP,nrB*nzB,1);
     phiCP=reshape(phiCP,nrC*nzC,1);
     phiEP=reshape(phiEP,nrE*nzE,1);
        %HzA=reshape(HzA,nrA,nzA);
 
     
%      HA=(HrA.^2+HzA.^2).^0.5;
%      HB=(HrB.^2+HzB.^2).^0.5;
%      HC=(HrC.^2+HzC.^2).^0.5;
%      HD=(HrD.^2+HzD.^2).^0.5;
%      HE=(HrE.^2+HzE.^2).^0.5;
     
%      HA=phiA;
%       HB=phiB;
%           HC=phiC;
%              HD=phiD;
%          HE=phiE;
%      
     
     H=[phiAP;phiBP;phiCP;phiDP;phiEP];
     H=real(H);
    z1=[min(z):(max(z)-min(z))/999:max(z)];
    r1=[0:max(r)/999:max(r)];
      [R,Z] = meshgrid(r1,z1);  
      H1 = griddata(r,z,H,R,Z);
       contourf(R,Z,H1)
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