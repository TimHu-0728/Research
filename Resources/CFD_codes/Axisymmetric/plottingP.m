nbi = 0;
for j=1:nbl-1
    bl = list_block{j};
    NVAR = eval(['NV' bl]);
    nunkn = eval(['nt' bl]); %number of gridpoints in block
    for i=1:NVAR
        inic = nbi+(i-1)*nunkn+1;
        ifin = nbi+i*nunkn;
        lv = eval(['list_var_' bl '{' num2str(i) '}']);
        lv=[lv 'P'];
        order = [lv '= Vf(' num2str(inic) ':' num2str(ifin) ',' num2str(lt) ',' num2str(is) ');'];
        evalc(order);
    end
    nbi = nbi + nunkn*NVAR; 
end

%readingbasic
xotoredeable

%energy of the perturbation
Eda=(0.5*(abs(uAP).^2+abs(wAP).^2));


Eda=reshape(Eda,nrA,nzA);

fAP=reshape(fAP,nrA,nzA)
gAP=reshape(gAP,nrA,nzA)
wAP=reshape(wAP,nrA,nzA);
uAP=reshape(uAP,nrA,nzA);
pAP=reshape(pAP,nrA,nzA);
phiAP=reshape(phiAP,nrA,nzA);
phiBP=reshape(phiBP,nrB,nzB);

rA1=reshape(rA,nrA,nzA)
zA1=reshape(zA,nrA,nzA)
 hold on
 %contourf(rA,zA,real(uAP),'-')
 hold on
% 
% contourf(zc,rc,real(Edc),'-')
% hold on



  plot(rA(nrA,:),fAP(nrA,:),'r-','LineWidth',2)
stop  
axis('equal')

%%%%%
figure
plot(zA(nrA,:),FAP(nrA,:),'r-','LineWidth',2)
figure
plot(zA(nrA,:),FAP(nrA,:)/max(FAP(nrA,:)),'r-','LineWidth',2)

[aux1, aux2]=max(abs(FAP(nrA,:)))
zA(nrA, aux2)
FA(nrA, aux2)
%hold on
%plot(zA(nrA,:),FA(nrA,:),'b-','LineWidth',2)