
rA=reshape(rA,nrA,nzA);
zA=reshape(zA,nrA,nzA);

%plotting inteface
plot(fA(nrA,:),gA(nrA,:),'r-','LineWidth',3)
hold on
contourf(rA,zA,wA); title ('axial velocity field') 
xlabel('r') 
ylabel('z')