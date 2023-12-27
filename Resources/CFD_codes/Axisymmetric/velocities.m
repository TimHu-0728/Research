
function a=velocities(za,ra,zb,rb,zc,rc,zd,rd,ze,re,pA,pB,pC,pD,pE)
a21=max(max(pA));
a22=max(max(pB));
a23=max(max(pC));
a24=max(max(pD));
a25=max(max(pE));
a11=min(min(pA));
a12=min(min(pB));
a13=min(min(pC));
a14=min(min(pD));
a15=min(min(pE));



a2=max(a21,a22);
a2=max(a2,a23);
a2=max(a2,a24);
a2=max(a2,a25);

a1=min(a11,a12);
a1=min(a1,a13);
a1=min(a1,a14);
a1=min(a1,a15);

%v1=0:2500:45000;%for paper
v1=[a1:(a2-a1)/30:a2];
v1=0:10000:2e5;
hold on
contourf(rb,zb,pB,v1)
contourf(ra,za,pA,v1)
contourf(rc,zc,pC,v1)
contourf(rd,zd,pD,v1)
contourf(re,ze,pE,v1)

nrA=length(ra(:,1));
plot(ra(nrA,:),za(nrA,:),'k-')
%plot(za(1,:),-ra(nrA,:),'r-')
axis('equal')