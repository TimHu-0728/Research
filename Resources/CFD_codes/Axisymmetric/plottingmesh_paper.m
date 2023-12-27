% Create new figure
figure

%real points
rA=reshape(rA,nrA*nzA,1);
zA=reshape(zA,nrA*nzA,1);
%real points
rB=reshape(rB,nrB*nzB,1);
zB=reshape(zB,nrB*nzB,1);
%real points
rC=reshape(rC,nrC*nzC,1);
zC=reshape(zC,nrC*nzC,1);
%real points
rE=reshape(rE,nrE*nzE,1);
zE=reshape(zE,nrE*nzE,1);
%mapping points
r1E=reshape(r1E,nrE*nzE,1);
z1E=reshape(z1E,nrE*nzE,1);

rD=FD;
zD=GD;

%% PLOT

% Parameters
linewidth     = 2;
linewidthmesh = 0.1;
hold on
colormap = [228,26,28
255,127,0
77,175,74
55,126,184
152,78,163];

% Mesh
%real points
rA=reshape(rA,nrA,nzA);
zA=reshape(zA,nrA,nzA);
%real points
rB=reshape(rB,nrB,nzB);
zB=reshape(zB,nrB,nzB);
%real points
rC=reshape(rC,nrC,nzC);
zC=reshape(zC,nrC,nzC);
%real points
rE=reshape(rE,nrE,nzE);
zE=reshape(zE,nrE,nzE);
%mapping points
r1E=reshape(r1E,nrE,nzE);
z1E=reshape(z1E,nrE,nzE);

%  %ploting mesh A
for j=1:nzA
    plot(rA(:,j),zA(:,j),'-','linewidth',linewidthmesh,'color',colormap(1,:)/255)
end
for i=1:nrA
    plot(rA(i,:),zA(i,:),'-','linewidth',linewidthmesh,'color',colormap(1,:)/255)
end

%     %ploting mesh B
for j=1:nzB
    plot(rB(:,j),zB(:,j),'-','linewidth',linewidthmesh,'color',colormap(2,:)/255)
end
for i=1:nrB  %manual
    plot(rB(i,:),zB(i,:),'-','linewidth',linewidthmesh,'color',colormap(2,:)/255)
end

%      %ploting mesh C
for j=1:nzC
    plot(rC(:,j),zC(:,j),'-','linewidth',linewidthmesh,'color',colormap(3,:)/255)
end
for i=1:nrC
    plot(rC(i,:),zC(i,:),'-','linewidth',linewidthmesh,'color',colormap(3,:)/255)
end

%      %ploting mesh D
for j=1:nzD
    plot(FD(:,j),GD(:,j),'-','linewidth',linewidthmesh,'color',colormap(4,:)/255)
end
for i=1:nrD
    plot(FD(i,:),GD(i,:),'-','linewidth',linewidthmesh,'color',colormap(4,:)/255)
end

%      %ploting mesh E
for j=1:nzE
    plot(rE(:,j),zE(:,j),'-','linewidth',linewidthmesh,'color',colormap(5,:)/255)
end
for i=1:nrE
    plot(rE(i,:),zE(i,:),'-','linewidth',linewidthmesh,'color',colormap(5,:)/255)
end

% Contours
plot(rA(Linea20A),zA(Linea20A),'k-','Linewidth',linewidth)
plot(rA(Linea3A),zA(Linea3A),'k-','Linewidth',linewidth)
plot(rA(Linea11A),zA(Linea11A),'k-','Linewidth',linewidth)
plot(rA(Linea9A),zA(Linea9A),'k-','Linewidth',linewidth)
plot(rB(Linea12B),zB(Linea12B),'k-','Linewidth',linewidth)
plot(rB(Linea4B),zB(Linea4B),'k-','Linewidth',linewidth)
plot(rC(Linea5C),zC(Linea5C),'k-','Linewidth',linewidth)
plot(rD([Linea10D,Linea13D]),zD([Linea10D,Linea13D]),'k-','Linewidth',linewidth)
plot(rC(Linea12C),zC(Linea12C),'k-','Linewidth',linewidth)
plot(rC(Linea6C),zC(Linea6C),'k-','Linewidth',linewidth)
plot(rE(Linea2E),zE(Linea2E),'k-','Linewidth',linewidth)
plot(rE(Linea1E),zE(Linea1E),'k-','Linewidth',linewidth)
plot(rE(Linea19E),zE(Linea19E),'k-','Linewidth',linewidth)
plot(rD([Linea7D,fliplr(Linea8D),fliplr(Linea21D)]),zD([Linea7D,fliplr(Linea8D),fliplr(Linea21D)]),'k-','Linewidth',linewidth)
plot(rB(VertexBCD),zB(VertexBCD),'k.','Linewidth',linewidth)
plot(rA(VertexABD),zA(VertexABD),'k.','Linewidth',linewidth)
plot(rA(VertexADE),zA(VertexADE),'k.','Linewidth',linewidth)
plot([RminD,RmaxD,RmaxD,RminD,RminD],[ZminD,ZminD,ZmaxD,ZmaxD,ZminD],'k--','linewidth',linewidth/2)

% Format
axis equal
xlabel('r [m]')
ylabel('z [m]')
xlim([0,0.15])
ylim([-0.1,0.15])