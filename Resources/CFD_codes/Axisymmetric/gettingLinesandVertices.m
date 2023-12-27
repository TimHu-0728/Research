clear GD
clear FD
%A*+++++++++++++++++++++++++++++++++++++++++
%getting lines
%pointers to axis
Linea3A=zeros(1,nrA);  %left
for  i=1:nrA
    l=sub2ind([nrA,nzA],i,1);
    Linea3A(i)=l;
end
Linea20A=[];
for i=2:nrA-1
    l=sub2ind([nrA,nzA],i,nzA);
    Linea20A=[Linea20A,l];
end

Linea9A=[];  %bottom
Linea11A=[];  %top
for i=2:nzA-1
    l=sub2ind([nrA,nzA],1,i);
    Linea9A(i-1)=l;
    l=sub2ind([nrA,nzA],nrA,i);
    Linea11A(i-1)=l;
end

%vertexAD
l=sub2ind([nrA,nzA],1,nzA);
VertexADE=l;

%verxtexABD
l=sub2ind([nrA,nzA],nrA,nzA);
VertexABD=l;

%getting id
ndA           = zeros(1,ntA);
ndA(Linea9A)  = 9;
ndA(Linea20A) = 20;
ndA(Linea3A)  = 3;
ndA(Linea11A) = 11;
ndA(VertexADE)= 145;
ndA(VertexABD)= 124;

%B*+++++++++++++++++++++++++++++++++++++++++

%getting lines
%pointers to axis
Linea4B=zeros(1,nrB);  %left

for i=1:nrB
    l=sub2ind([nrB,nzB],i,1);
    Linea4B(i)=l;
    l=sub2ind([nrB,nzB],i,nzB);
    Linea10B(i)=l;
end

Linea10B=zeros(1,nrB-2);  %left

for i=2:nrB-1
    l=sub2ind([nrB,nzB],i,nzB);
    Linea10B(i-1)=l;
end
l=sub2ind([nrB,nzB],1,nzB);
VertexBDA=l;
l=sub2ind([nrB,nzB],nrB,nzB);
VertexBCD=l;

Linea11B=zeros(1,nzB-2);  %bottom
Linea12B=zeros(1,nzB-2);  %top
for i=2:nzB-1
    l=sub2ind([nrB,nzB],1,i);
    Linea11B(i-1)=l;
    l=sub2ind([nrB,nzB],nrB,i);
    Linea12B(i-1)=l;
end
%getting id
ndB=zeros(1,ntB);
ndB(Linea10B)=10;
ndB(Linea12B)=12;
ndB(Linea4B)=4;
ndB(Linea11B)=11;
ndB(VertexBDA)=241;
ndB(VertexBCD)=234;

%C*+++++++++++++++++++++++++++++++++++++++++

Linea5C=zeros(1,nrC);  %left

for i=1:nrC
    l=sub2ind([nrC,nzC],i,1);
    Linea5C(i)=l;
    l=sub2ind([nrC,nzC],i,nzC);
end
Linea13C=[];  %left
for i=2:nrC
    l=sub2ind([nrC,nzC],i,1);
    l=sub2ind([nrC,nzC],i,nzC);
    Linea13C=[Linea13C,l];
end
l=sub2ind([nrC,nzC],1,nzC);
VertexCDB=l;


Linea12C=zeros(1,nzC-2);  %bottom
Linea6C=[];  %top
for i=2:nzC-1
    l=sub2ind([nrC,nzC],1,i);
    Linea12C(i-1)=l;
end
for i=2:nzC-1
    l=sub2ind([nrC,nzC],nrC,i);
    Linea6C=[Linea6C,l];
end

%getting id
ndC=zeros(1,ntC);
ndC(Linea5C)=5;
ndC(Linea13C)=13;
ndC(Linea12C)=12;
ndC(Linea6C)=6;
ndC(VertexCDB)=342;
Linea2E=[];

%E+++++++++++++++++++++++++++++++++++++++++
for i=1:nrE
    l=sub2ind([nrE,nzE],i,1);
    Linea2E=[Linea2E,l];
end

Linea19E=[];
for i=1:nrE-1
    l=sub2ind([nrE,nzE],i,nzE);
    Linea19E=[Linea19E,l];
end
VertexEAD=sub2ind([nrE,nzE],nrE,nzE);

Linea9E=[];
Linea1E=[];
for i=2:nzE-1
    l=sub2ind([nrE,nzE],nrE,i);
    Linea9E= [Linea9E,l];
end
for i=2:nzE-1
    l=sub2ind([nrE,nzE],1,i);
    Linea1E= [Linea1E,l];
end

%getting id
ndE=zeros(1,ntE);

ndE(Linea2E)=2;
ndE(Linea1E)=1;
ndE(Linea9E)=9;
ndE(Linea19E)=19;
ndE(VertexEAD)=514;

%     hold on
%       plot(rA(Linea20A),zA(Linea20A),'bx')
%   plot(rA(Linea3A),zA(Linea3A),'bx')
%   plot(rA(Linea11A),zA(Linea11A),'bo')
%   plot(rA(Linea9A),zA(Linea9A),'bo')
%    plot(rA(VertexABD),zA(VertexABD),'bdiamond')
%
%   plot(rA(VertexADE),zA(VertexADE),'bdiamond')
%
%
%
%   plot(rB(Linea11B),zB(Linea11B),'mx')
%   plot(rB(Linea12B),zB(Linea12B),'mx')
%   plot(rB(Linea4B),zB(Linea4B),'mo')
%  plot(rB(Linea10B),zB(Linea10B),'mo')
%   plot(rB(VertexBDA),zB(VertexBDA),'msquare')
%    plot(rB(VertexBCD),zB(VertexBCD),'msquare')
%
%
% %
%   plot(rC(Linea5C),zC(Linea5C),'co')
%   plot(rC(Linea13C),zC(Linea13C),'co')
%   plot(rC(Linea12C),zC(Linea12C),'cx')
%   plot(rC(Linea6C),zC(Linea6C),'cx')
%   plot(rC(VertexCDB),zC(VertexCDB),'cdiamond')
%
%
%
%
%  plot(rE(Linea2E),zE(Linea2E),'cx')
%  plot(rE(Linea9E),zE(Linea9E),'co')
%  plot(rE(Linea1E),zE(Linea1E),'co')
%  plot(rE(Linea19E),zE(Linea19E),'cx')
% plot(rE(VertexEAD),zE(VertexEAD),'csquare')
zAa=reshape(zA,nrA*nzA,1);
zEa=reshape(zE,nrE*nzE,1);
zBa=reshape(zB,nrB*nzB,1);
zCa=reshape(zC,nrC*nzC,1);

%creating line in z for block D [zout;
GD0=[zEa(Linea19E);zEa(VertexEAD);zAa(Linea20A);zAa(VertexABD);...
    zBa(Linea10B);zBa(VertexBCD);zCa(Linea13C)];
nrD=length(GD0);

if(lcompute==1)
    
    % Block D+ (Coil)
    
    %locating closer point to lower coil
    [a,imin]=min(abs(GD0-ZminD));
    %locating closer point to lower coil
    [a,imax]=min(abs(GD0-ZmaxD));
    
    %     hold on
    %     plot(GD0,GD0,'o')
    %     plot(GD0(imin),GD0(imin),'x')
    %     plot(ZminD,ZminD,'square')
    %     plot(GD0(imax),GD0(imax),'x')
    %     plot(ZmaxD,ZmaxD,'square')
    
    % Vertical nodes of coil
    GD1=GD0;
    GD1(imax)=ZmaxD;
    GD1(imin)=ZminD;
    GD1(1:imin)   = fliplr(zout - (zout-ZminD)*tanh(2*(1-linspace(0,1,imin)))/tanh(2)); 
    GD1((imin+1):(imax-1)) = linspace(ZminD + (ZmaxD-ZminD)/(imax-imin-1),...
        ZmaxD - (ZmaxD-ZminD)/(imax-imin-1),imax-imin-1);
    GD1(imax:nrD) = H2 - (H2-ZmaxD)*tanh(1.5*(1-linspace(0,1,nrD-imax+1)))/tanh(1.5); 
    
    
    %in nzD (radial) %sólo hasta la bobina
    
    nzD=nz1;
    
    %inital mesh in domain D
    [z0D,dz0D,dzz0D]=finitasdoblegood(nzD,1);    %%%%%%%% 12/03/2021
    [r0D,dr0D,drr0D]=finitas4ordentotal(nrD,1);
    
    ntD=nrD*nzD;
    
    %domain D
    FD=zeros(nrD*nzD,1);
    GD=zeros(nrD*nzD,1);
    
    Linea19D=[]
    for is=1:nrE-1
        i=is;
        l=sub2ind([nrD,nzD],i,1);
        Linea19D=[Linea19D,l];
    end
    FD(Linea19D)=rE(1:nrE-1,nzE);
    GD(Linea19D)=zE(1:nrE-1,nzE);
    
    %vertice DEA
    VertexDEA=sub2ind([nrD,nzD],nrE,1);
    FD(VertexDEA)=rE(nrE,nzE);
    GD(VertexDEA)=zE(nrE,nzE);
    Linea20D=[]
    for is=2:nrA-1
        i=nrE+is-1;
        l=sub2ind([nrD,nzD],i,1);
        Linea20D=[Linea20D,l];
    end
    FD( Linea20D)=rA(2:nrA-1,nzA);
    GD( Linea20D)=zA(2:nrA-1,nzA);
    VertexDAB=sub2ind([nrD,nzD],nrE+nrA-1,1);
    FD(VertexDAB)=rA(nrA,nzA);
    GD(VertexDAB)=zA(nrA,nzA);
    Linea10D=[]
    for is=2:nrB-1
        i=nrE+nrA+is-2;
        l=sub2ind([nrD,nzD],i,1);
        Linea10D=[Linea10D,l];
    end
    FD( Linea10D)=rB(2:nrB-1,nzB);
    GD( Linea10D)=zB(2:nrB-1,nzB);
    VertexDBC=sub2ind([nrD,nzD],nrE+nrA+nrB-2,1);
    FD(VertexDBC)=rB(nrB,nzB);
    GD(VertexDBC)=zB(nrB,nzB);
    Linea13D=[]
    for is=2:nrC
        i=nrE+nrA+nrB+is-3;
        l=sub2ind([nrD,nzD],i,1);
        Linea13D=[Linea13D,l];
    end
    FD( Linea13D)=rC(2:nrC,nzC);
    GD( Linea13D)=zC(2:nrC,nzC);
    
    hold on
    
    
    Linea7D=[];
    Linea21D=[];
    for j=2:nzD-1 %bottom wall
        l=sub2ind([nrD,nzD],nrD,j);
        Linea7D=[Linea7D,l];
        l=sub2ind([nrD,nzD],1,j);
        Linea21D=[Linea21D,l];
    end
    FD(Linea7D)=rC(nrC,nzC)+(RminD-rC(nrC,nzC))*z0D(2:nzD-1);
    FD(Linea21D)=rE(nrE,nzE)+(RminD-rE(nrE,nzE))*z0D(2:nzD-1);
    GD(Linea7D)=zC(nrC,nzC);
    GD(Linea21D)=zE(1,nzC);
    Linea8D=[];
    for i=1:nrD %top wall
        l=sub2ind([nrD,nzD],i,nzD);
        Linea8D=[Linea8D,l];
    end
    % GD(Linea8D)=GD(:,1);
    FD(Linea8D)=RminD;
    GD(Linea8D)=GD1;
    %getting id
    
    FD=reshape(FD,nrD,nzD);
    GD=reshape(GD,nrD,nzD);
    
    for i=1:nrD
        for j=1:nzD
            GD(i,j)=GD0(i)+(GD1(i)-GD0(i))*z0D(j);
        end
    end
    FD(nrD,1)=R;
    for j=2:nzD
        FD(:,j)= FD(1,j);
    end
    %Aumentamos el el nzD.
    %discretizacion bobina
    
    
    LineaCoi1a=Linea8D;
    
    a = linspace(0,1,nzbobina);
    a = RminD + (RmaxD-RminD)*(a-0.1*sin(2*pi*a)); % Concentration of points (instead of a = linspace(RminD,RmaxD,nzbobina);)
    for j=2:nzbobina
        FD(:,nzD+j-1)= a(j) ;
        GD(:,nzD+j-1)= GD(:,nzD) ;
    end
    nzD=length(FD(1,:));
    
    a = linspace(0,1,nz2);
    a = Rout - (Rout-RmaxD)*tanh(2*(1-a))/tanh(2); % Concentration of points (instead of a = linspace(RmaxD,Rout,nz2);)
    for j=2:nz2
        FD(:,nzD+j-1)= a(j) ;
        GD(:,nzD+j-1)= GD(:,nzD) ;
    end
    nzD=length(FD(1,:));
    
    
    %  plot(rE(Linea19E),zE(Linea19E),'x')
    % hold on
    %
    % plot(FD(Linea19D), GD(Linea19D),'o')
    %
    % stop
    
    [z0D,dz0D,dzz0D]=finitas4ordentotal(nzD,1);    %%%%%%%% 12/03/2021
    %[z0D,dz0D,dzz0D]=finitas2(nzD,1);
    %[z0D,dz0D,dzz0D]=Chevitanh2th(nzD,0,1,2);
    
    
    
    
    %     hold on
    %     %      %ploting mesh E
    %     for j=1:nzD
    %         plot(FD(:,j),GD(:,j),'b-')
    %     end
    %     for i=1:nrD
    %         plot(FD(i,:),GD(i,:),'b-')
    %     end
    
    %     plot(RminD,ZminD,'x')
    %     plot(RminD,ZmaxD,'x')
    %     plot(RmaxD,ZmaxD,'x')
    %     plot(RmaxD,ZminD,'x')
    
    
    
    %      %ploting mesh E
    
    
    
    
    
    %Coil1 dentro y fuera de la bobina
    CoilD1=[]; %external (in z)
    for j=nz1:nz1+nzbobina-1
        for i=2:imin-1
            CoilD1=[CoilD1,sub2ind([nrD,nzD],i,j)];
        end
        for i=imax+1:(nrD-1)
            CoilD1=[CoilD1,sub2ind([nrD,nzD],i,j)];
        end
    end
    
    CoilD_contour=[]; % contour
    for j=nz1:nz1+nzbobina-1
        for i=imin:imax
            if i==imin || i==imax || j==nz1 || j==nz1+nzbobina-1
                CoilD_contour=[CoilD_contour,sub2ind([nrD,nzD],i,j)];
            end
        end
    end
    
    CoilD=[]; % internal (in z)
    for j=nz1+1:nz1+nzbobina-2
        for i=imin+1:imax-1
            CoilD=[CoilD,sub2ind([nrD,nzD],i,j)];
        end
    end
    
    % Region at the right of the coil
    RightCoil=[];
    for j=(nz1+nzbobina):(nzD-1)
        for i=2:(nrD-1)
            RightCoil=[RightCoil,sub2ind([nrD,nzD],i,j)];
        end
    end
    
   
    for j=nz1:nzD-1 %nz1+nzbobina:nzD-1
        Linea7D=[Linea7D,sub2ind([nrD,nzD],nrD,j)];
        Linea21D=[Linea21D,sub2ind([nrD,nzD],1,j)];
    end
    
    Linea8D=[];
    for i=1:nrD %top wall
        l=sub2ind([nrD,nzD],i,nzD);
        Linea8D=[Linea8D,l];
    end
    
    % plot(FD(CoilD1),GD(CoilD1),'gdiamond')
    % plot(FD(CoilD),GD(CoilD),'rdiamond')
    
    ndD=zeros(nrD,nzD);
    ndD(Linea7D)=7;
    ndD(Linea8D)=8;
    ndD(Linea21D)=21;
    
    
    ndD(Linea10D)=10;
    ndD(Linea13D)=13;
    
    ndD(Linea20D)=20;
    ndD(Linea19D)=19;
    ndD(CoilD)     = 33;
    ndD(CoilD1)    = 34;
    ndD(RightCoil) = 35;
    ndD(CoilD_contour) = 36;
    
    ndD(VertexDBC)=423;
    ndD(VertexDEA)=451;
    ndD(VertexDAB)=412;
    
    %ndD(nrD,1)=7;
    
    % Save preliminary mesh
    save('meshgenerator/initmesh','FD','GD', 'r0D','z0D', 'dr0D','drr0D',...
        'dz0D','dzz0D', 'nrD','nzD', 'ndD', 'CoilD', 'CoilD1' ,'Linea10D',...
        'Linea7D','Linea8D','Linea20D','Linea13D','Linea21D','Linea19D',...
        'VertexDBC','VertexDEA','VertexDAB','zout','H2','R','Rout','imax',...
        'imin','nz1','nz2' ,'nzbobina')
    
    % Compute final mesh
    MainmeshD
    
    % Load final mesh
    load('meshgenerator/finalmesh' )
else
    load('meshgenerator/finalmesh' )
end

%    %C*+++++++++++++++++++++++++++++++++++++++++
%
%         Linea5C=zeros(1,nrC);  %left
%
%     for i=1:nrC
%         l=sub2ind([nrC,nzC],i,1);
%         Linea5C(i)=l;
%         l=sub2ind([nrC,nzC],i,nzC);
%     end
%     Linea13C=[];  %left
%     for i=2:nrC
%         l=sub2ind([nrC,nzC],i,1);
%         l=sub2ind([nrC,nzC],i,nzC);
%         Linea13C=[Linea13C,l];
%     end
%        l=sub2ind([nrC,nzC],1,nzC);
%     VertexCDB=l;
%
%
%     Linea12C=zeros(1,nzC-2);  %bottom
%     Linea6C=[];  %top
%    for i=2:nzC-1
%         l=sub2ind([nrC,nzC],1,i);
%         Linea12C(i-1)=l;
%         end
%      for i=2:nzC-1
%         l=sub2ind([nrC,nzC],nrC,i);
%         Linea6C=[Linea6C,l];
%     end
%
%      %getting id
%     ndC=zeros(1,ntC);
%     ndC(Linea5C)=5;
%     ndC(Linea13C)=13;
%     ndC(Linea12C)=12;
%     ndC(Linea6C)=6;
%     ndC(VertexCDB)=342;
%



