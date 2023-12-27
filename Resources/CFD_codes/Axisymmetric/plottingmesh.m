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
hold on




plot(rA(Linea20A),zA(Linea20A),'bx')
plot(rA(Linea3A),zA(Linea3A),'bx')
plot(rA(Linea11A),zA(Linea11A),'bo')
plot(rA(Linea9A),zA(Linea9A),'bo')
plot(rA(VertexABD),zA(VertexABD),'bdiamond')

plot(rA(VertexADE),zA(VertexADE),'bdiamond')



plot(rB(Linea11B),zB(Linea11B),'mx')
plot(rB(Linea12B),zB(Linea12B),'mx')
plot(rB(Linea4B),zB(Linea4B),'mo')
plot(rB(Linea10B),zB(Linea10B),'mo')
plot(rB(VertexBDA),zB(VertexBDA),'msquare')
plot(rB(VertexBCD),zB(VertexBCD),'msquare')


%
plot(rC(Linea5C),zC(Linea5C),'co')
plot(rC(Linea13C),zC(Linea13C),'co')
plot(rC(Linea12C),zC(Linea12C),'cx')
plot(rC(Linea6C),zC(Linea6C),'cx')
plot(rC(VertexCDB),zC(VertexCDB),'cdiamond')




plot(rE(Linea2E),zE(Linea2E),'cx')
plot(rE(Linea9E),zE(Linea9E),'co')
plot(rE(Linea1E),zE(Linea1E),'co')
plot(rE(Linea19E),zE(Linea19E),'cx')
plot(rE(VertexEAD),zE(VertexEAD),'csquare')



%plotting mesh D
%  plot(rD,zD,'x')
hold on
rD=FD;
zD=GD;
plot(rD(Linea19D),zD(Linea19D),'gdiamond')

plot(rD(Linea10D),zD(Linea10D),'cdiamond')

plot(rD(Linea20D),zD(Linea20D),'rdiamond')
plot(rD(Linea13D),zD(Linea13D),'mdiamond')
plot(rD(Linea21D),zD(Linea21D),'bdiamond')
%external BC


plot(rD(Linea8D),zD(Linea8D),'square')
plot(rD(Linea7D),zD(Linea7D),'o')

% plot(rD(Linea78D),zD(Linea78D),'diamond')

plot(rD(Linea10D),zD(Linea10D),'o')

plot(rD(Linea13D),zD(Linea13D),'x')




plot(rD(VertexDBC),zD(VertexDBC),'bsquare')
plot(rD(VertexDEA),zD(VertexDEA),'bsquare')
plot(rD(VertexDAB),zD(VertexDAB),'bsquare')

lplot=1

if(lplot==1)
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
        plot(rA(:,j),zA(:,j),'g-')
    end
    for i=1:nrA
        plot(rA(i,:),zA(i,:),'g-')
    end
    
    
    
    %     %ploting mesh B
    for j=1:nzB
        plot(rB(:,j),zB(:,j),'r-')
    end
    for i=1:nrB  %manual
        plot(rB(i,:),zB(i,:),'r-')
    end
    
    %      %ploting mesh C
    for j=1:nzC
        plot(rC(:,j),zC(:,j),'c-')
    end
    for i=1:nrC
        plot(rC(i,:),zC(i,:),'c-')
    end
    %      %ploting mesh E
    for j=1:nzE
        plot(rE(:,j),zE(:,j),'m-')
    end
    for i=1:nrE
        plot(rE(i,:),zE(i,:),'m-')
    end
    
    
    %      %ploting mesh E
    for j=1:nzD
        plot(FD(:,j),GD(:,j),'b-')
    end
    for i=1:nrD
        plot(FD(i,:),GD(i,:),'b-')
    end
    
    
end
