
%A*+++++++++++++++++++++++++++++++++++++++++
%getting lines
    %pointers to axis
    Linea3A=zeros(1,nrA);  %left
    Linea9Ar=zeros(1,nrA-2);  %left
 for  i=1:nrA   
        l=sub2ind([nrA,nzA],i,1);
        Linea3A(i)=l;
end
    
     for i=2:nrA-1
     l=sub2ind([nrA,nzA],i,nzA);
        Linea9Ar(i-1)=l;
     end 
   
    
    Linea9Ab=zeros(1,nzA-2);  %bottom
    Linea11A=zeros(1,nzA-2);  %top
    for i=2:nzA-1
        l=sub2ind([nrA,nzA],1,i);
        Linea9Ab(i-1)=l;
        l=sub2ind([nrA,nzA],nrA,i);
        Linea11A(i-1)=l;
    end
    
    %vertexAD
    l=sub2ind([nrA,nzA],1,nzA);
    VertexAD=l; 
    %verxtexABD
    l=sub2ind([nrA,nzA],nrA,nzA);
    VertexABD=l; 
    
    %getting id
    ndA=zeros(1,ntA);
    ndA(Linea9Ab)=98;
    ndA(Linea9Ar)=99;
    ndA(Linea3A)=3;
    ndA(Linea11A)=11;
    ndA(VertexAD)=14;
    ndA(VertexABD)=124;

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
    Linea13C=zeros(1,nrC-1);  %left
    for i=2:nrC
        l=sub2ind([nrC,nzC],i,1);
        l=sub2ind([nrC,nzC],i,nzC);
        Linea13C(i-1)=l;
    end
       l=sub2ind([nrC,nzC],1,nzC);
    VertexCDB=l;
    
    
    Linea12C=zeros(1,nzC-2);  %bottom
    Linea6C=zeros(1,nzC-2);  %top
    for i=2:nzC-1
        l=sub2ind([nrC,nzC],1,i);
        Linea12C(i-1)=l;
        l=sub2ind([nrC,nzC],nrC,i);
        Linea6C(i-1)=l;
    end
    
     %getting id
    ndC=zeros(1,ntC);
    ndC(Linea5C)=5;
    ndC(Linea13C)=13;
    ndC(Linea12C)=12;
    ndC(Linea6C)=6;
    ndC(VertexCDB)=342;
    
    %E+++++++++++++++++++++++++++++++++++++++++   
       for i=2:nzE-1
        l=sub2ind([nrE,nzE],1,i);
        Linea15E(i-1)=l;
        l=sub2ind([nrE,nzE],nrE,i);
        Linea17E(i-1)=l;
    end
    
    for i=1:nrE
        l=sub2ind([nrE,nzE],i,1);
        Linea14E(i)=l;
        l=sub2ind([nrE,nzE],i,nzE);
        Linea16E(i)=l;
    end 
    
     %getting id
    ndE=zeros(1,ntE);
    ndE(Linea14E)=14;
    ndE(Linea15E)=15;
    ndE(Linea16E)=16;
    ndE(Linea17E)=17;
    
    
    
    
    
    
    
    
    
    
    
    