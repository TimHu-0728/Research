
function w=filter6a(w,nz)
   
%central
     alpha=0.45;
    % alpha=0.1;
	a0=11/16+5*alpha/8;
    a1=15/32+17*alpha/16;
    a2=-3/16+3*alpha/8;
    a3=1/32-alpha/16;
    ac(1)=a0;
    ac(2)=a1;
    ac(3)=a2;
    ac(4)=a3;
%point 2
a12=1/64+31*alpha/32;
a22=29/32+3*alpha/16;
a32=15/64+17*alpha/32;
a42=-5/16+5*alpha/8;
a52=15/64-15*alpha/32;
a62=-3/32+3*alpha/16;
a72=1/64-alpha/32;
%point 3
a13=-1/64+alpha/32;
a23=3/32+13*alpha/16;
a33=49/64+15*alpha/32;
a43=5/16+3*alpha/8;
a53=-15/64+15*alpha/32;
a63=3/32-3*alpha/16;
a73=-1/64+alpha/32;

 
	nz2=nz-2;
  
     
          
     
     r=zeros(1,nz);
     a=zeros(1,nz);
     b=zeros(1,nz);
     c=zeros(1,nz);
	for j=1:nz
     
	b(j)=1.d0;
	a(j)=alpha;
	c(j)=alpha;
    r(j)=w(j);
  
    
     if( j<nz-2)&&(j>3) %centrals..... 
     r(j)=0;
      for l=0:3
        
      r(j)=r(j)+0.5*ac(l+1)*(w(j+l)+w(j-l));
      end
	  end  
    if( j==2)
    r(j)=a12*w(1)+a22*w(2)+a32*w(3)+a42*w(4)+a52*w(5)+a62*w(6)+a72*w(7);   
    end  
    if( j==3)
    r(j)=a13*w(1)+a23*w(2)+a33*w(3)+a43*w(4)+a53*w(5)+a63*w(6)+a73*w(7);   
    end 
    if (j==nz-1);
    r(j)=a12*w(nz)+a22*w(nz-1)+a32*w(nz-2)+a42*w(nz-3)+a52*w(nz-4)+a62*w(nz-5)+a72*w(nz-6);   
    end
    if( j==nz-2);
     r(j)=a13*w(nz)+a23*w(nz-1)+a33*w(nz-2)+a43*w(nz-3)+a53*w(nz-4)+a63*w(nz-5)+a73*w(nz-6); 
    end
    end

     A = spdiags([a' b' c'], -1:1, nz, nz);
     A(1,2)=0;
     A(nz,nz-1)=0;
 
    w= A\r';
    
	


    
