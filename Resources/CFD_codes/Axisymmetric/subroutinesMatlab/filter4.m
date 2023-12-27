
function w=filter4(w,nz)
      alpha=0.45;
	c1=(5+6*alpha)/8.;
	c2=0.5*(1+2*alpha);
	c3=-(1-2*alpha)/8.;
 
	nz1=nz-1;
    nz2=nz-2;
    nz3=nz-3;
    nz4=nz-4;
     
          
        f1=15*w(1)/16+(4*w(2)-6*w(3)+4*w(4)-w(5))/16;
        fn=15*w(nz)/16+(4*w(nz-1)-6*w(nz-2)+4*w(nz-3)-w(nz-4))/16;
        f2=3./4.*w(2)+1./16.*(w(1)+6*w(3)-4*w(4)+w(5));
	    fn1=3./4.*w(nz1)+1./16.*(w(nz)+6*w(nz2)-4*w(nz3)+w(nz4));
	
     r=zeros(1,nz4);
     a=zeros(1,nz4);
     b=zeros(1,nz4);
     c=zeros(1,nz4);
	for j=1:nz-4
	b(j)=1.d0;
	a(j)=alpha;
	c(j)=alpha;
	r(j)=c1*w(j+2)+0.5*c2*(w(j+3)+w(j+1))+0.5*c3*(w(j+4)+w(j));
    end
% 	b=2*ones(1,nz4);
%     a=ones(1,nz4);
%     c=-ones(1,nz4);
    r(1)=r(1)-alpha*f2;
	r(nz4)=r(nz4)-alpha*fn1;
     A = spdiags([a' b' c'], -1:1, nz4, nz4);

 
    w(3:nz2)= A\r';
    w(1)=f1;
	w(2)=f2;
	w(nz)=fn;
	w(nz1)=fn1;


    
