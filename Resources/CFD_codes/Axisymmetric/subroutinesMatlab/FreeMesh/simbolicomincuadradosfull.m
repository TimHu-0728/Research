clear all    
%number of points 
    M=16;

%block 1++++++++
for k=1:M  %variable
for i=1:M
xa(k)=sym( ['xa',num2str(k)],'real');
ya(k)=sym( ['ya',num2str(k)],'real');
va(k)=sym( ['va',num2str(k)],'real');
end
end


xa=xa';
ya=ya';
va=va';
% tic
% c=a\xa
% toc
       
    % coordenates
    x0=xa(1) ; y0=ya(1); f0=va(1); 
    x1r=xa(2:M)-x0;
    y1r=ya(2:M)-y0;
    f1r=va(2:M)-f0;
    % distances (al cuadrado) to the center
    % dmin=minimo de esas distancias
    
   
  

    x1r=x1r ; y1r=y1r ;
    
    N=5;%sencond order aproximation
    N0=M-1
    X=zeros(N0,N)*x0;
   x1r=zeros(1,N0)*x0;
   y1r=zeros(1,N0)*x0;
   x2r=zeros(1,N0)*x0;
   y2r=zeros(1,N0)*x0;
   x1ry1r=zeros(1,N0)*x0;
   x1ry2r=zeros(1,N0)*x0;
for k=1:N0
x1r(k)=sym( ['x1r',num2str(k)],'real');
y1r(k)=sym( ['y1r',num2str(k)],'real');
x2r(k)=sym( ['x2r',num2str(k)],'real');
y2r(k)=sym( ['y2r',num2str(k)],'real');
x1ry1r(k)=sym( ['x1ry1r',num2str(k)],'real');
end
    
 X(:,1)=x1r';
 X(:,2)=y1r';
 X(:,3)=x2r';
 X(:,4)=y2r';
 X(:,5)=x1ry1r';
 Y=f1r;
 
 %creating matrix A
 
 for k=1:N
 for i=1:N
 A(k,i)=sym( ['A',num2str(k),num2str(i)],'real');
 end
 end
 
%A=inv(X'*X);
 
  b=X'*Y; 
 
 lambda=A*b;
 
 dux=lambda(1);
 duy=lambda(2);
 duxx=2.0*lambda(3);
 duyy=2.0*lambda(4);
 duxy=lambda(5);
 va=va';
 
 dux=simplify(dux);
 duy=simplify(duy);
 duxx=simplify(duxx);
 duyy=simplify(duyy);
 duxy=simplify(duxy);
 
 %getting colocation matrix
 for k=1:M
 duxva(:,k)=dux;
 duxxva(:,k)=duxx;
 duyva(:,k)=duy;
 duyyva(:,k)=duyy;
 duxyva(:,k)=duxy;
 for l=1:M
   if(k==l)  
 duxva(:,k)=subs(duxva(:,k),va(l),1);
 duxxva(:,k)=subs(duxxva(:,k),va(l),1);
 duyva(:,k)=subs(duyva(:,k),va(l),1);
 duyyva(:,k)=subs(duyyva(:,k),va(l),1);
 duxyva(:,k)=subs(duxyva(:,k),va(l),1);
   else
 duxva(:,k)=subs(duxva(:,k),va(l),0);
 duyva(:,k)=subs(duyva(:,k),va(l),0);
 duxxva(:,k)=subs(duxxva(:,k),va(l),0);
 duyyva(:,k)=subs(duyyva(:,k),va(l),0);
 duxyva(:,k)=subs(duxyva(:,k),va(l),0);
   end
 end
 end
%  dux0=0;
% for k=1:6
% 
%  dux0=dux0+duyva(:,k)*va(k);
% end
% 
% 
% 
%  c=simplify(duy-dux0)
 %explicity
 A=reshape(A,1,5*5);
%  for k=1:5
%  for i=1:5
%  dux=subs(dux,A(k,i),1);
%  end



 tic
 matlabFunction(dux,duy,duxx,duyy,duxy,'file','derivativesN16.m','vars',{va,x1r,y1r,x2r,y2r,x1ry1r,A});
 toc
 
 tic
 %collocation matrix
% matlabFunction(duxva,duyva,duxxva,duyyva,'file','eigen/jacobians3D/collocationN7a.m','vars',{x1r,y1r,x2r,y2r,x1ry1r,A},'Optimize',false);
 matlabFunction(duxva,duyva,duxxva,duyyva,duxyva,'file','collocationN16.m','vars',{x1r,y1r,x2r,y2r,x1ry1r,A});

 toc
 

 
 
 
 
 