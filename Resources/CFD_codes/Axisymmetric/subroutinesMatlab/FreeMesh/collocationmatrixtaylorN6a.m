function [dux,duy,duxx,duyy]=collocationmatrixtaylorN6a(xa,ya)


   
    
    %number of points 
    M=length(xa);

    
    % coordenates
    x0=xa(1) ; y0=ya(1); 
    x1r=xa(2:M)-x0;
    y1r=ya(2:M)-y0;
   
    % distances (al cuadrado) to the center
    % dmin=minimo de esas distancias
    
   
  

    x2r=x1r.^2 ; y2r=y1r.^2;
    x1ry1r=y1r.*x1r;
    
    N=5;%sencond order aproximation
    X=zeros(M-1,N);
    C= con
 X(:,1)=x1r;
 X(:,2)=y1r;
 X(:,3)=x2r;
 X(:,4)=y2r;
 X(:,5)=x1ry1r;
 
 
 A=inv(X);
 
  A=reshape(A,1,5*5);  

 [dux,duy,duxx,duyy]=collocationtaylorN6(x1r',y1r',x2r',y2r',x1ry1r',A);

   
    
  


    
    
        


