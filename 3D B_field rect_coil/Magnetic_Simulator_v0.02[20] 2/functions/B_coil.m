function B=B_coil(x,y,z,a,b,mu,I,N,h)
% Function calculate the magnetic field around a rectangular coil with input parameters:
% Position r = [x,y,z] (m),
% Coil length a (m), width b (m), total height h (m), number of stack N, current I (A), mu (N*A^(-2))
% 
% B=B_coil(x,y,z,a,b,mu,I,N,h)
%
% INPUTS:
%
% OUTPUTS:
%

elvt=h/N;    
B=zeros(3,1);
    for i=1:N
        B=B+B_rect(x,y,z+(-1)^i*floor(i/2)*elvt,a,b,mu,I);
    end
end