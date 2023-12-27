

%This code generate  collocation matrices in a domain (r,z)

%ddr the matrix  containing the first derivate in r 
%ddrr the matrix  containing the second derivate in r 
%ddz the matrix  containing the first derivate in z
%ddzz the matrix  containing the second derivate in z 

%We needd 2 vector+ 1 matrix 
%r vector containing the radial position 1:N
%z vector containing the radial position 1:N
%la matrix MxN of pointers containing the position of the M closer points to (r,z).It could be: 7,8 or 9  
clear all
M=7
%select the approximation 
switch M
    case 7
[ddr,ddz,ddr2,ddz2]=collocationmatrixN7(la,r,z);
   case 8
[ddr,ddz,ddr2,ddz2]=collocationmatrixN8(la,r,z);
    case 9
[ddr,ddz,ddr2,ddz2]=collocationmatrixN9(la,r,z);
end