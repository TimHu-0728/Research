%correcting volumen
V=0;
for j=2:nzA
    ap = (gA(nrA,j)-gaxisA(nrA,j))*fA(nrA,j);
    am = (gA(nrA,j-1)-gaxisA(nrA,j-1))*fA(nrA,j-1);
    V  = V+0.5*(fA(nrA,j)-fA(nrA,j-1))*(ap+am);
end
V = V*2*pi;

% Nominal height
hnominal=V/(pi*R^2);

% gA=gA*H/ hnominal;

 

 
