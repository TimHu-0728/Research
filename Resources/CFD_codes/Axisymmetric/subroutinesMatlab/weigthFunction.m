function [w]=weigthFunction(d)
if (d<1)
   w=2/3-9/3*d.^2+19/192*d.^3-5/512*d.^5;
else
w=0.*d;
end