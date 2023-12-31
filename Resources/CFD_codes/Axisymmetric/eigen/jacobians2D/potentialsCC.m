function [dphir,dphiz,dphirr,dphizz,dphirz,dpsir,dpsiz,dpsirr,dpsizz,dpsirz] = potentialsCC(z1,r1,in3,in4)
%POTENTIALSCC
%    [DPHIR,DPHIZ,DPHIRR,DPHIZZ,DPHIRZ,DPSIR,DPSIZ,DPSIRR,DPSIZZ,DPSIRZ] = POTENTIALSCC(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    22-Apr-2021 14:56:50

pa_5 = in4(5,:);
pa_7 = in4(7,:);
yf1v1d2 = in3(2,:);
yf1v1d3 = in3(3,:);
yf1v1d4 = in3(4,:);
yf1v1d5 = in3(5,:);
yf1v1d6 = in3(6,:);
yf1v2d3 = in3(10,:);
yf1v2d5 = in3(12,:);
yf1v3d2 = in3(16,:);
yf1v3d3 = in3(17,:);
yf1v3d4 = in3(18,:);
yf1v3d5 = in3(19,:);
yf1v3d6 = in3(20,:);
t2 = -pa_7;
t3 = 1.0./yf1v2d3;
dphir = t3.*yf1v1d3;
if nargout > 1
    t4 = t3.^2;
    t5 = t3.^3;
    t6 = pa_5+t2;
    t7 = 1.0./t6;
    dphiz = -t7.*yf1v1d2;
end
if nargout > 2
    dphirr = t4.*yf1v1d5-t5.*yf1v1d3.*yf1v2d5;
end
if nargout > 3
    t8 = t7.^2;
    dphizz = t8.*yf1v1d4;
end
if nargout > 4
    dphirz = -t3.*t7.*yf1v1d6;
end
if nargout > 5
    dpsir = t3.*yf1v3d3;
end
if nargout > 6
    dpsiz = -t7.*yf1v3d2;
end
if nargout > 7
    dpsirr = t4.*yf1v3d5-t5.*yf1v2d5.*yf1v3d3;
end
if nargout > 8
    dpsizz = t8.*yf1v3d4;
end
if nargout > 9
    dpsirz = -t3.*t7.*yf1v3d6;
end
