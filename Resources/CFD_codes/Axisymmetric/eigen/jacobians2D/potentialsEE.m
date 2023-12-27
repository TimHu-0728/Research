function [dphir,dphiz,dphirr,dphizz,dphirz,dpsir,dpsiz,dpsirr,dpsizz,dpsirz] = potentialsEE(z1,r1,in3,in4)
%POTENTIALSEE
%    [DPHIR,DPHIZ,DPHIRR,DPHIZZ,DPHIRZ,DPSIR,DPSIZ,DPSIRR,DPSIZZ,DPSIRZ] = POTENTIALSEE(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    09-May-2021 10:18:51

pa_30 = in4(30,:);
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
t2 = 1.0./pa_30;
t4 = 1.0./yf1v2d3;
dphir = t4.*yf1v1d3;
if nargout > 1
    dphiz = -t2.*yf1v1d2;
end
if nargout > 2
    t3 = t2.^2;
    t5 = t4.^2;
    t6 = t4.^3;
    dphirr = t5.*yf1v1d5-t6.*yf1v1d3.*yf1v2d5;
end
if nargout > 3
    dphizz = t3.*yf1v1d4;
end
if nargout > 4
    dphirz = -t2.*t4.*yf1v1d6;
end
if nargout > 5
    dpsir = t4.*yf1v3d3;
end
if nargout > 6
    dpsiz = -t2.*yf1v3d2;
end
if nargout > 7
    dpsirr = t5.*yf1v3d5-t6.*yf1v2d5.*yf1v3d3;
end
if nargout > 8
    dpsizz = t3.*yf1v3d4;
end
if nargout > 9
    dpsirz = -t2.*t4.*yf1v3d6;
end
