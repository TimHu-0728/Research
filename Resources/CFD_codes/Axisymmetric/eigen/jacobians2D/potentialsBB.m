function [dphir,dphiz,dphirr,dphizz,dphirz,dpsir,dpsiz,dpsirr,dpsizz,dpsirz] = potentialsBB(z1,r1,in3,in4)
%POTENTIALSBB
%    [DPHIR,DPHIZ,DPHIRR,DPHIZZ,DPHIRZ,DPSIR,DPSIZ,DPSIRR,DPSIZZ,DPSIRZ] = POTENTIALSBB(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    22-Apr-2021 14:56:30

pa_5 = in4(5,:);
yf1v1d2 = in3(2,:);
yf1v1d3 = in3(3,:);
yf1v1d4 = in3(4,:);
yf1v1d5 = in3(5,:);
yf1v1d6 = in3(6,:);
yf1v2d3 = in3(10,:);
yf1v2d5 = in3(12,:);
yf1v3d1 = in3(15,:);
yf1v3d3 = in3(17,:);
yf1v3d5 = in3(19,:);
yf1v4d2 = in3(23,:);
yf1v4d3 = in3(24,:);
yf1v4d4 = in3(25,:);
yf1v4d5 = in3(26,:);
yf1v4d6 = in3(27,:);
t2 = yf1v3d3.^2;
t3 = r1-1.0;
t4 = -yf1v3d1;
t5 = 1.0./yf1v2d3;
t6 = t5.^2;
t7 = t5.^3;
t8 = pa_5+t4;
t9 = 1.0./t8;
dphir = t5.*yf1v1d3+t3.*t5.*t9.*yf1v1d2.*yf1v3d3;
if nargout > 1
    dphiz = t9.*yf1v1d2;
end
if nargout > 2
    t10 = t9.^2;
    t11 = t3.*t5.*t9.*yf1v3d5;
    t12 = t3.*t6.*t9.*yf1v2d5.*yf1v3d3;
    t13 = t2.*t3.*t5.*t10;
    t14 = t2.*t3.*t6.*t10;
    t15 = -t12;
    t16 = t11+t13+t15;
    t17 = t5.*t16;
    t18 = t14+t17;
    dphirr = t5.*(t5.*yf1v1d5+t3.*t5.*t9.*yf1v1d6.*yf1v3d3)+t18.*yf1v1d2-t7.*yf1v1d3.*yf1v2d5+t3.*t5.*t9.*yf1v3d3.*(t5.*yf1v1d6+t3.*t5.*t9.*yf1v1d4.*yf1v3d3);
end
if nargout > 3
    dphizz = t10.*yf1v1d4;
end
if nargout > 4
    dphirz = t5.*t9.*yf1v1d6+t5.*t10.*yf1v1d2.*yf1v3d3+t3.*t5.*t10.*yf1v1d4.*yf1v3d3;
end
if nargout > 5
    dpsir = t5.*yf1v4d3+t3.*t5.*t9.*yf1v3d3.*yf1v4d2;
end
if nargout > 6
    dpsiz = t9.*yf1v4d2;
end
if nargout > 7
    dpsirr = t5.*(t5.*yf1v4d5+t3.*t5.*t9.*yf1v3d3.*yf1v4d6)+t18.*yf1v4d2-t7.*yf1v2d5.*yf1v4d3+t3.*t5.*t9.*yf1v3d3.*(t5.*yf1v4d6+t3.*t5.*t9.*yf1v3d3.*yf1v4d4);
end
if nargout > 8
    dpsizz = t10.*yf1v4d4;
end
if nargout > 9
    dpsirz = t5.*t9.*yf1v4d6+t5.*t10.*yf1v3d3.*yf1v4d2+t3.*t5.*t10.*yf1v3d3.*yf1v4d4;
end