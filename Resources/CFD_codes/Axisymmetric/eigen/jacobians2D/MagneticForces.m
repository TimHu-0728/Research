function [Forcer,Forcez] = MagneticForces(z1,r1,in3,in4)
%MAGNETICFORCES
%    [FORCER,FORCEZ] = MAGNETICFORCES(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    03-Oct-2021 10:13:10

yf1v4d2 = in3(23,:);
yf1v4d3 = in3(24,:);
yf1v4d4 = in3(25,:);
yf1v4d5 = in3(26,:);
yf1v4d6 = in3(27,:);
yf1v5d1 = in3(29,:);
yf1v5d3 = in3(31,:);
yf1v5d5 = in3(33,:);
yf1v6d1 = in3(36,:);
yf1v6d3 = in3(38,:);
yf1v6d5 = in3(40,:);
yf1v7d1 = in3(43,:);
yf1v9d2 = in3(58,:);
yf1v9d3 = in3(59,:);
yf1v9d4 = in3(60,:);
yf1v9d5 = in3(61,:);
yf1v9d6 = in3(62,:);
yf1v11d1 = in3(71,:);
t2 = 1.0./yf1v5d3;
t3 = yf1v6d1-yf1v11d1;
t4 = 1.0./t3;
t5 = 1.0./yf1v5d1;
t6 = 1.0./t3.^2;
t7 = 1.0./yf1v5d3.^2;
t8 = yf1v6d3.^2;
t9 = t2.*yf1v4d3;
t28 = r1.*t2.*t4.*yf1v4d2.*yf1v6d3;
t10 = t9-t28;
t11 = t5.*t10;
t12 = t2.*t6.*yf1v4d2.*yf1v6d3;
t13 = r1.*t2.*t6.*yf1v4d4.*yf1v6d3;
t14 = t12+t13-t2.*t4.*yf1v4d6;
t15 = t5.*t14;
t16 = t11-t4.*yf1v9d2;
t17 = t2.*yf1v9d3;
t18 = t4.*t5.*yf1v4d2;
t19 = t17+t18-r1.*t2.*t4.*yf1v6d3.*yf1v9d2;
t20 = r1.*t2.*t6.*t8;
t21 = r1.*t4.*t7.*yf1v5d5.*yf1v6d3;
t22 = t20+t21-r1.*t2.*t4.*yf1v6d5;
t23 = t2.*t22;
t24 = r1.*t6.*t7.*t8;
t25 = t23+t24;
t26 = 1.0./yf1v5d3.^3;
t27 = 1.0./yf1v5d1.^2;
t29 = t2.*t4.*yf1v9d6;
Forcer = t19.*yf1v7d1.*(t15-t2.*(t2.*yf1v9d5-r1.*t2.*t4.*yf1v6d3.*yf1v9d6)-t25.*yf1v9d2+t4.*t27.*yf1v4d2+t26.*yf1v5d5.*yf1v9d3+r1.*t2.*t4.*yf1v6d3.*(t2.*yf1v9d6-r1.*t2.*t4.*yf1v6d3.*yf1v9d4)).*(-1.256637061435917e-6)-t16.*yf1v7d1.*(t29+t5.*t6.*yf1v4d4-t2.*t6.*yf1v6d3.*yf1v9d2-r1.*t2.*t6.*yf1v6d3.*yf1v9d4).*1.256637061435917e-6;
if nargout > 1
    Forcez = t16.*yf1v7d1.*(t15+t6.*yf1v9d4).*(-1.256637061435917e-6)-t19.*yf1v7d1.*(-t29-t10.*t27+t5.*(t2.*(t2.*yf1v4d5-r1.*t2.*t4.*yf1v4d6.*yf1v6d3)+t25.*yf1v4d2-t26.*yf1v4d3.*yf1v5d5-r1.*t2.*t4.*yf1v6d3.*(t2.*yf1v4d6-r1.*t2.*t4.*yf1v4d4.*yf1v6d3))+t2.*t6.*yf1v6d3.*yf1v9d2+r1.*t2.*t6.*yf1v6d3.*yf1v9d4).*1.256637061435917e-6;
end
