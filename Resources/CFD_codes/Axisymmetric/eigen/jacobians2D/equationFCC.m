function [FCC,dFCC] = equationFCC(z1,r1,in3,in4)
%EQUATIONFCC
%    [FCC,DFCC] = EQUATIONFCC(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    18-Apr-2021 19:04:27

pa_5 = in4(5,:);
pa_7 = in4(7,:);
yf1v1d2 = in3(2,:);
yf1v1d3 = in3(3,:);
yf1v1d4 = in3(4,:);
yf1v1d5 = in3(5,:);
yf1v2d1 = in3(8,:);
yf1v2d3 = in3(10,:);
yf1v2d4 = in3(11,:);
yf1v2d5 = in3(12,:);
yf1v3d3 = in3(17,:);
yf1v3d4 = in3(18,:);
yf1v3d5 = in3(19,:);
t2 = -pa_7;
t3 = 1.0./yf1v2d1;
t6 = 1.0./yf1v2d3;
t4 = t3.^2;
t5 = t3.^3;
t7 = t6.^2;
t8 = t6.^3;
t10 = pa_5+t2;
t17 = t6.*yf1v3d3.*1.256637061435917e-6;
t9 = t7.^2;
t11 = t7.*yf1v1d5;
t12 = t8.*yf1v1d3.*yf1v2d5;
t14 = 1.0./t10;
t13 = -t12;
t15 = t14.^2;
t18 = t3.*t14.*yf1v1d2.*1.256637061435917e-6;
t16 = t11+t13;
t19 = -t18;
t20 = t17+t19;
FCC = [-t3.*t20-t7.*yf1v3d5.*1.256637061435917e-6-t15.*yf1v3d4.*1.256637061435917e-6-t4.*t14.*yf1v1d2.*1.256637061435917e-6+t8.*yf1v2d5.*yf1v3d3.*1.256637061435917e-6,yf1v2d4,-t3.*t16+t4.*t6.*yf1v1d3-t3.*t15.*yf1v1d4];
if nargout > 1
    dFCC = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t4.*t6+t3.*t8.*yf1v2d5,0.0,0.0,-t3.*t15,0.0,0.0,-t3.*t7,0.0,0.0,0.0,0.0,0.0,0.0,t4.*t20+t5.*t14.*yf1v1d2.*1.256637061435917e-6,0.0,t4.*t16-t5.*t6.*yf1v1d3.*2.0+t4.*t15.*yf1v1d4,0.0,0.0,0.0,t8.*yf1v3d5.*2.513274122871835e-6+t3.*t7.*yf1v3d3.*1.256637061435917e-6-t9.*yf1v2d5.*yf1v3d3.*3.769911184307752e-6,0.0,t3.*(t8.*yf1v1d5.*2.0-t9.*yf1v1d3.*yf1v2d5.*3.0)-t4.*t7.*yf1v1d3,0.0,1.0,0.0,t8.*yf1v3d3.*1.256637061435917e-6,0.0,t3.*t8.*yf1v1d3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t3.*t6.*(-1.256637061435917e-6)+t8.*yf1v2d5.*1.256637061435917e-6,0.0,0.0,t15.*(-1.256637061435917e-6),0.0,0.0,t7.*(-1.256637061435917e-6),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,21]);
end
