function [FDAlinea20,dFDAlinea20] = equationFDAlinea20(z1,r1,in3,in4)
%EQUATIONFDALINEA20
%    [FDALINEA20,DFDALINEA20] = EQUATIONFDALINEA20(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    03-Oct-2021 10:16:02

yf1v4d2 = in3(23,:);
yf1v4d3 = in3(24,:);
yf1v5d1 = in3(29,:);
yf1v5d3 = in3(31,:);
yf1v6d1 = in3(36,:);
yf1v6d3 = in3(38,:);
yf1v7d1 = in3(43,:);
yf1v9d2 = in3(58,:);
yf1v9d3 = in3(59,:);
yf1v11d1 = in3(71,:);
t2 = 1.0./yf1v5d3;
t3 = yf1v6d1-yf1v11d1;
t4 = 1.0./t3;
t5 = 1.0./yf1v5d1;
t6 = yf1v7d1.*1.256637061435917e-6;
t7 = t6+1.256637061435917e-6;
t8 = 1.0./yf1v5d3.^2;
t9 = 1.0./t3.^2;
t10 = t5.*t9.*yf1v4d2;
t11 = t10-r1.*t2.*t9.*yf1v6d3.*yf1v9d2;
t12 = t7.*t11;
t13 = 1.0./yf1v5d1.^2;
t14 = t2.*yf1v4d3;
t15 = t14-r1.*t2.*t4.*yf1v4d2.*yf1v6d3;
FDAlinea20 = [-t7.*(t2.*yf1v9d3+t4.*t5.*yf1v4d2-r1.*t2.*t4.*yf1v6d3.*yf1v9d2),t5.*t15-t4.*yf1v9d2,yf1v5d1,r1.*t3,0.0];
if nargout > 1
    t16 = t9.*yf1v9d2;
    t17 = r1.*t2.*t5.*t9.*yf1v4d2.*yf1v6d3;
    dFDAlinea20 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t4.*t5.*t7,-r1.*t2.*t4.*t5.*yf1v6d3,0.0,0.0,0.0,0.0,t2.*t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t4.*t7.*t13.*yf1v4d2,-t13.*t15,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t7.*(t8.*yf1v9d3-r1.*t4.*t8.*yf1v6d3.*yf1v9d2),-t5.*(t8.*yf1v4d3-r1.*t4.*t8.*yf1v4d2.*yf1v6d3),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t12,t16+t17,0.0,r1,0.0,0.0,0.0,0.0,0.0,0.0,r1.*t2.*t4.*t7.*yf1v9d2,-r1.*t2.*t4.*t5.*yf1v4d2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2.*yf1v9d3.*(-1.256637061435917e-6)-t4.*t5.*yf1v4d2.*1.256637061435917e-6+r1.*t2.*t4.*yf1v6d3.*yf1v9d2.*1.256637061435917e-6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,r1.*t2.*t4.*t7.*yf1v6d3,-t4,0.0,0.0,0.0,-t2.*t7,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t12,-t16-t17,0.0,-r1,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[5,77]);
end
