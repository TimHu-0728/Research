function out1 = curlH_B(z1,r1,in3,in4)
%CURLH_B
%    OUT1 = CURLH_B(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    22-Apr-2021 17:45:45

pa_5 = in4(5,:);
yf1v1d2 = in3(2,:);
yf1v1d3 = in3(3,:);
yf1v1d4 = in3(4,:);
yf1v1d5 = in3(5,:);
yf1v1d6 = in3(6,:);
yf1v2d1 = in3(8,:);
yf1v2d3 = in3(10,:);
yf1v2d5 = in3(12,:);
yf1v3d1 = in3(15,:);
yf1v3d3 = in3(17,:);
yf1v3d5 = in3(19,:);
t2 = yf1v3d3.^2;
t3 = r1-1.0;
t4 = -yf1v3d1;
t5 = 1.0./yf1v2d1;
t6 = 1.0./yf1v2d3;
t7 = t6.^2;
t8 = pa_5+t4;
t9 = 1.0./t8;
t10 = t9.^2;
out1 = t5.^2.*(t6.*yf1v1d3+t3.*t6.*t9.*yf1v1d2.*yf1v3d3)-t5.*(t6.*(t6.*yf1v1d5+t3.*t6.*t9.*yf1v1d6.*yf1v3d3)+yf1v1d2.*(t6.*(t2.*t3.*t6.*t10+t3.*t6.*t9.*yf1v3d5-t3.*t7.*t9.*yf1v2d5.*yf1v3d3)+t2.*t3.*t7.*t10)-t6.^3.*yf1v1d3.*yf1v2d5+t3.*t6.*t9.*yf1v3d3.*(t6.*yf1v1d6+t3.*t6.*t9.*yf1v1d4.*yf1v3d3))-t5.*t10.*yf1v1d4;
