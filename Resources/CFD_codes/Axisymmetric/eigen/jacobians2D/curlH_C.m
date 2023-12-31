function out1 = curlH_C(z1,r1,in3,in4)
%CURLH_C
%    OUT1 = CURLH_C(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    22-Apr-2021 17:45:58

pa_5 = in4(5,:);
pa_7 = in4(7,:);
yf1v1d3 = in3(3,:);
yf1v1d4 = in3(4,:);
yf1v1d5 = in3(5,:);
yf1v2d1 = in3(8,:);
yf1v2d3 = in3(10,:);
yf1v2d5 = in3(12,:);
t2 = 1.0./yf1v2d1;
out1 = -t2.*(yf1v1d5.*1.0./yf1v2d3.^2-yf1v1d3.*1.0./yf1v2d3.^3.*yf1v2d5)-t2.*yf1v1d4.*1.0./(pa_5-pa_7).^2+(t2.^2.*yf1v1d3)./yf1v2d3;
