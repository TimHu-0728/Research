function [FEElinea1,dFEElinea1] = equationFEElinea1(z1,r1,in3,in4)
%EQUATIONFEELINEA1
%    [FEELINEA1,DFEELINEA1] = EQUATIONFEELINEA1(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    16-Feb-2021 10:52:42

pa_5 = in4(5,:);
pa_7 = in4(7,:);
pa_15 = in4(15,:);
pa_16 = in4(16,:);
pa_17 = in4(17,:);
pa_18 = in4(18,:);
yf1v1d1 = in3(1,:);
yf1v2d1 = in3(8,:);
yf1v2d2 = in3(9,:);
yf1v3d2 = in3(16,:);
t2 = yf1v2d1.^2;
t3 = 1.0./pi;
t4 = -pa_7;
t5 = r1.*(3.0./5.0);
t6 = -t2;
t7 = pa_5+t4;
t8 = -t5;
t9 = 1.0./t7;
t10 = pa_16+t8+3.0./5.0;
t11 = pa_18+t8+3.0./5.0;
t12 = t10.^2;
t13 = t11.^2;
t14 = t12.*2.0;
t15 = t13.*2.0;
t16 = t2+t12;
t17 = t2+t13;
t18 = 1.0./t16.^(5.0./2.0);
t19 = 1.0./t17.^(5.0./2.0);
FEElinea1 = [-t9.*yf1v3d2-(pa_15.*t3.*t18.*(t2-t14))./4.0-(pa_17.*t3.*t19.*(t2-t15))./4.0,yf1v2d2,yf1v1d1];
if nargout > 1
    dFEElinea1 = reshape([0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,pa_15.*t3.*t18.*yf1v2d1.*(-1.0./2.0)-(pa_17.*t3.*t19.*yf1v2d1)./2.0+pa_15.*t3.*1.0./t16.^(7.0./2.0).*yf1v2d1.*(t2-t14).*(5.0./4.0)+pa_17.*t3.*1.0./t17.^(7.0./2.0).*yf1v2d1.*(t2-t15).*(5.0./4.0),0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t9,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,21]);
end
