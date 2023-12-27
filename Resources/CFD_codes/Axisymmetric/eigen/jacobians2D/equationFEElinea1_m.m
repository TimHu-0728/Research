function [FEElinea1_m,dFEElinea1_m] = equationFEElinea1_m(z1,r1,in3,in4)
%EQUATIONFEELINEA1_M
%    [FEELINEA1_M,DFEELINEA1_M] = EQUATIONFEELINEA1_M(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    09-May-2021 10:18:53

pa_15 = in4(15,:);
pa_16 = in4(16,:);
pa_17 = in4(17,:);
pa_18 = in4(18,:);
pa_30 = in4(30,:);
yf1v1d1 = in3(1,:);
yf1v2d1 = in3(8,:);
yf1v2d2 = in3(9,:);
yf1v3d1 = in3(15,:);
t2 = pa_30.*r1;
t3 = yf1v2d1.^2;
t4 = 1.0./pi;
t5 = -pa_30;
t6 = pa_16+t2+t5;
t7 = pa_18+t2+t5;
t8 = t6.^2;
t9 = t7.^2;
t10 = t3+t8;
t11 = t3+t9;
FEElinea1_m = [yf1v3d1+(pa_15.*t4.*t6.*1.0./t10.^(3.0./2.0))./4.0+(pa_17.*t4.*t7.*1.0./t11.^(3.0./2.0))./4.0,yf1v2d2,yf1v1d1];
if nargout > 1
    dFEElinea1_m = reshape([0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,pa_15.*t4.*t6.*1.0./t10.^(5.0./2.0).*yf1v2d1.*(-3.0./4.0)-pa_17.*t4.*t7.*1.0./t11.^(5.0./2.0).*yf1v2d1.*(3.0./4.0),0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,21]);
end
