function [FEE,dFEE] = equationFEE(z1,r1,in3,in4,phi_ref)
%EQUATIONFEE
%    [FEE,DFEE] = EQUATIONFEE(Z1,R1,IN3,IN4,PHI_REF)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    09-May-2021 10:18:52

pa_30 = in4(30,:);
yf1v1d1 = in3(1,:);
yf1v1d2 = in3(2,:);
yf1v2d1 = in3(8,:);
yf1v2d3 = in3(10,:);
yf1v2d4 = in3(11,:);
yf1v2d5 = in3(12,:);
yf1v3d3 = in3(17,:);
yf1v3d4 = in3(18,:);
yf1v3d5 = in3(19,:);
t2 = 1.0./pa_30;
t4 = 1.0./yf1v2d1;
t6 = 1.0./yf1v2d3;
t3 = t2.^2;
t5 = t4.^2;
t7 = t6.^2;
t8 = t6.^3;
t9 = t6.*yf1v3d3.*1.256637061435917e-6;
t10 = t2.*t4.*yf1v1d2.*1.256637061435917e-6;
t11 = -t10;
t12 = t9+t11;
FEE = [-t4.*t12-t3.*yf1v3d4.*1.256637061435917e-6-t7.*yf1v3d5.*1.256637061435917e-6-t2.*t5.*yf1v1d2.*1.256637061435917e-6+t8.*yf1v2d5.*yf1v3d3.*1.256637061435917e-6,yf1v2d4,-phi_ref+yf1v1d1];
if nargout > 1
    dFEE = reshape([0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t5.*t12+t4.^2.*t10,0.0,0.0,0.0,0.0,0.0,t8.*yf1v3d5.*2.513274122871835e-6-t7.^2.*yf1v2d5.*yf1v3d3.*3.769911184307752e-6+t4.*t7.*yf1v3d3.*1.256637061435917e-6,0.0,0.0,0.0,1.0,0.0,t8.*yf1v3d3.*1.256637061435917e-6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t4.*t6.*(-1.256637061435917e-6)+t8.*yf1v2d5.*1.256637061435917e-6,0.0,0.0,t3.*(-1.256637061435917e-6),0.0,0.0,t7.*(-1.256637061435917e-6),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,21]);
end
