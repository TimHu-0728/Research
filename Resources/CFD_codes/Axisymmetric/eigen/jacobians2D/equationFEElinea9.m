function [FEElinea9,dFEElinea9] = equationFEElinea9(z1,r1,in3,in4)
%EQUATIONFEELINEA9
%    [FEELINEA9,DFEELINEA9] = EQUATIONFEELINEA9(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    09-May-2021 10:18:53

pa_30 = in4(30,:);
yf1v1d2 = in3(2,:);
yf1v1d3 = in3(3,:);
yf1v2d1 = in3(8,:);
yf1v2d3 = in3(10,:);
yf1v3d2 = in3(16,:);
yf1v3d3 = in3(17,:);
t2 = 1.0./pa_30;
t3 = 1.0./yf1v2d1;
t5 = 1.0./yf1v2d3;
FEElinea9 = [t5.*yf1v3d3-t2.*t3.*yf1v1d2,-yf1v2d1,t2.*yf1v3d2.*(-1.256637061435917e-6)-t3.*t5.*yf1v1d3.*1.256637061435917e-6];
if nargout > 1
    t4 = t3.^2;
    t6 = t5.^2;
    dFEElinea9 = reshape([0.0,0.0,0.0,-t2.*t3,0.0,0.0,0.0,0.0,t3.*t5.*(-1.256637061435917e-6),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2.*t4.*yf1v1d2,-1.0,t4.*t5.*yf1v1d3.*1.256637061435917e-6,0.0,0.0,0.0,-t6.*yf1v3d3,0.0,t3.*t6.*yf1v1d3.*1.256637061435917e-6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2.*(-1.256637061435917e-6),t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,21]);
end
