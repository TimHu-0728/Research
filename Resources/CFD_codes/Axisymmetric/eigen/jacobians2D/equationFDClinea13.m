function [FDClinea13,dFDClinea13] = equationFDClinea13(z1,r1,in3,in4)
%EQUATIONFDCLINEA13
%    [FDCLINEA13,DFDCLINEA13] = EQUATIONFDCLINEA13(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    27-Apr-2021 15:31:14

pa_5 = in4(5,:);
pa_7 = in4(7,:);
yf1v1d3 = in3(3,:);
yf1v2d1 = in3(8,:);
yf1v2d3 = in3(10,:);
yf1v3d3 = in3(17,:);
t2 = 1.0./yf1v2d3;
FDClinea13 = [t2.*yf1v1d3,t2.*yf1v3d3,yf1v2d1,pa_5-r1.*(pa_5-pa_7),0.0];
if nargout > 1
    t3 = t2.^2;
    dFDClinea13 = reshape([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t3.*yf1v1d3,-t3.*yf1v3d3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[5,21]);
end