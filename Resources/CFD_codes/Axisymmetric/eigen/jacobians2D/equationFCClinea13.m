function [FCClinea13,dFCClinea13] = equationFCClinea13(z1,r1,in3,in4)
%EQUATIONFCCLINEA13
%    [FCCLINEA13,DFCCLINEA13] = EQUATIONFCCLINEA13(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    18-Apr-2021 19:04:28

pa_6 = in4(6,:);
yf1v1d1 = in3(1,:);
yf1v2d1 = in3(8,:);
yf1v3d1 = in3(15,:);
FCClinea13 = [yf1v1d1,-pa_6+yf1v2d1,yf1v3d1];
if nargout > 1
    dFCClinea13 = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[3,21]);
end
