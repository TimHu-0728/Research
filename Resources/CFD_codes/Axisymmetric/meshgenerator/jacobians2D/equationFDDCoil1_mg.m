function [FDDCoil1,dFDDCoil1] = equationFDDCoil1_mg(z0,r0,in3,in4)
%EQUATIONFDDCOIL1_MG
%    [FDDCOIL1,DFDDCOIL1] = EQUATIONFDDCOIL1_MG(Z0,R0,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    19-Apr-2021 11:28:18

pa_4 = in4(4,:);
pa_5 = in4(5,:);
yfDv1d1 = in3(1,:);
yfDv2d1 = in3(8,:);
FDDCoil1 = [-pa_4+yfDv1d1,-pa_5+yfDv2d1];
if nargout > 1
    dFDDCoil1 = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[2,14]);
end