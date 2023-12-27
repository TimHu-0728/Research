function [FDDa,dFDDa] = equationFDDa_mg(z0,r0,in3,in4)
%EQUATIONFDDA_MG
%    [FDDA,DFDDA] = EQUATIONFDDA_MG(Z0,R0,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    19-Apr-2021 11:28:18

yfDv1d4 = in3(4,:);
yfDv1d5 = in3(5,:);
yfDv2d4 = in3(11,:);
yfDv2d5 = in3(12,:);
FDDa = [yfDv1d4+yfDv1d5,yfDv2d4+yfDv2d5];
if nargout > 1
    dFDDa = reshape([0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,1.0,0.0,0.0,0.0,0.0],[2,14]);
end