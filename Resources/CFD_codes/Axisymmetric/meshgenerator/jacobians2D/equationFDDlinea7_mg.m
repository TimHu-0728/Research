function [FDDlinea7,dFDDlinea7] = equationFDDlinea7_mg(z0,r0,in3,in4)
%EQUATIONFDDLINEA7_MG
%    [FDDLINEA7,DFDDLINEA7] = EQUATIONFDDLINEA7_MG(Z0,R0,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    19-Apr-2021 11:28:19

pa_2 = in4(2,:);
pa_4 = in4(4,:);
yfDv1d1 = in3(1,:);
yfDv2d1 = in3(8,:);
FDDlinea7 = [-pa_4+yfDv1d1,-pa_2+yfDv2d1];
if nargout > 1
    dFDDlinea7 = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[2,14]);
end
