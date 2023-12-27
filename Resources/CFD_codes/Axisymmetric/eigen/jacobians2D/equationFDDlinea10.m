function [FDDlinea10,dFDDlinea10] = equationFDDlinea10(z0,r0,in3,in4)
%EQUATIONFDDLINEA10
%    [FDDLINEA10,DFDDLINEA10] = EQUATIONFDDLINEA10(Z0,R0,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    27-Apr-2021 08:58:33

yfDv1d2 = in3(2,:);
yfDv1d3 = in3(3,:);
yfDv2d2 = in3(9,:);
yfDv2d3 = in3(10,:);
yfDv3d1 = in3(15,:);
yfDv3d2 = in3(16,:);
yfDv3d3 = in3(17,:);
yfDv4d1 = in3(22,:);
yfDv4d2 = in3(23,:);
yfDv4d3 = in3(24,:);
yfDv5d1 = in3(29,:);
t2 = yfDv3d2.*yfDv4d3;
t3 = yfDv3d3.*yfDv4d2;
t4 = yfDv4d2.^2;
t5 = yfDv4d3.^2;
t6 = -t3;
t7 = t2+t6;
t8 = 1.0./t7;
t9 = t8.^2;
t10 = t8.*yfDv4d2;
t11 = t8.*yfDv4d3;
t12 = -t11;
FDDlinea10 = [t10.*yfDv1d3+t12.*yfDv1d2,t10.*yfDv2d3+t12.*yfDv2d2,-yfDv3d1,-yfDv4d1,yfDv5d1];
if nargout > 1
    dFDDlinea10 = reshape([0.0,0.0,0.0,0.0,0.0,t12,0.0,0.0,0.0,0.0,t10,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t12,0.0,0.0,0.0,0.0,t10,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,0.0,t5.*t9.*yfDv1d2-t9.*yfDv1d3.*yfDv4d2.*yfDv4d3,t5.*t9.*yfDv2d2-t9.*yfDv2d3.*yfDv4d2.*yfDv4d3,0.0,0.0,0.0,t4.*t9.*yfDv1d3-t9.*yfDv1d2.*yfDv4d2.*yfDv4d3,t4.*t9.*yfDv2d3-t9.*yfDv2d2.*yfDv4d2.*yfDv4d3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,t8.*yfDv1d3+t3.*t9.*yfDv1d3-t9.*yfDv1d2.*yfDv3d3.*yfDv4d3,t8.*yfDv2d3+t3.*t9.*yfDv2d3-t9.*yfDv2d2.*yfDv3d3.*yfDv4d3,0.0,0.0,0.0,-t8.*yfDv1d2+t2.*t9.*yfDv1d2-t9.*yfDv1d3.*yfDv3d2.*yfDv4d2,-t8.*yfDv2d2+t2.*t9.*yfDv2d2-t9.*yfDv2d3.*yfDv3d2.*yfDv4d2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[5,35]);
end
