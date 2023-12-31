function [FDDlinea20,dFDDlinea20] = equationFDDlinea20(z0,r0,in3,in4)
%EQUATIONFDDLINEA20
%    [FDDLINEA20,DFDDLINEA20] = EQUATIONFDDLINEA20(Z0,R0,IN3,IN4)

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
t4 = yfDv3d2.^2;
t5 = yfDv3d3.^2;
t6 = yfDv4d2.^2;
t7 = yfDv4d3.^2;
t8 = 1.0./yfDv3d1;
t9 = t8.^2;
t10 = -t3;
t11 = t2+t10;
t12 = 1.0./t11;
t13 = t12.^2;
t14 = t12.*yfDv1d2;
t15 = t12.*yfDv1d3;
t16 = t14.*yfDv3d3;
t17 = t15.*yfDv3d2;
t18 = t14.*yfDv4d3;
t19 = t15.*yfDv4d2;
t20 = t13.*yfDv1d3.*yfDv3d2.*yfDv4d2;
t21 = t13.*yfDv1d2.*yfDv3d3.*yfDv4d3;
t22 = -t17;
t23 = -t19;
t24 = t16+t22;
t25 = t18+t23;
FDDlinea20 = [t8.*t24.*(-1.256637061435917e-6)+t12.*yfDv2d2.*yfDv4d3.*1.256637061435917e-6-t12.*yfDv2d3.*yfDv4d2.*1.256637061435917e-6,-t8.*t25-t12.*yfDv2d2.*yfDv3d3+t12.*yfDv2d3.*yfDv3d2,-yfDv3d1,-yfDv4d1,yfDv5d1];
if nargout > 1
    dFDDlinea20 = reshape([0.0,0.0,0.0,0.0,0.0,t8.*t12.*yfDv3d3.*(-1.256637061435917e-6),-t8.*t12.*yfDv4d3,0.0,0.0,0.0,t8.*t12.*yfDv3d2.*1.256637061435917e-6,t8.*t12.*yfDv4d2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t12.*yfDv4d3.*1.256637061435917e-6,-t12.*yfDv3d3,0.0,0.0,0.0,t12.*yfDv4d2.*(-1.256637061435917e-6),t12.*yfDv3d2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t9.*t24.*1.256637061435917e-6,t9.*t25,-1.0,0.0,0.0,t8.*(t15+t21-t2.*t13.*yfDv1d3).*1.256637061435917e-6-t7.*t13.*yfDv2d2.*1.256637061435917e-6+t13.*yfDv2d3.*yfDv4d2.*yfDv4d3.*1.256637061435917e-6,t8.*(t7.*t13.*yfDv1d2-t13.*yfDv1d3.*yfDv4d2.*yfDv4d3)+t12.*yfDv2d3-t2.*t13.*yfDv2d3+t13.*yfDv2d2.*yfDv3d3.*yfDv4d3,0.0,0.0,0.0,t8.*(t14-t20+t3.*t13.*yfDv1d2).*(-1.256637061435917e-6)-t6.*t13.*yfDv2d3.*1.256637061435917e-6+t13.*yfDv2d2.*yfDv4d2.*yfDv4d3.*1.256637061435917e-6,t8.*(t6.*t13.*yfDv1d3-t13.*yfDv1d2.*yfDv4d2.*yfDv4d3)-t12.*yfDv2d2+t10.*t13.*yfDv2d2+t13.*yfDv2d3.*yfDv3d2.*yfDv4d2,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-1.0,0.0,t8.*(t5.*t13.*yfDv1d2-t13.*yfDv1d3.*yfDv3d2.*yfDv3d3).*(-1.256637061435917e-6)-t12.*yfDv2d3.*1.256637061435917e-6-t3.*t13.*yfDv2d3.*1.256637061435917e-6+t13.*yfDv2d2.*yfDv3d3.*yfDv4d3.*1.256637061435917e-6,t8.*(t15-t21+t3.*t13.*yfDv1d3)-t5.*t13.*yfDv2d2+t13.*yfDv2d3.*yfDv3d2.*yfDv3d3,0.0,0.0,0.0,t8.*(t4.*t13.*yfDv1d3-t13.*yfDv1d2.*yfDv3d2.*yfDv3d3).*(-1.256637061435917e-6)+t12.*yfDv2d2.*1.256637061435917e-6-t2.*t13.*yfDv2d2.*1.256637061435917e-6+t13.*yfDv2d3.*yfDv3d2.*yfDv4d2.*1.256637061435917e-6,-t8.*(t14+t20-t2.*t13.*yfDv1d2)-t4.*t13.*yfDv2d3+t13.*yfDv2d2.*yfDv3d2.*yfDv3d3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[5,35]);
end
