function [FBBlinea11,dFBBlinea11] = equationFBBlinea11(z1,r1,in3,in4)
%EQUATIONFBBLINEA11
%    [FBBLINEA11,DFBBLINEA11] = EQUATIONFBBLINEA11(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    14-Mar-2021 11:45:35

pa_5 = in4(5,:);
yf1v1d2 = in3(2,:);
yf1v1d3 = in3(3,:);
yf1v2d1 = in3(8,:);
yf1v2d3 = in3(10,:);
yf1v3d1 = in3(15,:);
yf1v3d3 = in3(17,:);
yf1v4d2 = in3(23,:);
yf1v4d3 = in3(24,:);
t2 = abs(yf1v2d3);
t3 = abs(yf1v3d3);
t4 = sign(yf1v2d3);
t5 = sign(yf1v3d3);
t6 = yf1v3d3.^2;
t9 = r1-1.0;
t10 = -yf1v3d1;
t11 = 1.0./yf1v2d1;
t13 = 1.0./yf1v2d3;
t7 = t2.^2;
t8 = t3.^2;
t12 = t11.^2;
t14 = t13.^2;
t15 = pa_5+t10;
t16 = t13.*yf1v1d3;
t18 = t13.*yf1v4d3;
t17 = t14.*yf1v1d3;
t19 = 1.0./t15;
t21 = t7+t8;
t34 = t18.*1.256637061435917e-6;
t20 = t19.^2;
t22 = t19.*yf1v4d2;
t24 = 1.0./sqrt(t21);
t26 = t11.*t19.*yf1v1d2;
t27 = t9.*t13.*t19.*yf1v1d2.*yf1v3d3;
t28 = t9.*t14.*t19.*yf1v1d2.*yf1v3d3;
t23 = -t22;
t25 = t24.^3;
t29 = t9.*t13.*t22.*yf1v3d3;
t30 = t16+t27;
t31 = t17+t28;
t36 = t22.*1.256637061435917e-6;
t38 = t26.*1.256637061435917e-6;
t32 = t11.*t30;
t33 = t18+t26+t29;
t37 = -t36;
t39 = t29.*1.256637061435917e-6;
t35 = t23+t32;
t40 = t32.*1.256637061435917e-6;
t42 = t34+t38+t39;
t41 = t37+t40;
FBBlinea11 = [t24.*t33.*yf1v2d3+t24.*yf1v3d3.*(t22-t32),-yf1v2d1,t10,t24.*t41.*yf1v2d3+t24.*t42.*yf1v3d3];
if nargout > 1
    dFBBlinea11 = reshape([0.0,0.0,0.0,0.0,t11.*t19.*t24.*yf1v2d3-t6.*t9.*t11.*t13.*t19.*t24,0.0,0.0,t11.*t19.*t24.*yf1v3d3.*1.256637061435917e-6+t9.*t11.*t19.*t24.*yf1v3d3.*1.256637061435917e-6,-t11.*t13.*t24.*yf1v3d3,0.0,0.0,t11.*t24.*1.256637061435917e-6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t12.*t24.*t30.*yf1v3d3-t12.*t19.*t24.*yf1v1d2.*yf1v2d3,-1.0,0.0,t12.*t24.*t30.*yf1v2d3.*(-1.256637061435917e-6)-t12.*t19.*t24.*yf1v1d2.*yf1v3d3.*1.256637061435917e-6,0.0,0.0,0.0,0.0,t24.*t33-t24.*yf1v2d3.*(t14.*yf1v4d3+t9.*t14.*t22.*yf1v3d3)+t11.*t24.*t31.*yf1v3d3-t2.*t4.*t25.*t33.*yf1v2d3-t2.*t4.*t25.*yf1v3d3.*(t22-t32),0.0,0.0,t24.*t41-t24.*yf1v3d3.*(t14.*yf1v4d3.*1.256637061435917e-6+t9.*t14.*t36.*yf1v3d3)-t11.*t24.*t31.*yf1v2d3.*1.256637061435917e-6-t2.*t4.*t25.*t41.*yf1v2d3-t2.*t4.*t25.*t42.*yf1v3d3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t24.*yf1v2d3.*(t11.*t20.*yf1v1d2+t9.*t13.*t20.*yf1v3d3.*yf1v4d2)+t24.*yf1v3d3.*(t20.*yf1v4d2-t9.*t11.*t13.*t20.*yf1v1d2.*yf1v3d3),0.0,-1.0,-t24.*yf1v2d3.*(t20.*yf1v4d2.*1.256637061435917e-6-t9.*t11.*t13.*t20.*yf1v1d2.*yf1v3d3.*1.256637061435917e-6)+t24.*yf1v3d3.*(t11.*t20.*yf1v1d2.*1.256637061435917e-6+t9.*t13.*t20.*yf1v3d3.*yf1v4d2.*1.256637061435917e-6),0.0,0.0,0.0,0.0,t24.*(t22-t32)+t9.*t22.*t24-t3.*t5.*t25.*t33.*yf1v2d3-t9.*t13.*t24.*t26.*yf1v3d3-t3.*t5.*t25.*yf1v3d3.*(t22-t32),0.0,0.0,t24.*t39+t24.*t42+t9.*t24.*t38-t3.*t5.*t25.*t41.*yf1v2d3-t3.*t5.*t25.*t42.*yf1v3d3,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t19.*t24.*yf1v3d3+t9.*t19.*t24.*yf1v3d3,0.0,0.0,t19.*t24.*yf1v2d3.*(-1.256637061435917e-6)+t6.*t9.*t13.*t19.*t24.*1.256637061435917e-6,t24,0.0,0.0,t13.*t24.*yf1v3d3.*1.256637061435917e-6,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[4,28]);
end
