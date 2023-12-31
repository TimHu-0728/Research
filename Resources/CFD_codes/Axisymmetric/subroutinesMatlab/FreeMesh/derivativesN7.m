function [dux,duy,duxx,duyy] = derivativesN7(in1,in2,in3,in4,in5,in6,in7)
%DERIVATIVESN7
%    [DUX,DUY,DUXX,DUYY] = DERIVATIVESN7(IN1,IN2,IN3,IN4,IN5,IN6,IN7)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    31-Mar-2017 10:43:11

A11 = in7(:,1);
A12 = in7(:,6);
A13 = in7(:,11);
A14 = in7(:,16);
A15 = in7(:,21);
A21 = in7(:,2);
A22 = in7(:,7);
A23 = in7(:,12);
A24 = in7(:,17);
A25 = in7(:,22);
A31 = in7(:,3);
A32 = in7(:,8);
A33 = in7(:,13);
A34 = in7(:,18);
A35 = in7(:,23);
A41 = in7(:,4);
A42 = in7(:,9);
A43 = in7(:,14);
A44 = in7(:,19);
A45 = in7(:,24);
va1 = in1(:,1);
va2 = in1(:,2);
va3 = in1(:,3);
va4 = in1(:,4);
va5 = in1(:,5);
va6 = in1(:,6);
va7 = in1(:,7);
x1r1 = in2(:,1);
x1r2 = in2(:,2);
x1r3 = in2(:,3);
x1r4 = in2(:,4);
x1r5 = in2(:,5);
x1r6 = in2(:,6);
x2r1 = in4(:,1);
x2r2 = in4(:,2);
x2r3 = in4(:,3);
x2r4 = in4(:,4);
x2r5 = in4(:,5);
x2r6 = in4(:,6);
x1ry1r1 = in6(:,1);
x1ry1r2 = in6(:,2);
x1ry1r3 = in6(:,3);
x1ry1r4 = in6(:,4);
x1ry1r5 = in6(:,5);
x1ry1r6 = in6(:,6);
y1r1 = in3(:,1);
y1r2 = in3(:,2);
y1r3 = in3(:,3);
y1r4 = in3(:,4);
y1r5 = in3(:,5);
y1r6 = in3(:,6);
y2r1 = in5(:,1);
y2r2 = in5(:,2);
y2r3 = in5(:,3);
y2r4 = in5(:,4);
y2r5 = in5(:,5);
y2r6 = in5(:,6);
t2 = va1-va2;
t3 = va1-va3;
t4 = va1-va4;
t5 = va1-va5;
t6 = va1-va6;
t7 = va1-va7;
t8 = t2.*x1r1;
t9 = t3.*x1r2;
t10 = t4.*x1r3;
t11 = t5.*x1r4;
t12 = t6.*x1r5;
t13 = t7.*x1r6;
t14 = t8+t9+t10+t11+t12+t13;
t15 = t2.*x2r1;
t16 = t3.*x2r2;
t17 = t4.*x2r3;
t18 = t5.*x2r4;
t19 = t6.*x2r5;
t20 = t7.*x2r6;
t21 = t15+t16+t17+t18+t19+t20;
t22 = t2.*x1ry1r1;
t23 = t3.*x1ry1r2;
t24 = t4.*x1ry1r3;
t25 = t5.*x1ry1r4;
t26 = t6.*x1ry1r5;
t27 = t7.*x1ry1r6;
t28 = t22+t23+t24+t25+t26+t27;
t29 = t2.*y1r1;
t30 = t3.*y1r2;
t31 = t4.*y1r3;
t32 = t5.*y1r4;
t33 = t6.*y1r5;
t34 = t7.*y1r6;
t35 = t29+t30+t31+t32+t33+t34;
t36 = t2.*y2r1;
t37 = t3.*y2r2;
t38 = t4.*y2r3;
t39 = t5.*y2r4;
t40 = t6.*y2r5;
t41 = t7.*y2r6;
t42 = t36+t37+t38+t39+t40+t41;
dux = -A11.*t14-A13.*t21-A15.*t28-A12.*t35-A14.*t42;
if nargout > 1
    duy = -A21.*t14-A23.*t21-A25.*t28-A22.*t35-A24.*t42;
end
if nargout > 2
    duxx = A31.*t14.*-2.0-A33.*t21.*2.0-A35.*t28.*2.0-A32.*t35.*2.0-A34.*t42.*2.0;
end
if nargout > 3
    duyy = A41.*t14.*-2.0-A43.*t21.*2.0-A45.*t28.*2.0-A42.*t35.*2.0-A44.*t42.*2.0;
end
