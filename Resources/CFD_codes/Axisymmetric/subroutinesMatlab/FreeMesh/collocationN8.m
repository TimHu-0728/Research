function [duxva,duyva,duxxva,duyyva] = collocationN8(in1,in2,in3,in4,in5,in6)
%COLLOCATIONN8
%    [DUXVA,DUYVA,DUXXVA,DUYYVA] = COLLOCATIONN8(IN1,IN2,IN3,IN4,IN5,IN6)

%    This function was generated by the Symbolic Math Toolbox version 7.1.
%    31-Mar-2017 13:18:28

A11 = in6(:,1);
A12 = in6(:,6);
A13 = in6(:,11);
A14 = in6(:,16);
A15 = in6(:,21);
A21 = in6(:,2);
A22 = in6(:,7);
A23 = in6(:,12);
A24 = in6(:,17);
A25 = in6(:,22);
A31 = in6(:,3);
A32 = in6(:,8);
A33 = in6(:,13);
A34 = in6(:,18);
A35 = in6(:,23);
A41 = in6(:,4);
A42 = in6(:,9);
A43 = in6(:,14);
A44 = in6(:,19);
A45 = in6(:,24);
x1r1 = in1(:,1);
x1r2 = in1(:,2);
x1r3 = in1(:,3);
x1r4 = in1(:,4);
x1r5 = in1(:,5);
x1r6 = in1(:,6);
x1r7 = in1(:,7);
x2r1 = in3(:,1);
x2r2 = in3(:,2);
x2r3 = in3(:,3);
x2r4 = in3(:,4);
x2r5 = in3(:,5);
x2r6 = in3(:,6);
x2r7 = in3(:,7);
x1ry1r1 = in5(:,1);
x1ry1r2 = in5(:,2);
x1ry1r3 = in5(:,3);
x1ry1r4 = in5(:,4);
x1ry1r5 = in5(:,5);
x1ry1r6 = in5(:,6);
x1ry1r7 = in5(:,7);
y1r1 = in2(:,1);
y1r2 = in2(:,2);
y1r3 = in2(:,3);
y1r4 = in2(:,4);
y1r5 = in2(:,5);
y1r6 = in2(:,6);
y1r7 = in2(:,7);
y2r1 = in4(:,1);
y2r2 = in4(:,2);
y2r3 = in4(:,3);
y2r4 = in4(:,4);
y2r5 = in4(:,5);
y2r6 = in4(:,6);
y2r7 = in4(:,7);
t2 = x1r1+x1r2+x1r3+x1r4+x1r5+x1r6+x1r7;
t3 = x2r1+x2r2+x2r3+x2r4+x2r5+x2r6+x2r7;
t4 = x1ry1r1+x1ry1r2+x1ry1r3+x1ry1r4+x1ry1r5+x1ry1r6+x1ry1r7;
t5 = y1r1+y1r2+y1r3+y1r4+y1r5+y1r6+y1r7;
t6 = y2r1+y2r2+y2r3+y2r4+y2r5+y2r6+y2r7;
duxva = [-A11.*t2-A13.*t3-A12.*t5-A15.*t4-A14.*t6,A11.*x1r1+A13.*x2r1+A15.*x1ry1r1+A12.*y1r1+A14.*y2r1,A11.*x1r2+A13.*x2r2+A15.*x1ry1r2+A12.*y1r2+A14.*y2r2,A11.*x1r3+A13.*x2r3+A15.*x1ry1r3+A12.*y1r3+A14.*y2r3,A11.*x1r4+A13.*x2r4+A15.*x1ry1r4+A12.*y1r4+A14.*y2r4,A11.*x1r5+A13.*x2r5+A15.*x1ry1r5+A12.*y1r5+A14.*y2r5,A11.*x1r6+A13.*x2r6+A15.*x1ry1r6+A12.*y1r6+A14.*y2r6,A11.*x1r7+A13.*x2r7+A15.*x1ry1r7+A12.*y1r7+A14.*y2r7];
if nargout > 1
    duyva = [-A21.*t2-A23.*t3-A22.*t5-A25.*t4-A24.*t6,A21.*x1r1+A23.*x2r1+A25.*x1ry1r1+A22.*y1r1+A24.*y2r1,A21.*x1r2+A23.*x2r2+A25.*x1ry1r2+A22.*y1r2+A24.*y2r2,A21.*x1r3+A23.*x2r3+A25.*x1ry1r3+A22.*y1r3+A24.*y2r3,A21.*x1r4+A23.*x2r4+A25.*x1ry1r4+A22.*y1r4+A24.*y2r4,A21.*x1r5+A23.*x2r5+A25.*x1ry1r5+A22.*y1r5+A24.*y2r5,A21.*x1r6+A23.*x2r6+A25.*x1ry1r6+A22.*y1r6+A24.*y2r6,A21.*x1r7+A23.*x2r7+A25.*x1ry1r7+A22.*y1r7+A24.*y2r7];
end
if nargout > 2
    duxxva = [A31.*t2.*-2.0-A33.*t3.*2.0-A32.*t5.*2.0-A35.*t4.*2.0-A34.*t6.*2.0,A31.*x1r1.*2.0+A33.*x2r1.*2.0+A35.*x1ry1r1.*2.0+A32.*y1r1.*2.0+A34.*y2r1.*2.0,A31.*x1r2.*2.0+A33.*x2r2.*2.0+A35.*x1ry1r2.*2.0+A32.*y1r2.*2.0+A34.*y2r2.*2.0,A31.*x1r3.*2.0+A33.*x2r3.*2.0+A35.*x1ry1r3.*2.0+A32.*y1r3.*2.0+A34.*y2r3.*2.0,A31.*x1r4.*2.0+A33.*x2r4.*2.0+A35.*x1ry1r4.*2.0+A32.*y1r4.*2.0+A34.*y2r4.*2.0,A31.*x1r5.*2.0+A33.*x2r5.*2.0+A35.*x1ry1r5.*2.0+A32.*y1r5.*2.0+A34.*y2r5.*2.0,A31.*x1r6.*2.0+A33.*x2r6.*2.0+A35.*x1ry1r6.*2.0+A32.*y1r6.*2.0+A34.*y2r6.*2.0,A31.*x1r7.*2.0+A33.*x2r7.*2.0+A35.*x1ry1r7.*2.0+A32.*y1r7.*2.0+A34.*y2r7.*2.0];
end
if nargout > 3
    duyyva = [A41.*t2.*-2.0-A43.*t3.*2.0-A42.*t5.*2.0-A45.*t4.*2.0-A44.*t6.*2.0,A41.*x1r1.*2.0+A43.*x2r1.*2.0+A45.*x1ry1r1.*2.0+A42.*y1r1.*2.0+A44.*y2r1.*2.0,A41.*x1r2.*2.0+A43.*x2r2.*2.0+A45.*x1ry1r2.*2.0+A42.*y1r2.*2.0+A44.*y2r2.*2.0,A41.*x1r3.*2.0+A43.*x2r3.*2.0+A45.*x1ry1r3.*2.0+A42.*y1r3.*2.0+A44.*y2r3.*2.0,A41.*x1r4.*2.0+A43.*x2r4.*2.0+A45.*x1ry1r4.*2.0+A42.*y1r4.*2.0+A44.*y2r4.*2.0,A41.*x1r5.*2.0+A43.*x2r5.*2.0+A45.*x1ry1r5.*2.0+A42.*y1r5.*2.0+A44.*y2r5.*2.0,A41.*x1r6.*2.0+A43.*x2r6.*2.0+A45.*x1ry1r6.*2.0+A42.*y1r6.*2.0+A44.*y2r6.*2.0,A41.*x1r7.*2.0+A43.*x2r7.*2.0+A45.*x1ry1r7.*2.0+A42.*y1r7.*2.0+A44.*y2r7.*2.0];
end
