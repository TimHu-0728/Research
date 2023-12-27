function [duxva,duyva,duxxva,duyyva] = collocationtalylorN6(in1,in2,in3,in4,in5,in6)
%COLLOCATIONTALYLORN6
%    [DUXVA,DUYVA,DUXXVA,DUYYVA] = COLLOCATIONTALYLORN6(IN1,IN2,IN3,IN4,IN5,IN6)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    14-May-2020 16:41:20

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
duxva = [-A11-A12-A13-A14-A15,A11,A12,A13,A14,A15];
if nargout > 1
    duyva = [-A21-A22-A23-A24-A25,A21,A22,A23,A24,A25];
end
if nargout > 2
    duxxva = [A31.*-2.0-A32.*2.0-A33.*2.0-A34.*2.0-A35.*2.0,A31.*2.0,A32.*2.0,A33.*2.0,A34.*2.0,A35.*2.0];
end
if nargout > 3
    duyyva = [A41.*-2.0-A42.*2.0-A43.*2.0-A44.*2.0-A45.*2.0,A41.*2.0,A42.*2.0,A43.*2.0,A44.*2.0,A45.*2.0];
end