function B=B_rect(x,y,z,a,b,mu,I)
% Function calculate the magnetic field around a rectangular wire loop with input parameters:
% Position r = [x,y,z] (m),
% Coil length a (m), width b (m), total height h (m), number of stack N, current I (A), mu (N*A^(-2))
% 
% B=B_rect(x,y,z,a,b,mu,I)
%
% INPUTS:
%
% OUTPUTS:
%
% Calculate sin(alpha) function for each segment
sina12=(a/2-y)/sqrt((x-b/2)^2+(-z)^2+(a/2-y)^2);
sina11=-(a/2+y)/sqrt((x-b/2)^2+(-z)^2+(a/2+y)^2);
sina22=(b/2-x)/sqrt((-y-a/2)^2+(-z)^2+(b/2-x)^2);
sina21=-(b/2+x)/sqrt((-y-a/2)^2+(-z)^2+(b/2+x)^2);
sina32=(a/2+y)/sqrt((-x-b/2)^2+(-z)^2+(a/2+y)^2);
sina31=-(a/2-y)/sqrt((-x-b/2)^2+(-z)^2+(a/2-y)^2);
sina42=(b/2+x)/sqrt((y-a/2)^2+(-z)^2+(b/2+x)^2);
sina41=-(b/2-x)/sqrt((y-a/2)^2+(-z)^2+(b/2-x)^2);
% Calculate B field for each segment
B1=[mu*I/4/pi/sqrt((x-b/2)^2+(-z)^2)*(sina12-sina11)*(-sin(atan2(-z,x-b/2)));
mu*I/4/pi/sqrt((x-b/2)^2+(-z)^2)*(sina12-sina11)*0;
mu*I/4/pi/sqrt((x-b/2)^2+(-z)^2)*(sina12-sina11)*(-cos(atan2(-z,x-b/2)))];
B2=[mu*I/4/pi/sqrt((-y-a/2)^2+(-z)^2)*(sina22-sina21)*0;
mu*I/4/pi/sqrt((-y-a/2)^2+(-z)^2)*(sina22-sina21)*sin(atan2(-z,-y-a/2));
mu*I/4/pi/sqrt((-y-a/2)^2+(-z)^2)*(sina22-sina21)*-cos(atan2(-z,-y-a/2))];
B3=[mu*I/4/pi/sqrt((-x-b/2)^2+(-z)^2)*(sina32-sina31)*sin(atan2(-z,-x-b/2));
mu*I/4/pi/sqrt((-x-b/2)^2+(-z)^2)*(sina32-sina31)*0;
mu*I/4/pi/sqrt((-x-b/2)^2+(-z)^2)*(sina32-sina31)*(-cos(atan2(-z,-x-b/2)))];
B4=[mu*I/4/pi/sqrt((y-a/2)^2+(-z)^2)*(sina42-sina41)*0;
mu*I/4/pi/sqrt((y-a/2)^2+(-z)^2)*(sina42-sina41)*(-sin(atan2(-z,y-a/2)));
mu*I/4/pi/sqrt((y-a/2)^2+(-z)^2)*(sina42-sina41)*(-cos(atan2(-z,y-a/2)))];
% Calculate B field for one stack of rectangular coil
B=B1+B2+B3+B4;
end