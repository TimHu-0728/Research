function [FVertexBCD,dFVertexBCD] = equationFVertexBCD(z1,r1,in3,in4)
%EQUATIONFVERTEXBCD
%    [FVERTEXBCD,DFVERTEXBCD] = EQUATIONFVERTEXBCD(Z1,R1,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.5.
%    22-Apr-2021 14:56:29

pa_4 = in4(4,:);
pa_6 = in4(6,:);
yf1v1d1 = in3(1,:);
yf1v2d1 = in3(8,:);
yf1v2d3 = in3(10,:);
yf1v3d3 = in3(17,:);
yf1v4d1 = in3(22,:);
t2 = pi./2.0;
t3 = -t2;
t4 = pa_4+t3;
t5 = tan(t4);
FVertexBCD = [yf1v1d1,-pa_6+yf1v2d1,yf1v3d3-t5.*yf1v2d3,yf1v4d1];
if nargout > 1
    dFVertexBCD = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,-t5,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],[4,28]);
end