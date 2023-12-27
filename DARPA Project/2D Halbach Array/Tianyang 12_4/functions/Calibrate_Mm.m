function [Mm,B_xy] = Calibrate_Mm(B_ref,mu0,h,w,x,y,theta,r0)
%CALIBRATE_MM Summary of this function goes here
%   Detailed explanation goes here

Mm0 = 1.48/mu0;
fun = @(Mm) Find_Mm(B_ref,mu0,Mm,h,w,x,y,theta,r0);
options = optimoptions('fmincon','Algorithm','interior-point','StepTolerance',1e-20,'OptimalityTolerance',1e-20);
[Mm,fval] = fmincon(fun,Mm0,[],[],[],[],[],[],[],options);
Ke      = Mm;
[Bx,By] = calculatingB(mu0,Ke,h,w,x,y,theta,r0);            % Calculate B at [x,y]
B_xy    = sqrt(Bx.^2+By.^2);
end

function obj = Find_Mm(B_ref,mu0,Mm,h,w,x,y,theta,r0)
    Ke      = Mm;
    [Bx,By] = calculatingB(mu0,Ke,h,w,x,y,theta,r0);        % Calculate B at [x,y]
    B_xy    = sqrt(Bx.^2+By.^2);
    obj     = abs(B_xy-B_ref);
end
