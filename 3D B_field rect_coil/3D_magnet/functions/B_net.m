function [Bx_net,By_net,Bz_net]=B_net(x,y,z,a,b,N,I,h,mu,coil_positions,Euler_angles)
% Function calculates the net magnetic field B_net with coils centered at different positions, and different orientations in Euler angles.
% 
% [Bx_net,By_net,Bz_net]=B_net(x,y,z,a,b,N,I,h,mu,coil_positions,Euler_angles)
%
% INPUTS:
%
% OUTPUTS:
%

% Computations
    B=zeros(size(coil_positions,1),3);
% Calculate the magnetic field components Bx, By, Bz of each coil with shifted coil origin
    for i=1:size(coil_positions,1)
        B(i,:)=getB(x-coil_positions(i,1),y-coil_positions(i,2),z-coil_positions(i,3),a,b,N,I,h,mu,Euler_angles(i,1),Euler_angles(i,2),Euler_angles(i,3))';
    end
% Merge together to form the net magnetic field
    Bnet=sum(B);
    Bx_net=Bnet(1);
    By_net=Bnet(2);
    Bz_net=Bnet(3);
end
