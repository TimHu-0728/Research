function netplot(a,b,h,N,coil_positions,Euler_angles)
% Plot all the coils in one graph
% 
% netplot(a,b,h,N,coil_positions,Euler_angles)
% 
% INPUTS: 
% a = coil side 1 length
% b = coil side 2 length
% h

% Computations
    for i=1:size(coil_positions,1)
        color = [0.5+(cos(Euler_angles(i,2)))/2 0 0.5+(sin(Euler_angles(i,2)))/2]; % Coloring law
        plotcoil(a,b,h,N,coil_positions(i,1),coil_positions(i,2),coil_positions(i,3),Euler_angles(i,1),Euler_angles(i,2),Euler_angles(i,3),color);
        hold on
    end
    hold off
end