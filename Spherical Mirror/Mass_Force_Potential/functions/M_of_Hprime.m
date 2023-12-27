function  M = M_of_Hprime(H_prime,Ms,gamma)
%M_of_Hprime The integrand in PI_m, used for integration
%   URL: https://www.researchgate.net/publication/2148158_The_Surface_Topography_of_a_Magnetic_Fluid_--_a_Quantitative_Comparison_between_Experiment_and_Numerical_Simulation
M = Ms*(coth(gamma*H_prime)-1/gamma./H_prime);
end