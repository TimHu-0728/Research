function E=Energy(params)
%ENERGY

g = params.DesignVars(1).Value;
h = params.DesignVars(2).Value;
m = params.DesignVars(3).Value;
v = params.DesignVars(4).Value;

E = m*g*h+1/2*m*v^2;
end