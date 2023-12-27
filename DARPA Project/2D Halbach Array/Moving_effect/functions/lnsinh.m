function ls_out = lnsinh(input)
    %This function calculates ls_out = ln(sinh(input)/input) without floating-point errors. It is used in an integral in massForcePotential
    ls_out = abs(input) - log(abs(input)) + log(1 - exp(-2 * abs(input))) - log(2);
end

