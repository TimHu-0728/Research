function y = surf_func(x)
    %[IMPORTANT]: surf_func(x) is a parameter of the model and can be changed at will to get a different geometry!!!
    y = 0.*x;
%     y = 30 - sqrt(900 - x .^ 2);            %A sphere centered at z = 30, r = 0
%     y = x .* x * 0.01 + x .^4 / 1e4;        %custom geometry test!
end

