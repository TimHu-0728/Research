function [Br, Bz] = B(I, r, zold, co)%Calculates magnetic field for a single magnet
    %r, zold is the position in the old coordinate system. r is the same as in the new system, so is just r.
    %co = [rm, zm] are the coordinates of the coil.     I is the current through the coil.
    %Next, we convert to the "new" coordinate system, which is translated so the coil is at z = 0, r = rm

    z                  = zold - co(2);
    r                  = r + 1e-14;
    sum_of_squares     = co(1) ^ 2 + r .^ 2 + z .^ 2;                                                            %this is a group of variables that are used often
    al2                = sum_of_squares - 2 * co(1) .* r;                                                        %alpha^2
    be2                = sum_of_squares + 2 * co(1) .* r;                                                        %beta^2
    k2                 = 1 - al2 ./ be2;                                                                         %k^2
    [K, E]             = ellipke(k2);                                                                            %E(k^2) and K(k^2) elliptic integrals of 1st and 2nd kind
    C                  = 1.25663706212e-6 * I / pi;                                                              %constant C
                                                                                                                 %note that constants and functions are from Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop (on Mendeley)
    Br                 = (C .* z) ./ (2 .* al2 .* sqrt(be2) .* r) .* (sum_of_squares .* E - al2 .* K);           %B-field in r-direction
    Bz                 = C ./ (2 .* al2 .* sqrt(be2)) .* ((co(1) ^ 2 - r .^ 2 - z .^ 2) .* E + al2 .* K);        %B-field in z-direction
%     z   = zold - co(2);
%     r   = r + 1e-14;%prevents r from being zero since it will create a singulaity in cylindrical coordinates, and causes nan.
%     al2 = co(1) ^ 2 + r .^ 2 + z .^ 2 - 2 * co(1) .* r;%alpha^2
%     be2 = co(1) ^ 2 + r .^ 2 + z .^ 2 + 2 * co(1) .* r;%beta^2
%     k2  = abs(1 - al2 ./ be2);%k^2
%     [K, E] = ellipke(k2);%E(k^2) and K(k^2) elliptic integrals of 1st and 2nd kind
%     C = 1.25663706212e-6 * I / pi;%constant C
%     %note that constants and functions are from Simple Analytic Expressions for the Magnetic Field of a Circular Current Loop (on Mendeley)
%     Br = (C .* z) ./ (2 .* al2 .* sqrt(be2) .* r) .* ((co(1) ^ 2 + r .^ 2 + z .^ 2) .* E - al2 .* K);%B-field in r-direction
%     Bz = C ./ (2 .* al2 .* sqrt(be2)) .* ((co(1) ^ 2 - r .^ 2 - z .^ 2) .* E + al2 .* K);%B-field in z-direction
end