function [e, Ms, xi, keps] = find_E(lambda, hratio,RhovsMs_cfit)
    N      = 1000;
    Ms     = 28500;%initial guess: EMG700 ferrofluid
    
    dMs    = 100;
    [g, e, xi, keps]   = find_GandE(lambda, hratio, Ms,RhovsMs_cfit);
    i = 0;
    while (abs(g + 10) > 1e-9 & i < 500);
        [g_high, e, xi, keps] = find_GandE(lambda, hratio, Ms + dMs,RhovsMs_cfit);
        [g_low, e, xi, keps] = find_GandE(lambda, hratio, Ms - dMs,RhovsMs_cfit);
        dg = (g_high - g_low) / (2 * dMs);
        Ms = Ms - (g + 10) / dg;
        [g, e, xi, keps]   = find_GandE(lambda, hratio, Ms,RhovsMs_cfit);
        i = i+1;
    end
end