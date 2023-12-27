function phid_c = exactphiBCs(r,z,Ib,RminD, RmaxD, ZminD, ZmaxD)
% This function computes the exact phi potential at the points r,z for a
% given Ib parameter (current intensity times number of turns) and the
% shape of the coil
%
% V1.0 (18/04/2021)

% Parameters
magnet   = 0;      % Are we simulating a magnet? Override configuration with magnet = 1
nritz    = 20;     % Number of coil discretization units per side of magnet
M0       = 1500e3; % Magnet Magnetization (A/m)

% Numerical correction
r(r<0)   = 0;

% Surface discretization of the coil
phid_c    = 0;
if magnet
    RmaxD     = 0.0275;
    RminD     = 0.0025;
    ZmaxD     = -0.001;
    ZminD     = -0.004;
    Ib        = M0 * (ZmaxD-ZminD);
    for j = 1:nritz
        za     = ZminD + (ZmaxD-ZminD)/nritz * (1/2 + j-1);
        k_i      = sqrt(4 * RminD * r ./ ((RminD+r).^2 + (z-za).^2));
        [K_i,E_i]  = ellipke(k_i.^2);
        k_e      = sqrt(4 * RmaxD * r ./ ((RmaxD+r).^2 + (z-za).^2));
        [K_e,E_e]  = ellipke(k_e.^2);
        phid_c = phid_c - Ib*sqrt(RminD.*r)./(pi*k_i*nritz) .* ((1-k_i.^2/2).*K_i - E_i)...
            + Ib*sqrt(RmaxD.*r)./(pi*k_e*nritz) .* ((1-k_e.^2/2).*K_e - E_e);
    end
else
    for i = 1:nritz
        for j = 1:nritz
            ra     = RminD + (RmaxD-RminD)/nritz * (1/2 + i-1);
            za     = ZminD + (ZmaxD-ZminD)/nritz * (1/2 + j-1);
            k      = sqrt(4 * ra * r ./ ((ra+r).^2 + (z-za).^2));
            [K,E]  = ellipke(k.^2);
            phid_c = phid_c + Ib*sqrt(ra.*r)./(pi*k*nritz^2) .* ((1-k.^2/2).*K - E);
        end
    end
%     for i = 1:nritz
%         for j = 1:nritz
%             ra     = RminD + (RmaxD-RminD)/nritz * (1/2 + i-1);
%             za     = ZminD + (ZmaxD-ZminD)/nritz * (1/2 + j-1)  -  0.368;
%             k      = sqrt(4 * ra * r ./ ((ra+r).^2 + (z-za).^2));
%             [K,E]  = ellipke(k.^2);
%             phid_c = phid_c + Ib*sqrt(ra.*r)./(pi*k*nritz^2) .* ((1-k.^2/2).*K - E);
%         end
%     end
end

phid_c;