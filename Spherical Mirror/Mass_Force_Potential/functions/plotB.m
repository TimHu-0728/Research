function plotB(z1, r2, z2, N, Nq, I,co)
%plots (for now) - q is for quiver
    r = linspace(0, r2, N); z = linspace(z1, z2, N);
    pre_rq = linspace(0, r2, Nq + 1); pre_zq = linspace(z1, z2, Nq + 1);%must be changed to center things
    rq = (pre_rq(1:Nq) + pre_rq(2:Nq + 1)) / 2; zq = (pre_zq(1:Nq) + pre_zq(2:Nq + 1)) / 2;
    [R, Z] = meshgrid(r, z);
    [U, V] = multiB(R, Z, I,co);
    [Rq, Zq] = meshgrid(rq, zq);
    [Uq, Vq] = multiB(Rq, Zq, I,co);
    hold on
    plt3(R, Z, U, V, 0)
    plt3(Rq, Zq, Uq, Vq, 1)
    hold off
end