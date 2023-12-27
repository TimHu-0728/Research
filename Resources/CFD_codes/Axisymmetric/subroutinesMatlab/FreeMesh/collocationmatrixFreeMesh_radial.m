function [dPhidx,dPhidy,dPhidxx,dPhidyy,C]=collocationmatrixFreeMesh_radial(la,f1,g1)
% Implements the radial PIM collocation matrix approach (Liu 2003, Sec.
% 5.6)
% 
% la: Neighbours matrix (NxM);
% xa: Radial position of M neighbours (NxM)
% ya: Axial position of M neighbours (NxM)

% Parameters for Multiquadrics radial basis (see Liu 2003, Sec. 5.6.2, and 
% Wang 2002)
CF  = 1.42; 
qF  = 1.03; 

% Compute number of points N and neighbours M
N   = length(la(:,1));
M   = length(la(1,:));

% Preallocation of derivative & conditioning matrices
dPhidx  = zeros(N,N);
dPhidy  = zeros(N,N);
dPhidxx = zeros(N,N);
dPhidyy = zeros(N,N);
C   = zeros(N,1);

% Compute collocation matrix
for n=1:N
    % Radial and axial positions of M neighbours
    ra = f1(la(n,:));
    za = g1(la(n,:));
    
    % Local radial distance
    rl  = (ra - ra');
    zl  = (za - za');
    r2  = rl.^2 + zl.^2;
    
    % Symmetric radial MQ moment matrix (Eq. 5.119 of Liu 2003, or 9 of Wang 2002)
    RQ   = (r2 + CF^2).^qF;
    
    % Condition number
    C(n) = cond(RQ);
    
    % Local distances 2
    r2_1 = r2(1,:);
    rl_1 = rl(1,:);
    zl_1 = zl(1,:);
    
    % Derivative of MQ radial functions at point 1 (i rows, k cols)
    dRidx   = 2*qF*(r2_1 + CF^2).^(qF-1).*rl_1;
    dRidy   = 2*qF*(r2_1 + CF^2).^(qF-1).*zl_1;
    dRi2dx2 = 4*qF*(qF-1)*(r2_1 + CF^2).^(qF-2).*rl_1.^2 + 2*qF*(r2_1 + CF^2).^(qF-1);
    dRi2dy2 = 4*qF*(qF-1)*(r2_1 + CF^2).^(qF-2).*zl_1.^2 + 2*qF*(r2_1 + CF^2).^(qF-1);
    
    % Derivative matrices of shape functions
    dPhidx(n,la(n,:))  = dRidx / RQ;
    dPhidy(n,la(n,:))  = dRidy / RQ;
    dPhidxx(n,la(n,:)) = dRi2dx2 / RQ;
    dPhidyy(n,la(n,:)) = dRi2dy2 / RQ;    
end

% Sparse matrices
dPhidx  = sparse(dPhidx);
dPhidy  = sparse(dPhidy);
dPhidxx = sparse(dPhidxx);
dPhidyy = sparse(dPhidyy);














