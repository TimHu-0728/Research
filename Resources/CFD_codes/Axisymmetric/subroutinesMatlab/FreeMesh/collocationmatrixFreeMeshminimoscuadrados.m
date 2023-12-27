function [dPhidx,dPhidy,dPhidxx,dPhidyy,C]=collocationmatrixFreeMeshminimoscuadrados(la,f1,g1)
% Implements the polynomial PIM collocation matrix approach (Liu 2003, Sec.
% 5.5)
% 
% la: Neighbours matrix (NxM);
% xa: Radial position of M neighbours (NxM)
% ya: Axial position of M neighbours (NxM)

% Exponential indices for moment matrix (up to 10 neighbours)
%E   = [0,0; 1,0; 0,1; 1,1; 2,0; 0,2; 2,1; 1,2; 3,0; 0,3]'; 
%Second order aproximation
E   = [0,0; 1,0; 0,1; 1,1; 2,0; 0,2; 2,0]';
NP=legnth(E);

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
        
    % Moment matrix (Eq. 5.90 & 5.91 of Liu 2003)
    PQ      = ra.^E(1,1:M) .* za.^E(2,1:M);
    
    % Condition number
    C(n)    = cond(PQ);
    
    % Values of field variables at the M nodes in support domain
    % pt      = ra(1).^E(1,1:M) .* za(1).^E(2,1:M);
    
    % Single and twice derivation of basis function of monomials in the 
    % (r,z) coordinates
    dE       = (E(:,1:M)-1); dE(dE<0) = 0;   % 1st derivative exponents
    Ed       = E(:,1:M);                     % 1st derivative coefficents
    d2E      = (E(:,1:M)-2); d2E(d2E<0) = 0; % 2nd deirvative exponents
    Ed2      = dE(:,1:M) .* E(:,1:M);        % 2nd derivative coefficients
    
    ptdx     = Ed(1,1:M)  .* ra(1).^dE(1,:)  .* za(1).^dE(2,:);  % [0,1,0,za(1),2*ra(1),0];
    ptdy     = Ed(2,1:M)  .* ra(1).^dE(1,:)  .* za(1).^dE(2,:);  % [0,0,1,ra(1),0,2*za(1)];
    ptdxx    = Ed2(1,1:M) .* ra(1).^d2E(1,:) .* za(1).^d2E(2,:); % [0,0,0,0,2,0];
    ptdyy    = Ed2(2,1:M) .* ra(1).^d2E(1,:) .* za(1).^d2E(2,:); % [0,0,0,0,0,2];
    
    % Derivatives of matrices of PIM shape functions at central node
    dPhidx(n,la(n,:))  = ptdx/PQ;
    dPhidy(n,la(n,:))  = ptdy/PQ;
    dPhidxx(n,la(n,:)) = ptdxx/PQ;
    dPhidyy(n,la(n,:)) = ptdyy/PQ;
end

% Sparse matrices
dPhidx  = sparse(dPhidx);
dPhidy  = sparse(dPhidy);
dPhidxx = sparse(dPhidxx);
dPhidyy = sparse(dPhidyy);














