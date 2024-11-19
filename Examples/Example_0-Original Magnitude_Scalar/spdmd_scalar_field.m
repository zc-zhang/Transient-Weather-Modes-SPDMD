% spdmd_scalar_field.m

% load ('XScafull000.mat')
% load ('scaledata240805b.mat')

% Extract V0 and V1 from XScafull000
   V0 = XScafull000(:, 1:end-1);
   V1 = XScafull000(:, 2:end);
[U,S,V] = svd(V0, 'econ');
% matrix Vstar
Vstar = V';

% Determine matrix UstarX1
UstarV1 = U'*V1;
% the number of snopshots
N = size (Vstar,2);   

% Optimal DMD matrix resulting from Schmid's 2010 algorithm
%Fdmd = U'*X1*V*inv(S)
Fdmd = (UstarV1)/S;
% Determine the rank of Fdmd
r = rank(Fdmd); % set the number of modes
%r = 50;

Ur = U(:,1:r);
Sr = S(1:r,1:r);
Vr = V(:,1:r);

% Eigenvalue decomposition of Fdmd 
[Ydmd, Ddmd] = eig(Fdmd); 
Edmd = diag(Ddmd);  % eignvalues of the discrete-time system % lambda

dt=1; 
omega = log(Edmd)/dt; 

% Form Vandermonde matrix
Vand = zeros(r,N);
zdmd = Edmd;

for i = 1:N

    Vand(:,i) = zdmd.^(i-1);

end

% Determine optimal vector of amplitudes xdmd 
% Objective: minimize the least-squares deviation between 
% the matrix of snapshots X0 and the linear combination of the dmd modes
% Can be formulated as:    
% minimize || G - L*diag(xdmd)*R ||_F^2 + gamma*||xdmd||_1
%============================================================
L = Ydmd;            % modes
R = Vand;            % Vandermode matrix 
G = S*Vstar; 

Phi = Ur*Ydmd;



%Option 2 quadratic program + ADMM
% Form matrix P, vector q, and scalar s
% J = b'*P*b - q'*b - b'*q + s
% x - optimization variable (i.e., the unknown vector of amplitudes)

P = (L'*L).*conj(R*R');
q = conj(diag(R*G'*L));
s = trace(G'*G);


% Cholesky factorization of P
Pl = chol(P,'lower');

% Optimal vector of amplitudes xdmd
xdmd = (Pl')\(Pl\q);     % i.e., 


% Answer 
%gammaval=100;

%% Set a set of gamma values (sparsity level)
% % Define the parameters for generating gammaval values
 gamma_grd =350; % Number of gammaval values
min_gamma = 1e-1; % Minimum gammaval
 max_gamma =750; % Maximum gammaval
% Generate gammaval values using logspace
 gammaval = logspace(log10(min_gamma), log10(max_gamma), gamma_grd);


% Number of optimization variables
n = length(q);
% Identity matrix
I = eye(n);


% Initialize output structure
%answer = struct('gamma', [], 'Nz', [], 'Jsp', [], 'Jpol', [], 'Ploss', [], 'xsp', [], 'xpol', []);

num_values = length(gammaval);

% Allocate memory for gamma-dependent output variables
answer.gamma = gammaval;
answer.Nz    = zeros(1,length(gammaval)); % number of non-zero amplitudes 
answer.Jsp   = zeros(1,length(gammaval)); % square of Frobenius norm (before polishing)
answer.Jpol  = zeros(1,length(gammaval)); % square of Frobenius norm (after polishing)
answer.Ploss = zeros(1,length(gammaval)); % optimal performance loss (after polishing)
answer.xsp   = zeros(n,length(gammaval)); % vector of amplitudes (before polishing)
answer.xpol  = zeros(n,length(gammaval)); % vector of amplitudes (after polishing)

%% Sparsity-promoting using ADMM algorithm 
% set paramters for data processing

rho = 1;
maxiter = 10000;
Max_ADMM_Iter = maxiter;
eps_abs = 1.e-6;
eps_rel = 1.e-4;




% Cholesky factorization of matrix P + (rho/2)*I
Prho = (P + (rho/2)*I);
Plow = chol(Prho,'lower');
Plow_star = Plow';


% Loop over each gammaval value and compute the corresponding answer
for i = 1:length(gammaval)
    gamma = gammaval(i);

    % Initial conditions for ADMM
    y = zeros(n, 1); % Lagrange multiplier
    z = zeros(n, 1); % Copy of x

    % Use ADMM to solve the gamma-parameterized problem
    for ADMMstep = 1:maxiter
        % x-minimization step
        u = z - (1/rho) * y;
        xnew = Plow_star \ (Plow \ (q + (rho/2) * u));
        warning('off', 'MATLAB:nearlySingularMatrix');
        
        % z-minimization step
        a = (gamma / rho) * ones(n, 1);
        v = xnew + (1/rho) * y;
        znew = ((1 - a ./ abs(v)) .* v) .* (abs(v) > a);

        % Primal and dual residuals
        res_prim = norm(xnew - znew);
        res_dual = rho * norm(znew - z);

        % Lagrange multiplier update step
        y = y + rho * (xnew - znew);

        % Stopping criteria
        eps_prim = sqrt(n) * eps_abs + eps_rel * max([norm(xnew), norm(znew)]);
        eps_dual = sqrt(n) * eps_abs + eps_rel * norm(y);

        % Check convergence
        if (res_prim < eps_prim) && (res_dual < eps_dual)
            break;
        else
            z = znew;
        end
    end

    % Record output data
    answer.xsp(:, i) = z; % vector of amplitudes
    answer.Nz(i) = nnz(answer.xsp(:, i)); % number of non-zero amplitudes i.e., card(xdmd)
   % answer.Nz(i) = nnz(z); % number of non-zero amplitudes
    answer.Jsp(i) = real(z' * P * z) - 2 * real(q' * z) + s; % Frobenius norm (before polishing)

    % Polishing of the nonzero amplitudes
    ind_zero = find(abs(z) < 1.e-12); % find indices of zero elements of z
    m = length(ind_zero); % number of zero elements
    E = I(:, ind_zero);
    E = sparse(E);

    % Form KKT system for the optimality conditions
    KKT = [P, E; E', zeros(m, m)];
    rhs = [q; zeros(m, 1)];

    % Solve KKT system
    sol = KKT \ rhs;

    % Vector of polished (optimal) amplitudes
    xpol = sol(1:n);

    % Record output data
    answer.xpol(:, i) = xpol;
    % Polished (optimal) least-squares residual
    answer.Jpol(i) = real(xpol' * P * xpol) - 2 * real(q' * xpol) + s;
    % Polished (optimal) performance loss
    answer.Ploss(i) = 100 * sqrt(answer.Jpol(i) / s);

    i
end





%%==================
% [Norm_xdmd,Index_xdmd] = sort(xdmd,'ComparisonMethod','real');
% DEv_xdmd = Edmd(Index_xdmd);   %discrete-eigenvalues 
% DMDModes_xdmd = Phi(:,Index_xdmd);

% kk=340; rr=15 
% %sort_xsp = answer.xsp(:,kk);  
% %sort of the large amplitudes rather than rank of Eigvals
% [Norm_xsp,Index_xsp] = sort(answer.xsp(:,kk),'ComparisonMethod','real');  %return the value of order/peak of real
% DEv_xsp = Edmd(Index_xsp);   %discrete-eigenvalues 
% DMDModes_xsp = Phi(:,Index_xsp);
% 
% 
% [Norm_xpol,Index_xpol] = sort(answer.xpol(:,kk),'ComparisonMethod','real');
% DEv_xpol = Edmd(Index_xpol);   %discrete-eigenvalues 
% DMDModes_xpol = Phi(:,Index_xpol);
