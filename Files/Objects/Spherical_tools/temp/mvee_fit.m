function [A, c, iter] = mvee_fit(X, tol, max_iter)
%MVEE_FIT Minimal volume enclosing ellipsoid fit
%   Iterative implementation of minimal volume enclosing ellipsoid fit
%   algorithm. The algorithm uses QR decomposition.
%   [X]: oszlopai a pontok
%   [A, C, ITER] = MVEE_FIT(X, TOL, MAX_ITER) computes the minimal volume
%   enclosing ellipsoid represented by the transformation matrix A and the
%   offset C. ITER is the number of iterations. The default values for TOL
%   and MAX_ITER are 1e-8 and 1000, respectively.
%
%   The resulting matrix A gives the central form of the ellipsoid, such 
%   that (x - c).' * A * (x - c) <= 1, for all x vectors that are columns 
%   of the input matrix X. The minimal volume ellipse that fullfils this 
%   criterion is the one for which det(A) is minimal.

%   (c) AMORES Robotics, 2015.

% Default value for tolerance
if nargin < 2
    tol = 1e-8;
end

% Default value for max_iter
if nargin < 3
    max_iter = 1000;
end

N = size(X,2);          % Number of measurements
d = size(X,1);          % Dimension
n = d + 1;

% Build Q
Q = [X; ones(1, N)];
p = 1/N * ones(N, 1);

Q_tild = (Q * diag(sqrt(p))).';
[~, U]= qr(Q_tild, 0);

V = (U.') \ Q; 
kappa = dot(V, V, 1);

err = Inf;
iter = 0;

% Loop start
while (err > tol && iter < max_iter)

    iter = iter + 1;
    
    [kappa_max, j] = max(kappa);
    beta = (kappa_max - n)/(n*(kappa_max-1));
    
    p_plus = (1-beta) * p;
    p_plus(j) = p_plus(j) + beta;
    
    Q_plus = Q;
    
    Q_plus(:, j) = U \ ((U.') \ Q(:,j));
    
    % Cholesky update
    Q_tild = [U; sqrt(beta / (1-beta)) * Q(:,j).'];
    
    [~, U_plus] = qr(Q_tild, 0);
    U_plus = sqrt(1-beta) * U_plus;
    
    gamma = 1 / kappa_max;
    delta = n * (1-gamma) / (n-1);
    sigma = (1 - n*gamma)/(1-gamma);
    
    kappa = delta * (kappa - sigma*(Q_plus(:,j).' * Q).^2 / kappa_max);
    
    err = norm(p - p_plus)/norm(p);
  
    % Update for next step
    p = p_plus;
    U = U_plus;

end % Loop end
 
% Calculate offset
c = X * p;

P_tild = (X*diag(sqrt(p))).';
[~, V] = qr(P_tild, 0);

V_tild = V.' * V - c*c.';
A = 1/d * inv(V_tild.');            %#ok<MINV>

end % of mvee_fit funcion