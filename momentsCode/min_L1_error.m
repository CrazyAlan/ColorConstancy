% x = normal coefficients at a given location to be estimated, 
% A = the lighting location matrix
% b = the measurements, and
% lambda = the coefficient of the regularization term.


function x = min_L1_error(A,b)

MAX_ITERS   = 100;  % Maximum number of iterations
% MAX_ITERS   = 20;  % Maximum number of iterations
ERR_THR     = 1e-8; % For stopping criteria
GAMMA_THR   = 1e-8; % For numerical stability 
% GAMMA_THR   = 1e-4; % For numerical stability 
sigma     = 1e-6; % Source Variance

[n,m] = size(A);
gamma = ones(n,1);
x_old = 1000*ones(m,1);

for i = 1:MAX_ITERS
    
    W = diag(1./gamma);
    x = (sigma*eye(m) + A'*W*A)\(A'*W*b); % Tikh, in a wt'ed LS
    e = b - A*x; % not stacked, in this code.
    if (norm(x-x_old) < ERR_THR)
        break;
    end;
    
    x_old = x;   
    gamma = abs(e); % L1 !!!!!!!!!!!!!!!!!
    gamma = max(gamma,GAMMA_THR); 
    
end;
        %{
        [Ls, ii]=sort(b);
        plot(b(ii))
        ax = A*x;
        hold on;plot(ax(ii),'.r');hold off
        %}


return;
