function [Cfinal, Ctotal, Theta] = OptimizeBDRIS_MOU(Hd,F,G,Theta,Rxx,sigma2n,varargin)

% Description: This function optimizes a beyond-diagonal RIS (BDRIS) to maximize capacity
% in a MIMO link for a fixed transmit covariance. The BDRIS matrix is
% unitary and we optimize on the manifold of unitary matrices through an MO algorithm.
%
%
% Input parameters:
% H,F,G: (direct, RIS->Tx, Tx->RIS, resp.)
% Theta : initial unitary BD-RIS
% Rxx : Tx covariance matrix
% sigma2 : noise variance
% varargin: structure with the algoritm parameters
%
% I. Santamaria, UC, July 2025
%
% 11/07/2025: First version


[N,~] = size(Hd);   % Matrix F is N \times M
[~,~] = size(F);    % Matrix F is N \times M

%% Default values
opt_params = struct();
opt_params.maxiter = 100;        % Maximum number of iterations 
opt_params.threshold = 1e-10;    % To check convergence 
opt_params.mu = -1e-6;           % initial learning rate for the manifold optimization algorithm
opt_params.alpha = 1.01;         % learning rate update parameter
if nargin < 7
    error(message('TooFewInputs'));
elseif nargin == 7
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'maxiter'
                opt_params.maxiter  = param_value;
            case 'threshold'
                opt_params.threshold  = param_value;
            case 'mu'
                opt_params.mu  = param_value;
        end
    end
elseif nargin > 7
    error(message('TooManyInputs'));
end

maxiter = opt_params.maxiter;
mu = opt_params.mu;
threshold = opt_params.threshold;
alpha = opt_params.alpha;

Heq = Hd + F*Theta*G';
true = 1;
iter = 1;
compute_tangent = 1; % flag that establishes when we should compute directions again
Rxx = Rxx/sigma2n;
Ctotal = zeros(1,maxiter);
Ctotal(1) = log2(real(det(eye(N) + (Heq*Rxx*Heq'))));
while true == 1
    %iter = iter +1;
    %Q = sqrtm(Theta);  % we need the Takagi factor for the tangent plane
    if compute_tangent==1
        Gradunc = F'*((eye(N) + (Heq*Rxx*Heq'))\(Heq*Rxx*G));  % unconstrained gradient
        SkewHermitian = (Theta'*Gradunc - Gradunc'*Theta)/2;
        Cini = log2(real(det(eye(N) + (Heq*Rxx*Heq'))));
    end
    Update = expm(mu*SkewHermitian);
    Theta1new = Theta*Update;         % New unitary matrix 
    Heqnew = Hd + F*Theta1new*G';
    nFup = log2(real(det(eye(N) + (Heqnew*Rxx*Heqnew'))));
    if nFup>Cini
        iter = iter+1;
        Theta = Theta1new;
        Heq = Heqnew;
        mu = mu*alpha;      % increase mu
        compute_tangent = 1;
        Ctotal(iter) = nFup;
    else
        compute_tangent = 0;
        mu = mu/alpha;      % decrease mu
     end
    % Check convergence
    if iter >2   % 
        if abs(Ctotal(iter)-Ctotal(iter-1))<threshold || (iter>=maxiter)
            true = 0;
        end
    end
end
Ctotal = Ctotal(1:iter);
Cfinal = Ctotal(end);
