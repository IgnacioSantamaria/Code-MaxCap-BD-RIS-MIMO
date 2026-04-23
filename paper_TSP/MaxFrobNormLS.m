function [FrobNormFinal, FrobNorm, Us] = MaxFrobNormLS(Hd,F,G,Us,varargin)

% Description: This function optimizes a unitary amd symmetric matrix 
% (beyond-diagonal RIS) to maximize || Hd + F Us G^H||_F^2
% We optimize on the manifold of unitary + symmetric matrices through an MO 
% algorithm that uses a line search procedure
%
%
% Input parameters:
% Hd,F,G: (direct, RIS->Tx, Tx->RIS, resp.)
% Usa : The initial unitary+symmetric matrix
% varargin: structure with the algoritm parameters
%
% I. Santamaria, UC, Feb. 2026
%
% 11/03/2026: First version with optimal mu obtained via bisection


%[N,~] = size(Hd);   % Matrix F is N \times M
[~,M] = size(F);    % Matrix F is N \times M

%% Default values
opt_params = struct();
opt_params.maxiter = 1000;      % Maximum number of iterations
opt_params.threshold = 1e-2;    % To check convergence
opt_params.tolerance = 1e-3;    % tolerance for bisection

if nargin < 5
    error(message('TooFewInputs'));
elseif nargin == 5
    params = varargin{1};
    for arg = fieldnames(params)'
        parameter = arg{1};
        param_value = params.(parameter);
        switch parameter
            case 'maxiter'
                opt_params.maxiter  = param_value;
            case 'threshold'
                opt_params.threshold  = param_value;
            case 'tolerance'
                opt_params.tolerance  = param_value;
        end
    end
elseif nargin > 5
    error(message('TooManyInputs'));
end

maxiter = opt_params.maxiter;
threshold = opt_params.threshold;
tolerance = opt_params.tolerance;

FrobNorm = zeros(1,maxiter);
Q = sqrtm(Us);
FrobNorm(1) = norm(Hd + F*Us*G','fro')^2;

% redefine threshold as a percentage of the norm
threshold = threshold*FrobNorm(1);

true = 1;
iter = 1;

while true==1
    iter = iter+1;
    J = F'*(Hd + F*Us*G')*G;        % unconstrained gradient
    Js = (J+J.')/2;                 % symmetrize
    R = 1i*imag(Q'*Js*conj(Q)) +1i*eps*eye(M);
    % It is important to be sure R is pure imaginary and symmetric
    % the following lines are sanity checks for that.
    R = 1i*imag(R);
    R = (R+R.')/2;
    [Vr,Dr] = eig(R);
    %Draux = diag(exp(diag(Dr)));  % this is diag(e^(1i*theta1),...,e^(1i*thetaM)) like a "RIS" matrix
    Qr = Q*Vr;
    Fr = F*Qr;
    Gr = G*conj(Qr);
    
    %% find min and max mu's
    % we check the derivatives of the capacity wrt to mu
    %% find mu_min
    mu_min = 0;
    
    %% find mu_max
    mu_max = 1e-4;
    Drmax = diag(exp(diag(Dr)*mu_max));
    Hmax = Hd + Fr*Drmax*Gr';
    Derivative_max = real(trace((Fr*Dr*Drmax*Gr')'*Hmax));
    if Derivative_max>0
        truemax = 1;
        while truemax
            mu_max = mu_max*1.5;
            Drmax = diag(exp(diag(Dr)*mu_max));
            Hmax = Hd + Fr*Drmax*Gr';
            Derivative_max = real(trace((Fr*Dr*Drmax*Gr')'*Hmax));
            if Derivative_max<0
                truemax=0;
            end
        end
    end

    %% The optimal mu is in [mu_min,mu_max];
    % we find it through bisection

    true_bis = 1;
    while true_bis
        mu = (mu_min+mu_max)/2;  % mean point
        Drmu = diag(exp(diag(Dr)*mu));
        Hr = Hd + Fr*Drmu*Gr';
        Derivativer = real(trace((Fr*Dr*Drmu*Gr')'*Hr));
        if Derivativer>0
            mu_min = mu;
        else
            mu_max = mu;
        end
        if abs(mu_max-mu_min)<tolerance
            true_bis = 0;
        end
    end

    %% Apply muopt
    %muopt = mu;
    theta = exp(diag(Dr)*mu);
    Qnew = Qr*sqrtm(diag(theta));                 % new Takagi factor
    Usnew = Qnew*Qnew.';
    nFup = norm(Hd+F*Usnew*G','fro')^2;

    %if nFup>FrobNorm(iter-1) % we check if everything went well
        Us = Usnew;
        Q = Qnew;
        FrobNorm(iter) = nFup;
    %end
    % Check convergence
    if iter >2   %
        if abs(FrobNorm(iter) - FrobNorm(iter-1))<threshold || (iter>=maxiter)
            %FrobNorm(iter) - FrobNorm(iter-1)
            true = 0;
        end
    end

end
FrobNorm = FrobNorm(1:iter);
FrobNormFinal = FrobNorm(end);
