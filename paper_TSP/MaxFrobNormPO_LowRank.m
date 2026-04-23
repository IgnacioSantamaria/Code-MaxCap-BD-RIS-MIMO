function [FrobNormFinal, FrobNorm, Us] = MaxFrobNormPO_LowRank(Hd,F,G,Us,varargin)

% Description: This function optimizes a unitary amd symmetric matrix 
% (beyond-diagonal RIS) to maximize || Hd + F Us G^H||_F^2
% We optimize on the manifold of unitary + symmetric matrices through an MO 
% algorithm that optimizes phases one by one.
%
%
% Input parameters:
% Hd,F,G: (direct, RIS->Tx, Tx->RIS, resp.)
% Us : The initial unitary+symmetric matrix
% varargin: structure with the algoritm parameters
%
% I. Santamaria, UC, Feb. 2026
%
% 11/03/2026: First version with optimal mu obtained via bisection
% 13/03/26: Low-rank version (if M>Nrx+Ntx)

[N,~] = size(Hd);   % Matrix F is assumed to be N \times M (M is the number of BD-RIS elements)
[~,M] = size(F);    % Matrix F is N \times M

%% Default values
opt_params = struct();
opt_params.maxiter = 1000;      % Maximum number of iterations
opt_params.threshold = 1e-3;    % To check convergence

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
        end
    end
elseif nargin > 5
    error(message('TooManyInputs'));
end

maxiter = opt_params.maxiter;
threshold = opt_params.threshold;

FrobNorm = zeros(1,maxiter);
%Redefine F and G for low-rank version

Uz = orth([F' G.']);   % basis spanning F^H and G^T
F = F*Uz;
G = G*conj(Uz);
RankR = rank(Uz);      % This is min (M,Nrx+Ntx)
M = RankR;   % redefine M!!
Q = sqrtm(Us);  % Note Us must have correct dimensions for low-rank!!
FrobNorm(1) = norm(Hd + F*Us*G','fro')^2;
% redefine threshold as a percentage of the norm
threshold = threshold*FrobNorm(1);

true = 1;
iter = 1;

while true == 1
    iter = iter +1;
    J = F'*(Hd + F*Us*G')*G;        % unconstrained gradient
    R = 1i*imag(Q'*(J+J.')*conj(Q)/2) +1i*eps*eye(M); % R in Eq. (2). Tangent space is 1i*Q*R*Q^T. 
    % It is important to be sure R is pure imaginary and symmetric
    % the following lines are sanity checks for that.
    R = 1i*imag(R);
    R = (R+R.')/2;  
    [Vr,Dr] = eig(R);
    Draux = diag(exp(diag(Dr)));  % this is diag(e^(1i*theta1),...,e^(1i*thetaM)) like a "RIS" matrix
    Qr = Q*Vr;

    % The geodesic is Qr*(exp(Dr)^mu)*Qr^T
    %% Now we optimize the inner phases one by one 
    % Note: a standard gradient descent with line search over mu also does
    % the work
    Fr = F*Qr;
    Gr = G*conj(Qr);
    theta = ones(size(diag(Draux)));   % starting point
    for mm = 1:M  % loop to update the mth RIS element
        mindex = 1:M;
        thetam = theta;
        mindex(mm) = [];
        thetam(mm) = [];
        Thetam = diag(thetam);
        Fm = Fr(:,mindex);  % select the fixed columns
        Gm = Gr(:,mindex);
        fm = Fr(:,mm);      % select the column to update
        gm = Gr(:,mm);
        S = Hd + Fm*Thetam*Gm';    % fixed matrix (it does not depend on m)
        angleopt = -angle(gm'*S'*fm);  % optimal phase
        theta(mm) = exp(1i*angleopt);     
    end
    Us = Qr*(diag(theta))*Qr.';              % new BD-RIS
    Q = Qr*sqrtm(diag(theta));               % new Takagi factor
    Heq = Hd + F*Us*G';                   % equivalent channel
    FrobNorm(iter) = norm(Hd+F*Us*G','fro')^2;  % This is the final solution of the inner loop
   
    % Check convergence
    if iter >2   %
        if abs(FrobNorm(iter) - FrobNorm(iter-1))<threshold || (iter>=maxiter)
            true = 0;
        end
    end

end
FrobNorm = FrobNorm(1:iter);
FrobNormFinal = FrobNorm(end);

