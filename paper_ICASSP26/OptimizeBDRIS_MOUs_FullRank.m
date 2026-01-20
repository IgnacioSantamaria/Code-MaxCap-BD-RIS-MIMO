function [Cfinal, Ctotal, Theta] = OptimizeBDRIS_MOUs_FullRank(Hd,F,G,Theta,Rxx,sigma2n,varargin)

% Description: 
% 
% This function optimizes a beyond-diagonal RIS (BD-RIS) to maximize capacity
% in a MIMO link for a fixed transmit covariance. The BDRIS matrix is
% unitary and we optimize on the manifold of unitary + symmetric matrices 
% through an MO algorithm proposed in: 
%
% I. Santamaria, M. Soleymani, E. Jorswierck, J. Gutierrez, C. Beltran, "Riemannian optimization on the manifold of unitary and symmetric
% matrices with application to BD-RIS-assisted systems," ICASSP 2026,
% Barcelona, May 2026.
%
% Input parameters:
% H,F,G: (direct, RIS->Tx, Tx->RIS, resp.)
% Theta : The initial unitary BD-RIS
% Rxx : Tx covariance matrix
% sigma2 : noise variance
% varargin: structure with the algoritm parameters
%
% Output parameters:
% Cfinal: final rate 
% Ctotal: vector with rate vs iterations
% Theta: Final symmetric and unitary BD-RIS
%
% I. Santamaria, UC, July 2025
%
% 16/07/2025: Exploiting the geodesic parametrization, the MO algorithm optimizes the phases one by one


[N,~] = size(Hd);   % Matrix F is assumed to be N \times M (M is the number of BD-RIS elements)
[~,M] = size(F);    % Matrix F is N \times M

%% Default values
opt_params = struct();
opt_params.maxiter = 1000;      % Maximum number of iterations
opt_params.threshold = 1e-3;    % To check convergence

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
        end
    end
elseif nargin > 7
    error(message('TooManyInputs'));
end

maxiter = opt_params.maxiter;
threshold = opt_params.threshold;

Q = sqrtm(Theta);  % starting with Theta in Us, we have Theta = Q*Q^T
true = 1;
iter = 1;
Hd = Hd*sqrtm(Rxx/sigma2n);
G = sqrtm(Rxx/sigma2n)*G;
Heq = Hd + F*Theta*G';
Ctotal = zeros(1,maxiter);
Ctotal(1) = log2(real(det(eye(N) + (Heq*Heq'))));
while true == 1
    iter = iter +1;
    J = F'*((eye(N) + (Heq*Heq'))\(Heq*G));           % unconstrained gradient
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
        rm = S*gm;
        A = eye(N)+ (S*S'+ fm*fm'*(gm'*gm));
        angleopt = angle(fm'*(A\rm));  % optimal phase
        theta(mm) = exp(1i*angleopt);
     
    end
    Theta = Qr*(diag(theta))*Qr.';           % new BD-RIS
    Q = Qr*sqrtm(diag(theta));               % new Takagi factor
    Heq = Hd + F*Theta*G';                   % equivalent channel
    Ctotal(iter) = log2(real(det(eye(N) + Heq*Heq')));  % This is the final solution of the inner loop
    DeltaCap =  Ctotal(iter)- Ctotal(iter-1);
    %% Check convergence
    if (DeltaCap  < threshold) || (iter==maxiter)
        true = 0;
    end
end
Ctotal = Ctotal(1:iter);
Cfinal = Ctotal(end);


