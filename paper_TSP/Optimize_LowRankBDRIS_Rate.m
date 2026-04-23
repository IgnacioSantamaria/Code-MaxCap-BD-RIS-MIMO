function [Cfinal, Ctotal, Theta] = Optimize_LowRankBDRIS_Rate(Hd,F,G,Rxx,sigma2n,varargin)

% Description: This function optimizes a beyond-diagonal RIS (BDRIS) 
% to maximize the achievable rate in a MIMO link for a fixed transmit 
% covariance. The BDRIS matrix is symetric and unitary and we optimize 
% on the manifold Us through an MO algorithm.
% - It exploits Theorem 2 in [A] to optimize a rxr BD-RIS with r=Nt+Nr,
% instead of a MxM BD-RIS.
% - It applies the PO (Phase Optimization) version in which the phases are
% optmized one by one
%
% Input parameters:
% H,F,G: (direct, RIS->Tx, Tx->RIS, resp.)
% Rxx : Tx covariance matrix
% sigma2 : noise variance
% varargin: structure with the algoritm parameters
%
% Output parameters:
% Cfinal: final rate
% Ctotal: rate vs. iterations (Cfinal = Ctotal(end))
% Theta: final unitary + symmetric BD-RIS
%
% I. Santamaria, UC, March 2026
%
% [A] I. Santamaria et al "The manifold of unitary and symmetric matrices:
% characterization, Riemmanian optimization and application to BD-RIS
% design, submitted to IEEE TSP.


[Nrx,~] = size(Hd);   % Matrix Hd is Nrx \times Ntx
%[~,M] = size(F);    % Matrix F is N \times M
%% Default values
opt_params = struct();
opt_params.maxiter = 1000;      % Maximum number of iterations
opt_params.threshold = 1e-2;    % To check convergence

if nargin < 6
    error(message('TooFewInputs'));
elseif nargin == 6
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
elseif nargin > 6
    error(message('TooManyInputs'));
end

maxiter = opt_params.maxiter;
threshold = opt_params.threshold;

true = 1;
iter = 1;
Hd = Hd*sqrtm(Rxx/sigma2n);
G = sqrtm(Rxx/sigma2n)*G;
Q = orth([F' G.']);   % Initial Q (basis spanning F^H and G^T) : Takagi factor
RankR = rank(Q);      % This is r = Nrx+Ntx
Theta = Q*Q.';        % Initial Theta rank-deficient
Heq = Hd + F*Theta*G';
Ctotal = zeros(1,maxiter);
Ctotal(1) = log2(real(det(eye(Nrx) + (Heq*Heq'))));

while true == 1
    iter = iter +1;
    J = F'*((eye(Nrx) + (Heq*Heq'))\(Heq*G));  % unconstrained gradient
    R = 1i*imag(Q'*(J+J.')*conj(Q)/2) +1i*eps*eye(size(Q,2)); % This is (Ntx+Nrx x Ntx+Nrx)
    % Checks: It is important to be sure R is pure imaginary and symmetric
    R = 1i*imag(R);
    R = (R+R.')/2;      
    [Ur,Dr] = eig(R);
    Draux = diag(exp(diag(Dr)));  % this is diag(e^jtheta1,...,e^jthetaM) like a "RIS" matrix
    Qr = Q*Ur;
    Fr = F*Qr;
    Gr = G*conj(Qr);
    theta = ones(size(diag(Draux)));  % Initial phases for PO
    for mm = 1:RankR  % loop to update the mth phase
        mindex = 1:RankR;
        thetam = theta;
        mindex(mm) = [];
        thetam(mm) = [];
        Thetam = diag(thetam);
        Fm = Fr(:,mindex);  % select the fixed columns
        Gm = Gr(:,mindex);
        fm = Fr(:,mm);      % select the column to update
        gm = Gr(:,mm);
        S = Hd + Fm*Thetam*Gm'; % fixed matrix 
        rm = S*gm;
        A = eye(Nrx)+ (S*S'+ fm*fm'*(gm'*gm));
        angleopt = angle(fm'*(A\rm));  % optimal phase
        theta(mm) = exp(1i*angleopt);
     
    end
    Theta = Qr*(diag(theta))*Qr.';           % new BD-RIS
    Q = Qr*sqrtm(diag(theta));               % new Takagi factor
    Heq = Hd + F*Theta*G';                   % equivalent channel
    Ctotal(iter) = log2(real(det(eye(Nrx) + Heq*Heq')));  % This is the final solution of the inner loop
    DeltaCap =  Ctotal(iter)- Ctotal(iter-1);
    %% Check convergence
    if (DeltaCap  < threshold) || (iter==maxiter)
        true = 0;
    end
end
Ctotal = Ctotal(1:iter);
Cfinal = Ctotal(end);


