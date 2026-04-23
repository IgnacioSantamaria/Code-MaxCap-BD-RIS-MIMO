function [MSEfinal, MSEtotal, Theta] = Optimize_LowRankBDRIS_MSE(Hd,F,G,Rxx,sigma2n,varargin)

% Description: This function optimizes a beyond-diagonal RIS (BDRIS) to maximize capacity
% in a MIMO link for a fixed transmit covariance. The BDRIS matrix is
% unitary and we optimize on the manifold of unitary + symmetric matrices through an MO algorithm.
% The critrion is the MSE
%
% Input parameters:
% H,F,G: (direct, RIS->Tx, Tx->RIS, resp.)
% Rxx : Tx covariance matrix
% sigma2 : noise variance
% varargin: structure with the algoritm parameters
%
% I. Santamaria, UC, July 2025
%
% 16/07/2025: Optimizes the phases one by one


[Nrx,Ntx] = size(Hd);   % Matrix F is N \times M
N  = min(Nrx,Ntx);
[~,M] = size(F);    % Matrix F is N \times M
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
Q = orth([F' G.']);   % Initial Q (basis spanning F^H and G^T) : Takagi factor
RankR = rank(Q);      % This is Nrx+Ntx
Theta = Q*Q.';        % Initial Theta rank-deficient

true = 1;
iter = 1;
Hd = Hd*sqrtm(Rxx/sigma2n);
G = sqrtm(Rxx/sigma2n)*G;
Heq = Hd + F*Theta*G';
MSEtotal = zeros(1,maxiter);
E = eye(N) + Heq*Heq'; 
MSEtotal(1) = real(trace(inv(E)));
global A
global b
global c
while true == 1
    iter = iter +1;
    %Q = sqrtm(Theta);  % we need the Takagi factor for the tangent plane
    Ei = inv(eye(N) + Heq*Heq'); 
    J = F'*(Ei^2)*(Heq*G);  % unconstrained gradient (mse)
    R = 1i*imag(Q'*(J+J.')*conj(Q)/2) + 1i*eps*eye(size(Q,2));
    % Checks: It is important to be sure R is pure imaginary and symmetric
    R = 1i*imag(R);
    R = (R+R.')/2;  
    [Ur,Dr] = eig(R);
    Draux = diag(exp(diag(Dr)));  % this is diag(e^jtheta1,...,e^jthetaM) like a "RIS" matrix
    Qr = Q*Ur;
    Fr = F*Qr;
    Gr = G*conj(Qr);

    theta = ones(size(diag(Draux)));  % diag(Draux);   % This is the actual "RIS" matrix
    for mm = 1:RankR  % loop to update the mth RIS element
        mindex = 1:RankR;
        thetam = theta;
        mindex(mm) = [];
        thetam(mm) = [];
        Thetam = diag(thetam);
        Fm = Fr(:,mindex);  % select the fixed columns
        Gm = Gr(:,mindex);
        fm = Fr(:,mm);      % select the column to update
        gm = Gr(:,mm);
        S = Hd + Fm*Thetam*Gm';    % fixed matrix (it does not depend on m)
 
        A = eye(N)+ (S*S'+ fm*fm'*(gm'*gm));
        b = fm;
        c = S*gm;

        %% optimal phase
        myfun = @(x,A,b,c)  imag(exp(1i*x)*c'*(inv(A + exp(1i*x)*b*c' + exp(-1i*x)*c*b')^2)*b)^2;  % parameterized function
        fun = @(x) myfun(x,A,b,c);    % function of x alone
        angleopt = fminbnd(fun, -pi, pi);
        theta(mm) = exp(1i*angleopt);
     
    end
    Theta = Qr*(diag(theta))*Qr.';           % new BD-RIS
    Q = Qr*sqrtm(diag(theta));               % new Takagi factor
    Heq = Hd + F*Theta*G';                   % equivalent channel
    E = eye(N) + Heq*Heq'; 
    MSEtotal(iter) = trace(inv(E));  % This is the final solution of the inner loop
    Deltamse =  MSEtotal(iter-1) - MSEtotal(iter);
    %% Check convergence
    if (Deltamse  < threshold) || (iter==maxiter)
        true = 0;
    end
end
MSEtotal = MSEtotal(1:iter);
MSEfinal = MSEtotal(end);


