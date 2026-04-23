
% This script reproduces Fig. 5 in [A]. It considers a MIMO link assisted by a BD-RIS. 
% The BD-RIS is optimized to either minimize the MSE or maximize the rate. 
% For the otpmization we use the Manifold Optimization algorithm with Phase Optimization in [A].
% In particular, we apply the low-rank version.
% The results compare the Rate and MSE (assuming a fixed isotropic covariance matrix) 
% obtained by both solutions.
% 
%
% I. Santamaria, UC April 2026
% 
% [A] I. Santamaria, C. Beltran, E. Jorswieck, M. Soleymani, J. Gutierrez,
% "The manifold of unitary and symmetric matrices: characterization,
% Riemannian optimization and application to BD-RIS design", submitted to
% IEEE Trans. on Signal Processing.

format compact
clc; clear;
%% Parameters
Ntx = 2;                 % Number of transmit antennas
Nrx = 2;                 % Number of receive antennas
M = 4:16:128;                % Number of BD-RIS elements  
PmaxdBm = 0;                 % Pmax (in dBm) 
Pmax = 10.^(PmaxdBm/10);     % Pmax
P = Pmax*eye(Ntx);           % Tx Cov. matrix  (fixed)
NsimMC = 100;                % Number of Monte Carlo simulations (increase it to 500 for better results)

%% Parameters for the MO algorithm
opt_paramsBDRIS_MO = struct();
opt_paramsBDRIS_MO.maxiter = 100;       % Maximum number of iterations
opt_paramsBDRIS_MO.threshold = 1e-3;    % To check convergence

%% Parameters for figures
fs = 12;   % fontsize
lw = 1.5;  % linewidth
ms = 8;    % markersize
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Channel parameters
channelparams = struct;
channelparams.blocked = 0;         % Set to 1 if direct channel is blocked
channelparams.RiceRIS = 3;         % Rician factor for the channel btw RIS and BS (if Inf-> pure LoS channels)
channelparams.RiceDirect = 0;      % Rician factor for the channel btw RIS and UEs (if 0 -> Rayleigh fading)
channelparams.pl_0 = -28;          % Path loss at a reference distance (d_0)
channelparams.alpha_RIS = 2;       % Path loss exponent for the RIS links
channelparams.alpha_direct = 4;    % Path loss exponent for the direct links
channelparams.ray_fading = 0;      % Set to 1 if all channels Rayleigh

% ===== Position of the Tx/Rx/RIS (units in meters) ======
sqr_size = 50;                  % square of size sqr_size
PosTx_XYZ = [0 0 1.5];          % Position Tx
PosRx_XYZ = [sqr_size 0 1.5];   % Position Rx
posicion_RIS = 50; % We vary the x coordinate of RIS position

B = 20;   % Bandwidth MHz
NF = 0;   % Noise Factor in dBs
noiseVariancedBm = -174 + 10*log10(B*10^6) + NF;
sigma2n = 10^(noiseVariancedBm/10);       % additive noise variance

%% Variables to store  results wrt M
Rate_OptRate = zeros(size(M));  % BD-RIS that maximizes rate
MSE_OptRate = zeros(size(M));
MSE_OptMSE = zeros(size(M));    % BD-RIS that minimizes MSE
Rate_OptMSE = zeros(size(M));

for mm = 1:NsimMC  % You can alos use parfor 

    disp(['Simulation:', num2str(mm)])
    for dd = 1:length(M)
        %disp(['RIS elements:', num2str(M(dd))])
        %% Generate channels
        PosRIS_XYZ = [posicion_RIS, 3, 3];     % (x,y,z) coordinates for the RIS position
        [Hd,G,F] = ChannelsMIMO(M(dd),Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  
        
        %% Rate optimization 
        [~, ConvRate_OptRate, ThetaRate] = Optimize_LowRankBDRIS_Rate(Hd,F,G,P,sigma2n,opt_paramsBDRIS_MO);
        Heq_Rate = Hd + F*ThetaRate*G';
        Rate_OptRate(dd) =  Rate_OptRate(dd) + log2(real(det(eye(Nrx)+ (Heq_Rate*P*Heq_Rate')/sigma2n))); 
        MSE_OptRate(dd) = MSE_OptRate(dd) + real(trace(inv(eye(Nrx)+ (Heq_Rate*P*Heq_Rate')/sigma2n)));

        %% MSE optimization 
        [~, ConvMSE_OptMSE, ThetaMSE] = Optimize_LowRankBDRIS_MSE(Hd,F,G,P,sigma2n,opt_paramsBDRIS_MO);
        Heq_MSE = Hd + F*ThetaMSE*G';
        Rate_OptMSE(dd) =  Rate_OptMSE(dd) + log2(real(det(eye(Nrx)+ (Heq_MSE*P*Heq_MSE')/sigma2n))); 
        MSE_OptMSE(dd) = MSE_OptMSE(dd) + real(trace(inv(eye(Nrx)+ (Heq_MSE*P*Heq_MSE')/sigma2n)));
        
    end
end

Rate_OptRate = Rate_OptRate/NsimMC;
MSE_OptRate =  MSE_OptRate/NsimMC;
Rate_OptMSE = Rate_OptMSE/NsimMC;
MSE_OptMSE =  MSE_OptMSE/NsimMC;

%% Plot results
figure(30);clf; plot(M,Rate_OptRate,'k:o','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(M,Rate_OptMSE,'r-*','MarkerSize',ms,'LineWidth',lw);
legend('MO-PO (rate opt.)', 'MO-PO (MSE opt.)');
ylabel('Rate (b/s/Hz)');
xlabel('M');
grid on;
hold off

%%%%%%%%%%%%%%%%%%%%%%
figure(40);clf; plot(M,10*log10(MSE_OptRate),'k:o','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(M,10*log10(MSE_OptMSE),'r-*','MarkerSize',ms,'LineWidth',lw);
legend('MO-PO (Rate opt.)', 'MO-PO (MSE opt.)');
ylabel('10 log_{10}(MSE))');
xlabel('M');
grid on;
hold off
