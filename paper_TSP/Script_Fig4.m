
% This script reproduces Fig. 4 in [A]. It considers a MIMO link assisted by a BD-RIS. 
% The BD-RIS is optimized to maximize the rate. 
% For the otpmization we use the Low-Rank version of the 
% Manifold Optimization algorithm with Phase Optimization in [A].
% 
% Note: We have removed from the comparison the MM+MO algorithm (SPAWC24
% paper) since it's too slow!
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
Ntx = 4;                 % Number of transmit antennas
Nrx = 4;                 % Number of receive antennas
M = 8:16:136;                 % Number of BD-RIS elements  
PmaxdBm = 0;                  % Pmax (in dBm) 
Pmax = 10.^(PmaxdBm/10);      % Pmax
P = Pmax*eye(Ntx);            % Tx Cov. matrix  (fixed)

NsimMC = 100;            % Number of Monte Carlo simulations 

%% Optimization parameters BDRIS (MO algorithm in U)
opt_paramsBDRIS_MOU = struct();
opt_paramsBDRIS_MOU.maxiter = 1000;       % Maximum number of iterations
opt_paramsBDRIS_MOU.threshold = 1e-3;     % To check convergence
opt_paramsBDRIS_MOU.mu = 1e-2;            % initial learning rate
opt_paramsBDRIS_MOU.alpha = 1.01;         % parameter to update learning rate
opt_paramsBDRIS_MOU.tolerance = 1e-3;     % tolerance for bisection

%% Parameters for the MO algorithm in Us (with PO)
opt_paramsBDRIS_MO = struct();
opt_paramsBDRIS_MO.maxiter = 100;        % Maximum number of iterations
opt_paramsBDRIS_MO.threshold = 1e-3;     % To check convergence

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

%% Variables to store the rates (b/s/Hz)
CBDRIS_MOUp = zeros(size(M));          % BD-RIS (unitary MO + projection)
RunTime_MOUp = zeros(size(M));
CBDRIS_MOUs_PO = zeros(size(M));     % BD-RIS (unitary+ symmetric MO)
RunTime_MOUs_PO = zeros(size(M));
CBDRIS_Mao = zeros(size(M));         % BD-RIS (unitary+ symmetric closed form)
RunTime_Mao = zeros(size(M));

for mm = 1:NsimMC  
    disp(['Simulation:', num2str(mm)])
    for dd = 1:length(M)    % You can use a parfor here
        disp(['RIS elements:', num2str(M(dd))])

        %% Generate channels
        PosRIS_XYZ = [posicion_RIS, 3, 3];     % (x,y,z) coordinates for the RIS position
        [Hd,G,F] = ChannelsMIMO(M(dd),Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  

        %% Initialization (unitary+symmetric matrix)
        [Uf,Df,Vf] = svd(F);
        [Ug,Dg,Vg] = svd(G');
        Q = Vf*Ug';
        Thetaini = (Q + Q.')/2;  % symmetric projection
        [F_SVD,~,G_SVD] = svd(Thetaini);
        dummy = F_SVD'*conj(G_SVD);   % This matrix is unitary and diagonal with distinct sv's
        F_SVD = F_SVD*sqrtm(dummy);  
        Thetaini = F_SVD*F_SVD.';

        %% Max-Cap BDRIS (MO in U + projection)
        tic_ini = tic;
        [CfinalMOU, CtotalMOU, BDRIS_MOU] = OptimizeBDRIS_MOU(Hd,F,G,Thetaini,P,sigma2n,opt_paramsBDRIS_MOU);
        BDRIS_MOUp = (BDRIS_MOU+BDRIS_MOU.')/2;
        [Fsvd,Ksvd,Gsvd] = svd(BDRIS_MOUp);
        BDRIS_MOUp = Fsvd*Gsvd';
        RunTime_MOUp(dd) = RunTime_MOUp(dd) + toc(tic_ini);
        Heq_MOUp = Hd + F*BDRIS_MOUp*G';
        CBDRIS_MOUp(dd) = CBDRIS_MOUp(dd) + log2(real(det(eye(Nrx)+ (Heq_MOUp*P*Heq_MOUp')/sigma2n)));

        %% Closed-form Mao's solution
        tic_ini = tic;
        Q = F'*Hd*G;
        Thetasymm = (Q+Q.')/2;  % symmetric projection
        [F_SVD,D_SVD,G_SVD] = svd(Thetasymm);
        dummy = F_SVD'*conj(G_SVD);   % This matrix should be unitary and diagonal with distinct sv's
        F_SVD = F_SVD*sqrtm(dummy);
        ThetaMao = F_SVD*F_SVD.';
        RunTime_Mao(dd) = RunTime_Mao(dd) + toc(tic_ini);
        Heq_Mao = Hd + F*ThetaMao*G';
        CBDRIS_Mao(dd) =  CBDRIS_Mao(dd) + log2(real(det(eye(Nrx)+ (Heq_Mao*P*Heq_Mao')/sigma2n)));

        %% Max-Cap BDRIS (MO+PO)
        tic_ini = tic;
        [CfinalMOUs_new, CtotalMOUs_new, BDRIS_MOUs_new] = Optimize_LowRankBDRIS_Rate(Hd,F,G,P,sigma2n,opt_paramsBDRIS_MO);
        RunTime_MOUs_PO(dd) = RunTime_MOUs_PO(dd) + toc(tic_ini);
        Heq_MOUs_new = Hd + F*BDRIS_MOUs_new*G';
        CBDRIS_MOUs_PO(dd) =  CBDRIS_MOUs_PO(dd) + log2(real(det(eye(Nrx)+ (Heq_MOUs_new*P*Heq_MOUs_new')/sigma2n)));

    end
end

CBDRIS_MOUp = CBDRIS_MOUp/NsimMC;
RunTime_MOUp = RunTime_MOUp/NsimMC;
RunTime_Mao = RunTime_Mao/NsimMC;
CBDRIS_Mao = CBDRIS_Mao/NsimMC;
CBDRIS_MOUs_PO = CBDRIS_MOUs_PO/NsimMC;
RunTime_MOUs_PO = RunTime_MOUs_PO/NsimMC;

%% Plot results
figure(30);clf; plot(M,CBDRIS_MOUs_PO, 'r-*','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(M,CBDRIS_MOUp,'b-s','MarkerSize',ms,'LineWidth',lw);
plot(M,CBDRIS_Mao,'Color', '#7E2F8E', 'Marker', '+','MarkerSize',ms,'LineWidth',lw);

legend('MO+PO in $\mathcal{U}_s$', 'MO in $\mathcal{U}$ + Proj.',...
    'Low Cost','Location','best');
ylabel('Rate (b/s/Hz)');
grid on;
hold off
%%%%%%%%%%%%%%%%%%%%%%
figure(40);clf; plot(M,10*log10(RunTime_MOUs_PO), 'r-*','MarkerSize',ms,'LineWidth',lw);
hold on;
plot(M,10*log10(RunTime_MOUp), 'b-s','MarkerSize',ms,'LineWidth',lw);
plot(M, 10*log10(RunTime_Mao), 'Color', '#7E2F8E', 'Marker', '+','MarkerSize',ms,'LineWidth',lw);
legend('MO+PO in $\mathcal{U}_s$', 'MO in $\mathcal{U}$ + Proj.',...
    'Low Cost','Location','best');
ylabel('10 log_{10} (Run Time [s])');
grid on;
hold off

