
% This script allows you to reproduce some of the results of the paper [A]
% 
% [A] I. Santamaria, M. Soleymani, E. Jorswirck, J. Gutierrez, C. Beltran 
% "Riemannian optimization on the manifold of unitary and symmetric matrices 
% with application to BD-RIS-assisted systems", ICASSP 2026.
% 
% This program maximizes capacity in a MIMO link with the MO algorithm in [A] and 
% some alternative optimization. We assume an isotropic transmit covariance
% matrix.
%
%
% I. Santamaria, UC Jan. 2026
%

format compact
clc; clear;

%% Parameters
Ntx = 4;                    % Number of transmit antennas
Nrx = 4;                    % Number of receive antennas
M = 8:16:136;               % Number of RIS elements (size of the scattering matrix Theta)
PmaxdBm = 0;                % Pmax (in dBm) 
Pmax = 10.^(PmaxdBm/10);    % Pmax
Rxx = (Pmax/Ntx)*eye(Ntx);  % Tx Cov. matrix  (fixed)

NsimMC = 10;              % Number of Monte Carlo simulations (you should use at least 25 to get a reasonable smooth curve)
show = 0;                 % To plot intermediate convergence/rate results (only without the parfor)

%% Optimization parameters for the SPAWC'24 MO algorithm
% If used in the comparison. Note: This is very slow for large M
opt_paramsBDRIS = struct();
opt_paramsBDRIS.niterMM = 10000;       % Maximum number of iterations for MM (inner loop)
opt_paramsBDRIS.thresholdMM = 1e-5;    % To check convergence of the inner loop
opt_paramsBDRIS.muMO = 1e-2;           % initial learning rate for the manifold optimization algorithm
opt_paramsBDRIS.niterMO = 2000;        % maximum number of iterations for the manifold optimization (MO) algorithm
opt_paramsBDRIS.thresholdMO = 1e-5;    % convergence threshold for the manifold optimization algorithm

%% Optimization parameters for the MO algorithm on U(M) (unitary but not symmetric)
opt_paramsBDRIS_MOU = struct();
opt_paramsBDRIS_MOU.maxiter = 1000;       % Maximum number of iterations
opt_paramsBDRIS_MOU.threshold = 1e-3;     % To check convergence
opt_paramsBDRIS_MOU.mu = 1e-2;            % initial learning rate
opt_paramsBDRIS_MOU.alpha = 1.01;         % parameter to update learning rate
opt_paramsBDRIS_MOU.tolerance = 1e-3;     % tolerance for bisection

%% Optimization parameters for the ICASSP'26 MO algorithm on Us(M) (unitary and symmetric)
opt_paramsBDRIS_MOUs = struct();
opt_paramsBDRIS_MOUs.maxiter = 100;        % Maximum number of iterations
opt_paramsBDRIS_MOUs.threshold = 1e-3;     % To check convergence

%% Parameters for figures
fs = 12;   % fontsize
lw = 1.5;  % linewidth
ms = 8;    % markersize
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');

%% Channel parameters
channelparams = struct;
channelparams.blocked = 1;         % Set to 1 if direct channel is blocked
channelparams.RiceRIS = 3;         % Rician factor for the channel btw RIS and BS (if Inf-> pure LoS channels)
channelparams.RiceDirect = 0;      % Rician factor for the channel btw RIS and UEs (if 0 -> Rayleigh fading)
channelparams.pl_0 = -28;          % Path loss at a reference distance (d_0)
channelparams.alpha_RIS = 3;       % Path loss exponent for the RIS links
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
CBDRIS_MOUp = zeros(size(M));    % BD-RIS (unitary MO + projection onto the set of symmetric matrices)
RunTime_MOUp = zeros(size(M));
CBDRIS_MOUs = zeros(size(M));    % BD-RIS (unitary+symmetric MO, ICASSP'26)
RunTime_MOUs = zeros(size(M));
CBDRIS_Mao = zeros(size(M));     % BD-RIS (Low-cost algorithm by Fang & Mao ref [10])
RunTime_Mao = zeros(size(M));
CBDRIS_Spawc = zeros(size(M));   % BD-RIS (MO algorithm SPAWC'24 version)
RunTime_Spawc = zeros(size(M));

for mm = 1:NsimMC  
    disp(['Simulation:', num2str(mm)])
    for dd = 1:length(M)    % you can use parfor here
        disp(['BD-RIS elements:', num2str(M(dd))])

        %% Generate channels
        PosRIS_XYZ = [posicion_RIS, 3, 3];     % (x,y,z) coordinates for the RIS position
        [Hd,G,F] = ChannelsMIMO(M(dd),Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  % With this version the RiceanFactor also affects the direct links

        %% Initialization 
        % We get a unitary+symmetric matrix, to be used as initial point
        % for iterative algorithms
        [Uf,Df,Vf] = svd(F);
        [Ug,Dg,Vg] = svd(G');
        Qu = Vf*Ug';
        Thetaini = (Qu + Qu.')/2;       % symmetric projection
        [F_SVD,~,G_SVD] = svd(Thetaini);
        dummy = F_SVD'*conj(G_SVD);   
        F_SVD = F_SVD*sqrtm(dummy);   
        Thetaini = F_SVD*F_SVD.';

        %% Max-Cap BDRIS (MOU + projection)       
        tic_ini = tic;
        [CfinalMOU, CtotalMOU, BDRIS_MOU] = OptimizeBDRIS_MOU(Hd,F,G,Thetaini,Rxx,sigma2n,opt_paramsBDRIS_MOU);
        % Projection onto the set of unitary + symmetric matrices
        BDRIS_MOUp = (BDRIS_MOU+BDRIS_MOU.')/2;
        [Fsvd,Ksvd,Gsvd] = svd(BDRIS_MOUp);
        BDRIS_MOUp = Fsvd*Gsvd';
        RunTime_MOUp(dd) = RunTime_MOUp(dd) + toc(tic_ini);        
        Heq_MOUp = Hd + F*BDRIS_MOUp*G';
        CBDRIS_MOUp(dd) = CBDRIS_MOUp(dd) + log2(real(det(eye(Nrx)+ (Heq_MOUp*Rxx*Heq_MOUp')/sigma2n)));
        
        %% Closed-form suboptimal solution (it works only with an unobstructed direct link)
        % This is the method in T. Fang and Y. Mao "A low-complexity beamformin design for BD-RIS aided
        % multiuser networks," IEEE Comm. Letters, 2024.
        tic_ini = tic;
        Q = F'*Hd*G;
        Thetasymm = (Q+Q.')/2;        % symmetric projection
        [F_SVD,D_SVD,G_SVD] = svd(Thetasymm);
        dummy = F_SVD'*conj(G_SVD);   
        F_SVD = F_SVD*sqrtm(dummy);   
        ThetaMao = F_SVD*F_SVD.';
        RunTime_Mao(dd) = RunTime_Mao(dd) + toc(tic_ini);
        Heq_Mao = Hd + F*ThetaMao*G';
        CBDRIS_Mao(dd) =  CBDRIS_Mao(dd) + log2(real(det(eye(Nrx)+ (Heq_Mao*Rxx*Heq_Mao')/sigma2n)));

        %% Max-Cap BDRIS (SPAWC24 version) 
        % this algorithm is very slow for large M, so we have commented out the 
        % following lines. If you want to include this method in the comparison, you just have to
        % uncomment the following lines and the corresponding variables..
        
        % Q = sqrtm(Thetaini);
        % tic_ini = tic;
        % [Cfinal, Ctotal, BDRIS, Qt] = OptimizeBDRIS_SPAWC24(Hd,F,G,Q,Rxx,sigma2n,opt_paramsBDRIS);
        % RunTime_Spawc(dd) =  RunTime_Spawc(dd)+ toc(tic_ini);
        % Heq_spawc = Hd + F*BDRIS*G';
        % CBDRIS_Spawc(dd) = CBDRIS_Spawc(dd) + log2(real(det(eye(Nrx)+ (Heq_spawc*Rxx*Heq_spawc')/sigma2n)));

        %% Max-Cap BDRIS (MO algorithm in Us, ICASSP'26 ref[A])       
        tic_ini = tic;
        [CfinalMOUs_new, CtotalMOUs_new, BDRIS_MOUs_new] = OptimizeBDRIS_MOUs_FullRank(Hd,F,G,Thetaini,Rxx,sigma2n,opt_paramsBDRIS_MOUs);
        RunTime_MOUs(dd) = RunTime_MOUs(dd) + toc(tic_ini);
        Heq_MOUs_new = Hd + F*BDRIS_MOUs_new*G';
        CBDRIS_MOUs(dd) =  CBDRIS_MOUs(dd) + log2(real(det(eye(Nrx)+ (Heq_MOUs_new*Rxx*Heq_MOUs_new')/sigma2n)));
        
    end
end

CBDRIS_MOUp = CBDRIS_MOUp/NsimMC;
RunTime_MOUp = RunTime_MOUp/NsimMC;
RunTime_Mao = RunTime_Mao/NsimMC;
CBDRIS_Mao = CBDRIS_Mao/NsimMC;
% Uncomment if the SPAWC'24 method is included in the comparison
% CBDRIS_Spawc = CBDRIS_Spawc/NsimMC;
% RunTime_Spawc = RunTime_Spawc/NsimMC;
CBDRIS_MOUs = CBDRIS_MOUs/NsimMC;
RunTime_MOUs = RunTime_MOUs/NsimMC;

%% Plot results
figure(30);clf; plot(M,CBDRIS_MOUs, 'r-*','MarkerSize',ms,'LineWidth',lw);
hold on;
%plot(M,CBDRIS_Spawc,'k:o','MarkerSize',ms,'LineWidth',lw);
plot(M,CBDRIS_MOUs, 'r-*','MarkerSize',ms,'LineWidth',lw);
plot(M,CBDRIS_MOUp,'b-s','MarkerSize',ms,'LineWidth',lw);
plot(M,CBDRIS_Mao,'Color', '#7E2F8E', 'Marker', '+','MarkerSize',ms,'LineWidth',lw);
legend('MO in $\mathcal{U}_s$ (ICASSP26)', ...
    'MO in $\mathcal{U}$ + Proj.', 'Low Cost');
%legend('MO in $\mathcal{U}_s$ (ICASSP26)','MO in $\mathcal{U}_s$ (SPAWC24)', ...
%    'MO in $\mathcal{U}$ + Proj.', 'Low Cost');
ylabel('Rate (b/s/Hz)');
xlabel('Number of elements (M)');
grid on;
hold off

figure(40);clf; plot(M,10*log10(RunTime_MOUp), 'b-s','MarkerSize',ms,'LineWidth',lw);
hold on;
%plot(M,10*log10(RunTime_Spawc),'k:o','MarkerSize',ms,'LineWidth',lw);
plot(M,10*log10(RunTime_MOUs), 'r-*','MarkerSize',ms,'LineWidth',lw);
plot(M, 10*log10(RunTime_Mao), 'Color', '#7E2F8E', 'Marker', '+','MarkerSize',ms,'LineWidth',lw);
legend('MO in $\mathcal{U}_s$ (ICASSP26)', ...
    'MO in $\mathcal{U}$ + Proj.', 'Low Cost');
%legend('MO in $\mathcal{U}_s$ (ICASSP26)','MO in $\mathcal{U}_s$ (SPAWC24)', ...
%    'MO in $\mathcal{U}$ + Proj.', 'Low Cost');
ylabel('10 log_{10} (Run Time [s])');
xlabel('Number of elements (M)');
grid on;
hold off

