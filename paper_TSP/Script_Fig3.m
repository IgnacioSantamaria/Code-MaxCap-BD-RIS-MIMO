
% This script reproduces Fig. 2 in [A]. It considers a MIMO link assisted by a BD-RIS. 
% The BD-RIS is optimized to maximize the Frobenius norm of the equivalent MIMO channel. 
% For the optimization we use the Manifold Optimization algorithm with Phase Optimization in [A].
% It compares the Line Search and Phase optimization methods.
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
M = 5:5:50;                 % Number of BD-RIS elements  
PmaxdBm = 0;                  % Pmax (in dBm) 
Pmax = 10.^(PmaxdBm/10);      % Pmax
P = Pmax*eye(Ntx);            % Tx Cov. matrix  (fixed)

NsimMC = 50;               % Number of Monte Carlo simulations (you should use at least 25 to get a reasonable smooth curve)

%% Parameters for the MO algorithm
opt_params = struct();
opt_params.maxiter = 200;        % Maximum number of iterations
opt_params.threshold = 1e-4;     % To check convergence
opt_params.tolerance = 1e-4;     % To check convergence (bisection for LS)

maxiter = opt_params.maxiter;
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

%% Variables to store the Frobenius norms 
LSconv = zeros(1,length(M));
LSconvLR = zeros(1,length(M));
Phasesconv = zeros(1,length(M));
PhasesconvLR = zeros(1,length(M));

for mm = 1:NsimMC  
    disp(['Simulation:', num2str(mm)])
    for dd = 1:length(M)
        disp(['RIS elements:', num2str(M(dd))])
        %% Generate channels
        PosRIS_XYZ = [posicion_RIS, 3, 3];     % (x,y,z) coordinates for the RIS position
        [Hd,G,F] = ChannelsMIMO(M(dd),Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  

        %% Initialization
        Q = orth(randn(M(dd),M(dd))+1i*randn(M(dd),M(dd)));
        Usini = Q*Q.';
        Uz = orth([F' G.']);   % basis spanning F^H and G^T
        Mlr = rank(Uz);
        Qlr = orth(randn(Mlr,Mlr)+1i*randn(Mlr,Mlr));
        UsiniLR = Qlr*Qlr.';  

        %% MO with LS algorithm (Low-rank)     
        tic_ini = tic;
        [FrobNormFinalLSLR, FrobNormLSLR, UsLSLR] = MaxFrobNormLS_LowRank(Hd,F,G,UsiniLR,opt_params);
        LSconvLR(dd) = LSconvLR(dd) + toc(tic_ini);

        %% MO with LS algorithm     
        tic_ini = tic;
        [FrobNormFinalLS, FrobNormLS, UsLS] = MaxFrobNormLS(Hd,F,G,Usini,opt_params);
        LSconv(dd) = LSconv(dd) + toc(tic_ini);

        %% MO with Phase Opt. algorithm (Low-Rank)    
        tic_ini = tic;
        [FrobNormFinalPhasesLR, FrobNormPhasesLR, UsPhasesLR] = MaxFrobNormPO_LowRank(Hd,F,G,UsiniLR,opt_params);
        PhasesconvLR(dd) = PhasesconvLR(dd) + toc(tic_ini);

        %% MO with Phase Opt. algorithm (Full-Rank)     
        tic_ini = tic;
        [FrobNormFinalPhases, FrobNormPhases, UsPhases] = MaxFrobNormPO(Hd,F,G,Usini,opt_params);
        Phasesconv(dd) = Phasesconv(dd) + toc(tic_ini);
        
    end
    
end

LSconv = LSconv/NsimMC;
LSconvLR = LSconvLR/NsimMC;
Phasesconv = Phasesconv/NsimMC;
PhasesconvLR = PhasesconvLR/NsimMC;

%% Plot results
figure(1); clf
plot(M,10*log10(LSconv),'b:','MarkerSize',ms,'LineWidth',lw);
hold on
plot(M,10*log10(Phasesconv),'r','MarkerSize',ms,'LineWidth',lw);
plot(M,10*log10(LSconvLR),'b--','MarkerSize',ms,'LineWidth',lw);
plot(M,10*log10(PhasesconvLR),'r--','MarkerSize',ms,'LineWidth',lw);
hold off
xlabel('M');
ylabel('10 log_{10} (Run Time [s])');
grid on;
legend('MO Line-Search (full rank)', 'MO Opt. Phases (full rank)','MO Line-Search (low rank)','MO Opt. Phases (low rank)','Location','best')
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
