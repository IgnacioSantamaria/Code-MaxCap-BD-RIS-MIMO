
% This script reproduces Fig. 2 in [A]. It considers a MIMO link assisted by a BD-RIS. 
% The BD-RIS is optimized to maximize the rate. 
% The Script only illustrates the convergence of the MO+PO algorithm.
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
M = 64;                  % Number of BD-RIS elements  
PmaxdBm = 0;                   % Pmax (in dBm) 
Pmax = 10.^(PmaxdBm/10);       % Pmax
P = Pmax*eye(Ntx);             % Tx Cov. matrix  (fixed)

NsimMC = 100;             % Number of Monte Carlo simulations 

%% Parameters for the MO algorithm in Us (with PO)
opt_paramsBDRIS_MOUnew = struct();
opt_paramsBDRIS_MOUnew.maxiter = 100;       % Maximum number of iterations
opt_paramsBDRIS_MOUnew.threshold = 1e-3;     % To check convergence

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

%% Generate channels  (channel realization (fixed))
ConvRate = cell(1,NsimMC);
PosRIS_XYZ = [posicion_RIS, 3, 3];     % (x,y,z) coordinates for the RIS position
[Hd,G,F] = ChannelsMIMO(M,Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams);  % With this version the RiceanFactor also affects the direct links

for mm = 1:NsimMC    
    disp(['Simulation:', num2str(mm)])
    %% Initialization (unitary+symmetric matrix)
    Q = orth(randn(M,M) + 1i*rand(M,M));
    Thetaini = Q*Q.';
   %% Max-Cap BDRIS (MO+PO, Full-Rank version (ICASSP2026))
    tic_ini = tic;
    [CfinalMOUs_new, CtotalMOUs_new, BDRIS_MOUs_new] = Optimize_FullRankBDRIS_Rate(Hd,F,G,Thetaini,P,sigma2n,opt_paramsBDRIS_MOUnew);
    Heq_MOUs_new = Hd + F*BDRIS_MOUs_new*G';
    ConvRate{1,mm} = CtotalMOUs_new;  
end

% Assuming ConvRate is a 1xNsim cell array
Nsim = numel(ConvRate);

% Determine max length of all curves
maxLen = max(cellfun(@numel, ConvRate));

% Preallocate with NaNs for padding
allCurves = NaN(Nsim, maxLen);

for i = 1:Nsim
    curve = ConvRate{i};
    allCurves(i, 1:numel(curve)) = curve;
end

% Compute mean across simulations, ignoring NaNs
avgCurve = nanmean(allCurves, 1);

% Compute min and max curves at each iteration
minCurve = nanmin(allCurves, [], 1);
maxCurve = nanmax(allCurves, [], 1);

% Plotting
figure(1);clf;
hold on;
% Fill area between min and max curves
fill([0:maxLen-1, fliplr(0:maxLen-1)], [maxCurve, fliplr(minCurve)], ...
     [0.8, 0.9, 1], 'EdgeColor', 'none');  % Light blue shaded area

% Plot average convergence curve
plot(0:maxLen-1, avgCurve, 'b', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Rate');
grid on
legend({'Range', 'Average Curve'}, 'Location', 'Best');
hold off;
