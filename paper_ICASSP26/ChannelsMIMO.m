function [H,G,F] = ChannelsMIMO(M,Nrx,Ntx,PosTx_XYZ, PosRx_XYZ,PosRIS_XYZ,channelparams)


% Description: Generate channels for a RIS-assisted MIMO link
%
% Input parameters:
% M: Number of RIS elements
% Nrx, Ntx: number of transmit and receive antennas for each user
% PosTx_XYZ, PosRx_XYZ, PosRIS_XYZ: positions, txs, rxs and RIS
% channelparams: structure with the channel parameters
%
% Output parameters:
% H, G, F channels. H is Nrx x Ntx matrix, G is  NtxM matrix, F is Nrx x M
%
% Ignacio Santamaria, UC 2025


H = zeros(Nrx,Ntx);     % Direct MIMO link
G = zeros(Ntx,M);       % Channel from TX to BD-RIS
F = zeros(Nrx,M);       % Channels from STAR-RIS to RX's

% Extract positions
x_tx = PosTx_XYZ(1);
y_tx = PosTx_XYZ(2);
z_tx = PosTx_XYZ(3);
x_rx = PosRx_XYZ(:,1);
y_rx = PosRx_XYZ(:,2);
z_rx = PosRx_XYZ(:,3);
x_ris = PosRIS_XYZ(1);
y_ris = PosRIS_XYZ(2);
z_ris = PosRIS_XYZ(3);

d_RIS_rx  = sqrt((x_ris-x_rx).^2+(y_ris-y_rx).^2+(z_ris-z_rx).^2);
d_tx_RIS  = sqrt((x_ris-x_tx).^2+(y_ris-y_tx).^2+(z_ris-z_tx).^2);

%% Link Tx-RIS (G)
pl_tx_ris_db  = channelparams.pl_0-10*channelparams.alpha_RIS*log10(d_tx_RIS);
pl_tx_ris_eff = 10^((pl_tx_ris_db)/20);
if channelparams.ray_fading ==1
    G = pl_tx_ris_eff* ...
        (1/sqrt(2)*randn(Ntx,M)+1i*1/sqrt(2)*randn(Ntx,M));
else
    % ======================================================
    % =============== Modeling Rician Fading ===============
    % ======================================================
    % ---- Modeling the LOS links Tx-RIS-Rx ---
    phi_AoD1=2*pi*rand;
    phi_AoA1=2*pi*rand;

    a_D_r = exp(1i*2*pi*rand)*exp(1i*pi*(0:M-1)'*sin(phi_AoD1));
    a_D_t = exp(1i*2*pi*rand)*exp(1i*pi*(0:Ntx-1)'*sin(phi_AoA1));

    if channelparams.RiceRIS == Inf %pure LoS
        G = pl_tx_ris_eff*a_D_t*a_D_r';
    else
        RiceFactor = channelparams.RiceRIS;
        G = pl_tx_ris_eff* ...
            ((sqrt(RiceFactor)/sqrt(RiceFactor+1))*a_D_t*a_D_r' ...
            + (1/sqrt(RiceFactor+1))*(1/sqrt(2)*randn(Ntx,M)+1i*1/sqrt(2)*randn(Ntx,M)));
    end
end

%% Link RIS-Rx

pl_ris_rx_db  = channelparams.pl_0-10*channelparams.alpha_RIS*log10(d_RIS_rx);
pl_ris_rx_eff = 10^((pl_ris_rx_db)/20);

if channelparams.ray_fading ==1
    F = pl_ris_rx_eff* ...
        (1/sqrt(2)*randn(Nrx,M)+1i*1/sqrt(2)*randn(Nrx,M));
else
    % ======================================================
    % =============== Modeling Rician Fading ===============
    % ======================================================
    % ---- Modeling the LOS link Tx-RIS-Rx ---

    phi_AoD1=2*pi*rand;
    phi_AoA1=2*pi*rand;
    a_D_r = exp(1i*2*pi*rand)*exp(1i*pi*(0:M-1)'*sin(phi_AoD1));
    a_D_t = exp(1i*2*pi*rand)*exp(1i*pi*(0:Nrx-1)'*sin(phi_AoA1));
    if channelparams.RiceRIS == Inf %pure LoS
        F = pl_ris_rx_eff* a_D_t*a_D_r';
    else
        RiceFactor = channelparams.RiceRIS;
        F = pl_ris_rx_eff* ...
            ((sqrt(RiceFactor)/sqrt(RiceFactor+1))*a_D_t*a_D_r' ...
            +(1/sqrt(RiceFactor+1))*(1/sqrt(2)*randn(Nrx,M)+1i*1/sqrt(2)*randn(Nrx,M)));
    end
end

%% Direct Link Hd channel
d_tx_rx  = sqrt((x_tx-x_rx)^2+(y_tx-y_rx)^2);
pl_tx_rx_db  = channelparams.pl_0-10*channelparams.alpha_direct*log10(d_tx_rx);
pl_tx_rx_eff = 10^((pl_tx_rx_db)/20);

if channelparams.ray_fading ==1
    H = pl_tx_rx_eff*(1/sqrt(2)*randn(Nrx,Ntx)+...
        1i*1/sqrt(2)*randn(Nrx,Ntx));
else
    phi_AoD1=2*pi*rand;
    phi_AoA1=2*pi*rand;

    a_D_r = exp(1i*2*pi*rand)*exp(1i*pi*(0:Nrx-1)'*sin(phi_AoD1));
    a_D_t = exp(1i*2*pi*rand)*exp(1i*pi*(0:Ntx-1)'*sin(phi_AoA1));

    if channelparams.RiceDirect == Inf %pure LoS
        H = pl_tx_rx_eff*a_D_r*a_D_t';
    else
        RiceFactor = channelparams.RiceDirect;
        H = pl_tx_rx_eff* ...
            ((sqrt(RiceFactor)/sqrt(RiceFactor+1))*a_D_r*a_D_t' ...
            + (1/sqrt(RiceFactor+1))*(1/sqrt(2)*randn(Nrx,Ntx)+1i*1/sqrt(2)*randn(Nrx,Ntx)));
    end
end
if channelparams.blocked == 1
    H = 0*H;
end
