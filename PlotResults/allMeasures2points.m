
%%  Extract simulation reslults from *.mat file pfoduces by run_simulation.m
sim1 = load('Alpha0.1/LamSimK10-C3-M3-carRate0.1til25-loops50000000-samples50000000.mat');  
maxLam = sim1.lamSim(end);

%%  Tailor expansion is calculated around poins given in LamJac array
LamJac = [1, 20];  nPlots = length(LamJac);
LamJac1 = '1'; LamJac2 = '20'; % string values are specified for legend
PLOT_TERMS = [5, 15, 40]; % number of terms included in different plots

%%  Fetch terms of series expansions from *.dat files produced by calculateTermsLam.cpp
data = []; 
data(1).x = dlmread('Alpha0.1/LamTermsK10-M3-C3-lam1.0-carRate0.1-W0.9-stop1e-05.dat');
data(2).x = dlmread('Alpha0.1/LamTermsK10-M3-C3-lam20.0-carRate0.1-W0.9-stop1e-05.dat');

Terms = []; 
for i = 1:nPlots
    Terms(i).Q     = data(i).x(:,1); 
    Terms(i).Acc = data(i).x(:,2); 
    Terms(i).Tr = data(i).x(:,3); 
end
nTerms = min(length(Terms(1).Q), length(Terms(2).Q));

%%  Mean number of cars in the system (calculated in script derive_average_number_of_cars.m)
%EG = 3.6683; % Optimistic; alpha = 0.1, K = 10
%EG = 6.6714; % Optimistic; alpha = 0.2, K = 10
%EG = 3.5057; % Optimistic; alpha = 0.2, K=7
%EG = 8.1549; % Optimistic; alpha = 0.3, K = 10
EG = 2.3996; % Pessimistic; alpha = 0.1; K = 10

%% fetch simulation results
sim1.MEASURES.Tr = sim1.MEASURES.Tr/EG;
sim1.MEASURES.Acc = sim1.MEASURES.Acc/EG;
beta = [1/30,1/15,1/30];
R = sum(1./beta);	% Residence time
R2 =  sum(1./beta.^2)+(sum(1./beta)).^2
ResTime = R2/R/2;

nPlotSamples = 300;
PLOT_MODE = 'LOG'; % can take value LOG or LIN

minLamLog = -1; minLam = 10^minLamLog;
diap = [[minLamLog,log10(maxLam)];[minLamLog,log10(maxLam)];[minLamLog,log10(maxLam)]];
lam = [logspace(diap(1,1), diap(1,2), nPlotSamples)',...
             logspace(diap(2,1), diap(2,2), nPlotSamples)',...
             logspace(diap(3,1), diap(3,2), nPlotSamples)'];
RES = [];
for iPlot = 1:nPlots
    RES(iPlot).Q       = zeros(nPlotSamples,nTerms);
    RES(iPlot).Acc   = zeros(nPlotSamples,nTerms);
    RES(iPlot).Tr      =  zeros(nPlotSamples,nTerms);
end    
for i = 1:nPlotSamples
    for iPlot = 1:nPlots
        %% 0-term caclucation
        RES(iPlot).Q(i,1)      = Terms(iPlot).Q(1);
        RES(iPlot).Acc(i, 1) = lam(i,iPlot)*Terms(iPlot).Acc(1)/EG;
        RES(iPlot).Tr(i, 1)    = Terms(iPlot).Tr(1)/EG;
        for iTerm = 2:nTerms
            powTerm = (lam(i,iPlot)-LamJac(iPlot))^(iTerm-1);
            RES(iPlot).Q(i,iTerm)    = RES(iPlot).Q(i,iTerm-1)      + Terms(iPlot).Q(iTerm)*powTerm;
            RES(iPlot).Acc(i,iTerm) = RES(iPlot).Acc(i,iTerm-1)  + lam(i,iPlot)/EG*Terms(iPlot).Acc(iTerm)*powTerm;
            RES(iPlot).Tr(i,iTerm)    = RES(iPlot).Tr(i,iTerm-1)     + Terms(iPlot).Tr(iTerm)*powTerm/EG;
        end
    end
end
%%  Delay time of all packets and successfully transmitted packets
%% + Dicrarded & Rejected packets rate calculation
for iPlot = 1:nPlots
    RES(iPlot).ED = RES(iPlot).Q./RES(iPlot).Acc/EG;
    RES(iPlot).EDs = RES(iPlot).Q./EG./RES(iPlot).Tr - (RES(iPlot).Acc./RES(iPlot).Tr-1)*ResTime;
    RES(iPlot).DIS = RES(iPlot).Acc-RES(iPlot).Tr;
    RES(iPlot).REJ = lam(:,iPlot)*ones(1,nTerms)-RES(iPlot).Acc;
end
sim1.MEASURES.ED = sim1.MEASURES.EQ/EG./sim1.MEASURES.Acc;
sim1.MEASURES.DIS = sim1.MEASURES.Acc-sim1.MEASURES.Tr;

%% Plot parameters
FIG_POSISION = [32 3 12 22];
LineStyle1 = '--k'; LineStyle2 = '-k';
fontSize = 12;
markerSize = 13; simMarker = 5;
lw = 0.8;   % LineWidth
yPl = length(PLOT_TERMS); % Number of subplots in every plot  

%% Plot mean queue content EQ
fig = figure; set(fig, 'units', 'centimeters', 'Position', FIG_POSISION) % set up figure
titlePos = [0.5,2.5]; % position of the title showing number of terms
maxEQ = 3; % limit on y
DIAPdef = 1:nPlotSamples; % default diapasone includes all samples in the plot

% EQ PLOT1
DIAP = {DIAPdef,DIAPdef,DIAPdef}; 
subplot(yPl, 1, 1); plTerm = PLOT_TERMS(1); plot_EQ
lgnd = legend('sim', strcat(['\lambda_0 = ', LamJac1]),...
                                      strcat(['\lambda_1 = ', LamJac2])) 
% EQ PLOT 2
DIAP = {DIAPdef,DIAPdef,DIAPdef}; 
subplot(yPl, 1, 2); plTerm = PLOT_TERMS(2); plot_EQ
% EQ PLOT 3
DIAP = {DIAPdef,DIAPdef,DIAPdef}; 
subplot(yPl, 1, 3); plTerm = PLOT_TERMS(3); plot_EQ

%% Plot Delay 
fig = figure; set(fig, 'units', 'centimeters', 'Position', FIG_POSISION)
maxDelay = 2; titlePos = [0.5,maxDelay*0.9];
DIAP = {DIAPdef,DIAPdef,DIAPdef}; 

% Delay PLOT 1
subplot(3,1, 1); plTerm = PLOT_TERMS(1); plot_delay  % run script plotting delay
lgnd = legend('sim', strcat(['\lambda_0 = ', LamJac1]),...
                                      strcat(['\lambda_1 = ', LamJac2])) 
% Delay PLOT 2
subplot(3,1,2); plTerm = PLOT_TERMS(2); plot_delay
% Delay PLOT 3
subplot(3,1,3); plTerm = PLOT_TERMS(3); plot_delay

%% Plot rate of discarted packets
fig = figure; set(fig, 'units', 'centimeters', 'Position', FIG_POSISION)
maxDiscard = 0.15;
titlePos = [0.5,maxDiscard*0.9];
DIAP = {DIAPdef,DIAPdef,DIAPdef}; 

% Discard PLOT1
subplot(yPl,1,1); plTerm = PLOT_TERMS(1); plot_discarded
lgnd = legend('sim', strcat(['\lambda_0 = ', LamJac1]),...
                                      strcat(['\lambda_1 = ', LamJac2])) 
% Discard PLOT2
subplot(yPl,1,2); plTerm = PLOT_TERMS(2); plot_discarded
% Discard PLOT3
subplot(yPl,1,3); plTerm = PLOT_TERMS(3); plot_discarded

%% Plot rate of rejected packets
fig = figure ; set(fig, 'units', 'centimeters', 'Position', FIG_POSISION)
maxReject = 40; 
titlePos = [0.5,maxReject*0.9];
DIAP = {DIAPdef,DIAPdef,DIAPdef}; 

% Rejected PLOT 1
subplot(yPl,1,1); plTerm = PLOT_TERMS(1); plot_rejected
lgnd = legend('sim', strcat(['\lambda_0 = ', LamJac1]),...
                                      strcat(['\lambda_1 = ', LamJac2])) 
% Rejected PLOT 2
subplot(yPl,1,2); plTerm = PLOT_TERMS(2); plot_rejected;
% Rejected PLOT 3
subplot(yPl,1,3); plTerm = PLOT_TERMS(3); plot_rejected

