%% Synthetic test
clear;
close("all")
maxNumCompThreads(4);

%% Set parameters

% Perturbation amplitudes
perturbations = [1e-2, 1e-6, 1e-10];

% Number of independent runs
S = 6;

% Tensor Size
N = 10;  % Number of modes
d = 100;  % Size of each mode

%% Tensor rank setup

% Rank
r1 = 50;  

% Perturbation Rank
r2 = 50; 

% Ranks to test Every 5
test_ranks = 35:5:80;

tol = 1e-16; % Tolerance is small because we are specifying the maximum rank

% Mode dimensions and ranks
I = d*ones(N,1);
R1 = [1; r1*ones(N-1,1); 1];
R2 = [1; r2*ones(N-1,1); 1];

%% Run and time experiments
L = length(test_ranks);
errorsLR = zeros(L,S,3);
errorsRandLR = zeros(L,S,3);
errorsRandOrth = zeros(L,S,3);
errorsTwoSide = zeros(L,S,3);

timeLR = zeros(L,S,3);
timeRandLR = zeros(L,S,3);
timeRandOrth = zeros(L,S,3);
timeTwoSide = zeros(L,S,3);

for p=1:3
    perturbation = perturbations(p);
    % Computing the perturbed tensor
    TT1 = TTrand(I,R1);
    TT2 = TTrand(I,R2);

    norm_TT1 = TTnorm(TT1);
    norm_TT2 = TTnorm(TT2);

    TT = TTaxby(1/norm_TT1,TT1,perturbation/norm_TT2,TT2);
    norm_TT = TTnorm(TT);

    for i = 1:L

        rank = test_ranks(i);
        ranks = [1; rank*ones(N-1,1); 1];

        for s=1:S
            start = tic;
            Y_Det = TTrounding(TT, tol, rank);
            timeLR(i,s,p) = toc(start);
            errorsLR(i,s,p) = TTnorm(TTaxby(1,TT,-1,Y_Det), "OLR") / norm_TT;
        end

        for s=1:S
            start = tic;
            Y_OrthRand = TTrounding_Orthogonalize_then_Randomize(TT, ranks);
            timeRandLR(i,s,p) = toc(start);
            errorsRandLR(i,s,p) = TTnorm(TTaxby(1,TT,-1,Y_OrthRand), "OLR") / norm_TT;
        end

        for s=1:S
            start = tic;
            Y_RandOrth = TTrounding_Randomize_then_Orthogonalize(TT, ranks);
            timeRandOrth(i,s,p) = toc(start);
            errorsRandOrth(i,s,p) = TTnorm(TTaxby(1,TT,-1,Y_RandOrth), "OLR") / norm_TT;
        end

        for s=1:S
            start = tic;
            Y_TwoSide = TTrounding_Two_Sided_Randomization(TT, ranks);
            timeTwoSide(i,s,p) = toc(start);
            errorsTwoSide(i,s,p) = TTnorm(TTaxby(1,TT,-1,Y_TwoSide), "OLR") / norm_TT;
        end

    end
end


%% Plot Results

% Post-process errors

[errLR,       negLR,       posLR]       = computeError(errorsLR);
[errRandLR,   negRandLR,   posRandLR]   = computeError(errorsRandLR);
[errRandOrth, negRandOrth, posRandOrth] = computeError(errorsRandOrth);
[errTwoSide,  negTwoSide,  posTwoSide]  = computeError(errorsTwoSide);

f = figure(1);
f.Position(1:2) = [0,1050];
f.Position(3:4) = [1050, 700];

for p=1:3
    subplot(2,2,p)

    errorbar(test_ranks, errLR(:,p), negLR(:,p), posLR(:,p),     'bo','markersize',10,'linewidth',2)
    ax = gca;
    ax.YScale = 'log';

    hold on
    errorbar(test_ranks, errRandLR(:,p),   negRandLR(:,p),   posRandLR(:,p),    'ks','markersize',10,'linewidth',2)
    errorbar(test_ranks, errRandOrth(:,p), negRandOrth(:,p), posRandOrth(:,p),  'rx','markersize',10,'linewidth',2)
    errorbar(test_ranks, errTwoSide(:,p),  negTwoSide(:,p),  posTwoSide(:,p),   '+' ,'markersize',10,'linewidth',2)
    hold off
    title(['\epsilon = ',num2str(perturbations(p))])
    xlabel('Maximum Target Rank', 'FontSize', 18)
    ylabel('Relative Error', 'FontSize', 18)
    legend('TT-Rounding','Orth-then-Rand', 'Rand-then-Orth','Two-Sided-Rand')
    set(gca,'FontSize',16)
end

% Post-process timings

[speedupLR,         negLR,          posLR]          = computeSpeedup(timeLR,       timeLR);
[speedupRandLR,     negRandLR,      posRandLR]      = computeSpeedup(timeRandLR,   timeLR);
[speedupRandOrth,   negRandOrth,    posRandOrth]    = computeSpeedup(timeRandOrth, timeLR);
[speedupTwoSide,    negTwoSide,     posTwoSide]     = computeSpeedup(timeTwoSide,  timeLR);

subplot(2,2,4)

errorbar(test_ranks, speedupLR,       negLR,       posLR,       'bo-','markersize',10,'linewidth',2);
hold on
errorbar(test_ranks, speedupRandLR,   negRandLR,   posRandLR,   'ks-','markersize',10,'linewidth',2);
errorbar(test_ranks, speedupRandOrth, negRandOrth, posRandOrth, 'rx-','markersize',10,'linewidth',2);
errorbar(test_ranks, speedupTwoSide,  negTwoSide,  posTwoSide,  '+-' ,'markersize',10,'linewidth',2);
    
hold off
title("Observed speedup relative to TT-rounding")
xlabel('Maximum Target Rank', 'FontSize', 18)
ylabel('Speedup', 'FontSize', 18)
set(gca,'FontSize',16)


exportgraphics(f,'figures/synth.eps')


%% Custom function to post-process errors

function [error, neg, pos] = computeError(errors)
    errors = errors(:,2:end,:);
    error = squeeze(median(errors, 2));
    neg = error -  squeeze(min(errors,[],2));
    pos = squeeze(max(errors,[],2)) - error;
end

%% Custom function to post-process runtimes

function [speedup, neg, pos] = computeSpeedup(time, timeref)
    time = time(:,2:end,:);
    timeref = timeref(:,2:end,:);
    speedup = median(timeref, [2,3]) ./ median(time, [2,3]);
    neg = speedup - median(timeref,[2,3]) ./ max(time,[],[2,3]);
    pos = median(timeref,[2,3]) ./ min(time,[],[2,3]) - speedup;
end