%% Hilbert Tensor
clear

%% Tensor Size

% Number of modes
N = 10;  

% Modal dimensions
min_mode = 6;
mode_increment = 2;
max_mode = (N-1)*mode_increment + min_mode;
I = (min_mode:mode_increment:max_mode)';

%% Parameters and tensor rank setup

% number of independent runs
S = 6;

min_rank = 8;
rank_increment = 1;
max_rank = (N-2)*rank_increment + min_rank;

true_ranks = [1; (min_rank:rank_increment:max_rank)'; 1];

 % Ranks to test
test_ranks = 1:8;
tol = 1e-16; % Tolerance is small because we are specifying the maximum rank

%% Computing the rounded tensor for each rank
L = length(test_ranks);

errRounding      = zeros(L,S);
errOrthThenRand  = zeros(L,S);
errRandThenOrth  = zeros(L,S);
errTwoSided      = zeros(L,S);
timeRounding     = zeros(L,S);
timeOrthThenRand = zeros(L,S);
timeRandThenOrth = zeros(L,S);
timeTwoSided     = zeros(L, S);
OutRanks         = zeros(L,N+1);

power = 1;
for s = 1:S
    Y = TThilbert(I, true_ranks, power);
    normY = TTnorm(Y);

    for r = 1:L
        ranks = [1; test_ranks(r) * ones(N-1,1); 1];

        tic;
        X = TTrounding(Y, tol, test_ranks(r));
        timeRounding(r,s) = toc;
        errRounding(r,s) = TTnorm(TTaxby(1, Y, -1, X), "OLR") / normY;

        OutRanks(r,:) = TTranks(X);

        tic;
        X = TTrounding_Orthogonalize_then_Randomize(Y, ranks);
        timeOrthThenRand(r,s) = toc;
        errOrthThenRand(r,s) = TTnorm(TTaxby(1, Y, -1, X), "OLR") / normY;

        tic;
        X = TTrounding_Randomize_then_Orthogonalize(Y, ranks);
        timeRandThenOrth(r,s) = toc;
        errRandThenOrth(r,s) = TTnorm(TTaxby(1,Y,-1,X), "OLR") / normY;
        
        tic;
        X = TTrounding_Two_Sided_Randomization(Y, ranks);
        timeTwoSided(r,s) = toc;
        errTwoSided(r,s) = TTnorm(TTaxby(1,Y,-1,X), "OLR") / normY;
    end
end

%% Plot Results
%% Plot for errors
f = figure(1);
f.Position(1:2) = [0,525];
f.Position(3:4) = [525, 350];
clf()


[err, neg, pos] = processErrors(errRounding(:,:));
errorbar(test_ranks, err, neg, pos, 'bo','markersize',10,'linewidth',2)
set(gca, 'YScale', 'log')
hold on
[err, neg, pos] = processErrors(errOrthThenRand(:,:));
errorbar(test_ranks, err, neg, pos, 'ks','markersize',10,'linewidth',2)
[err, neg, pos] = processErrors(errRandThenOrth(:,:));
errorbar(test_ranks, err, neg, pos, 'rx','markersize',10,'linewidth',2)
[err, neg, pos] = processErrors(errTwoSided(:,:));
errorbar(test_ranks, err, neg, pos, '+-','markersize',10,'linewidth',2)
hold off

ylim([10^-15, 1])
xlabel('Maximum Target Rank')
ylabel('Error')
set(gca,'FontSize',16)
legend('TT-Rounding','Orth-then-Rand', 'Rand-then-Orth','location','best')

exportgraphics(f,'figures/Hilbert.eps')

%% Plot for time
f = figure(2);
f.Position(1:2) = [0,525];
f.Position(3:4) = [525, 350];
clf()

[err, neg, pos] = processTime(timeRounding(:,:));
errorbar(test_ranks, err, neg, pos, 'bo','markersize',10,'linewidth',2)
set(gca, 'YScale', 'log')
hold on
[err, neg, pos] = processTime(timeOrthThenRand(:,:));
errorbar(test_ranks, err, neg, pos, 'ks','markersize',10,'linewidth',2)
[err, neg, pos] = processTime(timeRandThenOrth(:,:));
errorbar(test_ranks, err, neg, pos, 'rx','markersize',10,'linewidth',2)
[err, neg, pos] = processTime(timeTwoSided(:,:));
errorbar(test_ranks, err, neg, pos, '+-','markersize',10,'linewidth',2)
hold off

xlabel('Maximum Target Rank')
ylabel('Time')
set(gca,'FontSize',16)
legend('TT-Rounding','Orth-then-Rand', 'Rand-then-Orth','Two-Sided-Rand','location','best')

exportgraphics(f,'figures/Hilbert_time.eps')


%%% Time speedup

% Post-process timings

f = figure(3);
f.Position(1:2) = [0,525];
f.Position(3:4) = [525, 350];
clf()

[speedupLR,         negLR,          posLR]          = computeSpeedup(timeRounding,       timeRounding);
[speedupRandLR,     negRandLR,      posRandLR]      = computeSpeedup(timeOrthThenRand,   timeRounding);
[speedupRandOrth,   negRandOrth,    posRandOrth]    = computeSpeedup(timeRandThenOrth, timeRounding);
[speedupTwoSide,    negTwoSide,     posTwoSide]     = computeSpeedup(timeTwoSided,  timeRounding);


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
legend('TT-Rounding','Orth-then-Rand', 'Rand-then-Orth','Two-Sided-Rand','location','best')

exportgraphics(f,'figures/Hilbert_speedup.eps')




%% Hilbert TT-tensor construction

function X = TThilbert(I,R,p)
% TTHILBERT  create a TT-tensor whose cores are Hilbert tensors.
%
% X = TTHILBERT(I, R, p)  returns a TT-tensor X,
%                   with mode sizes I(1), ..., I(N) and
%                   TT-ranks R(1), ..., R(N + 1).
%                   p is the power of the denominator in the Hilbert tensor.


N = length(I); % number of modes

assert(length(R) == N+1, "Please provide conforming sizes of the modes sizes and the TT ranks\n");
assert(R(1)==1 && R(N+1)==1, "Boundary TT ranks must be 1");

X = cell(N, 1);
for n = 1:N
    X{n} = zeros(R(n), I(n), R(n+1));
    for i=1:R(n)
        for j=1:I(n)
            for k=1:R(n+1)
                X{n}(i,j,k) = 1/(i+j+k-2)^p;
            end
        end
    end
    X{n} = reshape(X{n}, R(n)*I(n), R(n+1));
end

end

%% Custom function to post-process errors
function [error, neg, pos] = processErrors(errors)
    error = median(errors, 2);
    neg = min(errors, [], 2);
    pos = max(errors, [], 2);
end

%% Custom function to post-process time
function [error, neg, pos] = processTime(time)
    error = median(time, 2);
    neg = min(time, [], 2);
    pos = max(time, [], 2);
end

%% Custom function to post-process runtimes and speedup

function [speedup, neg, pos] = computeSpeedup(time, timeref)
    time = time(:,2:end,:);
    timeref = timeref(:,2:end,:);
    speedup = median(timeref, [2,3]) ./ median(time, [2,3]);
    neg = speedup - min(timeref,[],[2,3]) ./ max(time,[],[2,3]);
    pos = max(timeref,[],[2,3]) ./ min(time,[],[2,3]) - speedup;
end