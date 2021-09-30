%% TTGMRES test
clear;
clc;
maxNumCompThreads(4);
close("all")


%% Global parameters

% Number of independent runs
S = 6;

% Number of parameters ranges from 2^2 to 2^l
l = 8; 
N = 2.^(3:l);


%% Load FEM matrices

A0 = readcoomat('data/A0.txt');
A1 = readcoomat('data/A1.txt');
A2 = readcoomat('data/A2.txt');
A3 = readcoomat('data/A3.txt');
A4 = readcoomat('data/A4.txt');

n1 = length(A0);
A = speye(n1);

a0 = readb('data/b.txt');

%% Run and time experiments

sumTimeNormalGMRES = zeros(length(N),S);
sumTimeRandGMRES = zeros(length(N),S);
opTimeNormalGMRES = zeros(length(N),S);
opTimeRandGMRES = zeros(length(N),S);
precTimeNormalGMRES = zeros(length(N),S);
precTimeRandGMRES = zeros(length(N),S);
remTimeNormalGMRES = zeros(length(N),S);
remTimeRandGMRES = zeros(length(N),S);
runtimeNormalGMRES = zeros(length(N),S);
runtimeRandGMRES = zeros(length(N),S);
ranksrand = cell(length(N));
ranksnormal = cell(length(N));

for i = 1:length(N)
    n = N(i);

    % Distribution of parameters
    MinD = 1e0;
    MaxD = 1e1;

    d = linspace(MinD, MaxD, n);

    % Matrices and operator in Kronecker form
    D = sparse(1:n, 1:n, d);

    c = ones(n,1);
    B = speye(n);

    b = {a0;c;c;c;c};

    A = {A0, A1, A2, A3, A4;
         B,  D,  B,  B,  B;
         B,  B,  D,  B,  B;
         B,  B,  B,  D,  B;
         B,  B,  B,  B,  D};

    % TT-GMRES parameter setup
    tol      = 1e-8;
    maxit    = 50;
    [L,U,p]  = lu(A0 + 1*(A1 + A2 + A3 + A4), 'vector');
    prec     = @(x) Preconditioner(L,U,p,x);
    tols     = tol*1e-2;
    normb    = TTnorm(b);
    b        = TTscale(b, 1/normb);
    normb    = TTnorm(b);

    Op            = @(x)       TTsummandsKronOp(A, x);
    deterministic = @(W,a,tol) TTsum(W, a, tol);
    randomized    = @(W,a,tol) TTsum_Randomize_then_Orthogonalize(W, a, tol);

    for s=1:S
        fprintf("*********\nRandomized TT-GMRES run - n2 = %i\t", n);
        tStart = tic;
        [x, ranksR, sumtime, optime, prectime, remtime] = timed_TTGMRES(Op, b, tol, maxit, prec, tols, randomized);
        tEnd = toc(tStart);
        fprintf("Elapsed time is %.2f seconds.\n", tEnd);
        t = TTsummandsKronOp(A,x);
        r = t{1};
        for j=2:length(t)
            r = TTaxby(1., r, 1., t{j});
        end
        r = TTaxby(1., r, -1., b);
        fprintf("*********\nnorm real res = %e\t", TTnorm(r, "OLR")/normb);
        fprintf("rounding and sum time = %.2f\n********\n", sumtime);

        for j = 1 : size(ranksR,1)
            ranksrand{i}(j,s) = max(ranksR{j});
        end
        sumTimeRandGMRES(i,s)   = sumtime;
        opTimeRandGMRES(i,s)    = optime;
        precTimeRandGMRES(i,s)  = prectime;
        remTimeRandGMRES(i,s)   = remtime;
        runtimeRandGMRES(i,s)   = tEnd;

        fprintf("*********\nDeterministic TT-GMRES run - n2 = %i\t", n);
        tStart = tic;
        [x, ranksN, sumtime, optime, prectime, remtime] = timed_TTGMRES(Op, b, tol, maxit, prec, tols, deterministic);
        tEnd = toc(tStart);
        fprintf("Elapsed time is %.2f seconds.\n", tEnd);
        t = TTsummandsKronOp(A,x);
        r = t{1};
        for j=2:length(t)
            r = TTaxby(1., r, 1., t{j});
        end
        r = TTaxby(1., r, -1., b);
        fprintf("*********\nnorm real res = %e\t", TTnorm(r, "OLR")/normb);
        fprintf("rounding and sum time = %.2f\n*********\n", sumtime);
        
        for j = 1 : size(ranksN,1)
            ranksnormal{i}(j,s) = max(ranksN{j});
        end
        sumTimeNormalGMRES(i,s)  = sumtime;
        opTimeNormalGMRES(i,s)   = optime;
        precTimeNormalGMRES(i,s) = prectime;
        remTimeNormalGMRES(i,s)  = remtime;
        runtimeNormalGMRES(i,s)  = tEnd;
    end
end

% Discard results of first run
sumTimeNormalGMRES = sumTimeNormalGMRES(:,2:end);
otherTimeNormalGMRES = opTimeNormalGMRES+precTimeNormalGMRES+remTimeNormalGMRES;
otherTimeNormalGMRES = otherTimeNormalGMRES(:,2:end);
runtimeNormalGMRES = runtimeNormalGMRES(:,2:end);

sumTimeRandGMRES = sumTimeRandGMRES(:,2:end);
otherTimeRandGMRES = opTimeRandGMRES+precTimeRandGMRES+remTimeRandGMRES;
otherTimeRandGMRES = otherTimeRandGMRES(:,2:end);
runtimeRandGMRES = runtimeRandGMRES(:,2:end);


%% Plot Results

close("all")

f = figure(1);
f.Position(1:2) = [0, 0];
f.Position(3:4) = [1050,350];
X = categorical(N);

tiledlayout(1,2);
nexttile;
bar(X, [mean(sumTimeRandGMRES,2), mean(otherTimeRandGMRES,2)], 'stacked')
title('Randomized TT-GMRES')
xlabel('Number of parameter samples', 'FontSize', 18)
ylabel('Average Runtime (s)', 'FontSize', 18)
legend('RandRoundSum', 'Other', 'location', 'northwest')
set(gca,'FontSize',16)

nexttile;
bar(X, [mean(sumTimeNormalGMRES,2), mean(otherTimeNormalGMRES,2)], 'stacked')
title('Naive TT-GMRES')
xlabel('Number of parameter samples', 'FontSize', 18)
ylabel('Average Runtime (s)', 'FontSize', 18)
legend('TT-Sum + Round', 'Other', 'location', 'northwest')
set(gca,'FontSize',16)
exportgraphics(f, 'figures/TTGMRES_Runtime.eps')


f = figure(2);
f.Position(1:2) = [0,525];
f.Position(3:4) = [525, 350];
speedup = median(runtimeNormalGMRES, 2) ./ median(runtimeRandGMRES, 2);
neg = speedup - min(runtimeNormalGMRES,[],2)./ max(runtimeRandGMRES,[],2);
pos = max(runtimeNormalGMRES,[],2)./ min(runtimeRandGMRES,[],2) - speedup;
errorbar(N, speedup, neg, pos,'r-s','markersize',10,'linewidth',2);
hold on;
speedup = median(sumTimeNormalGMRES, 2) ./ median(sumTimeRandGMRES, 2);
neg = speedup - min(sumTimeNormalGMRES,[],2)./ max(sumTimeRandGMRES,[],2);
pos = max(sumTimeNormalGMRES,[],2)./ min(sumTimeRandGMRES,[],2) - speedup;
errorbar(N, speedup, neg, pos,'b-o','markersize',10,'linewidth',2);
hold off;
legend('Total speedup', 'TT-Sum + Round', 'location', 'northwest')

set(gca, 'XScale', 'log')
xlabel('Number of parameter samples (log scale)', 'FontSize', 18)
ylabel('Speedup', 'FontSize', 18)
xticks(N);
set(gca,'FontSize',16)
exportgraphics(f, 'figures/TTGMRES_Speedup.eps')

f = figure(3);
clf();
f.Position(1:2) = [525,525];
f.Position(3:4) = [525, 350];
ranks = max(ranksrand{1}(1:end-1,:),[],2);
p1 = plot(ranks, 'bo','markersize',10,'linewidth',2);
hold on
ranks = max(ranksnormal{1}(1:end-1,:),[],2);
p2 = plot(ranks, 'r+','markersize',10,'linewidth',2);

for j=2:length(N)
    ranks = max(ranksrand{j}(1:end-1,:),[],2);
    plot(ranks, 'bo', 'markersize',10,'linewidth',2)
end
for j=2:length(N)
    ranks = max(ranksnormal{j}(1:end-1,:),[],2);
    plot(ranks, 'r+', 'markersize',10,'linewidth',2)
end
hold off

legend([p1, p2], 'Randomized TT-GMRES', 'Naive TT-GMRES', 'location', 'southeast')
xlabel('TT-GMRES iteration number', 'FontSize', 18)
ylabel('TT-rank', 'FontSize', 18)
set(gca,'FontSize',16)
exportgraphics(f, 'figures/TTGMRES_Ranks.eps')

beep;
pause(1)
beep;


%% Preconditioner
function X = Preconditioner(L,U,p, Y)
    X = Y;
    X{1} = U\(L\Y{1}(p,:));
end
