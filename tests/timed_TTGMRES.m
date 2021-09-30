function [x, ranks, time_RoundSum, time_Operator, time_Prec, time_Other] = timed_TTGMRES(A, b, tol, maxit, prec, tols, rounding_sum)

time_RoundSum = 0;
time_Operator = 0;
time_Prec = 0;
time_Other = 0;

k = 0;
normb = TTnorm(b);
V = cell(maxit, 1);

tStart = tic;
    r = b;
    normres = TTnorm(r);
    V{1} = TTscale(r, 1/normres);
time_Other =  time_Other + toc(tStart);


tStart = tic;
    H = zeros(maxit + 1, maxit);
    res = zeros(maxit + 1, 1);
    res(1) = normres;

    rK = zeros(maxit + 1, 1);
    rK(1) = normres;
    
    ranks = cell(maxit + 2, 1);
    ranks{1} = TTranks(V{1});
time_Other =  time_Other + toc(tStart);


while(k < maxit && normres > tol * normb)
tStart = tic;
    k = k + 1;
    w = prec(V{k});
time_Prec = time_Prec + toc(tStart);

tStart = tic;
    w = A(w);
time_Operator = time_Operator + toc(tStart);

tStart = tic;
    delta = 0.01*tol/(normres * normb);
    w = rounding_sum(w, ones(length(w)), delta);
time_RoundSum = time_RoundSum + toc(tStart);
    
tStart = tic;
    h = zeros(k + 1, 1);
    for j = 1 : k
        h(j) = TTdot(V{j}, w);
    end
time_Other = time_Other + toc(tStart);
tStart = tic;
    V{k+1} = w;
    w = rounding_sum(V(1:k+1), [-h(1:k);1], delta);
time_RoundSum = time_RoundSum + toc(tStart);
tStart = tic;
    h(k + 1) = TTnorm(w);
    if(h(k + 1) < 2e-16)
        keyboard
    end
    V{k + 1} = TTscale(w, 1/h(k + 1));
    ranks{k + 1} = TTranks(V{k+1});
    
    H(1 : k + 1, k) = h;
    [Q, R] = qr(H(1 : k + 1, 1 : k));
    t = Q' * rK(1 : k + 1);
    res(k + 1) = abs(t(k + 1));
    normres = res(k + 1);
    
    t = Q' * rK(1 : k + 1);
time_Other = time_Other + toc(tStart);
end

tStart = tic;
    y = R(1 : k, 1 : k)\t(1 : k);
time_Other = time_Other + toc(tStart);

tStart = tic;
    x = rounding_sum(V(1:k), y(1:k), tols);
time_RoundSum = time_RoundSum + toc(tStart);

tStart = tic;
    x = prec(x);
time_Prec = time_Prec + toc(tStart);

tStart = tic;
    x = TTrounding(x, tols);
time_RoundSum = time_RoundSum + toc(tStart);

tStart = tic;
    ranks{k + 2} = TTranks(x);
    ranks = ranks(1 : k + 2);
time_Other = time_Other + toc(tStart);
