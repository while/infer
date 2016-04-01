% ==============================================================================
%  Statistical hypothesis tests
% ==============================================================================
-module(stattest).

-export([]).

% TODO: Remove export_all!
-compile([export_all]).

-define(SQR(X), ((X)*(X))).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

% ------------------------------------------------------------------------------
%  Simple helpers
% ------------------------------------------------------------------------------
mean(Xs) ->
        lists:sum(Xs)/length(Xs).

var(Xs) ->
        Mu = mean(Xs),
        N = length(Xs),
        lists:foldl(fun(X,Acc) -> Acc + ?SQR(X-Mu) end, 0, Xs)/(N-1).

% ------------------------------------------------------------------------------
%  Student's t test (Equal or unequal sample sizes, equal variance)
% ------------------------------------------------------------------------------
sttest(Xs, Ys) ->
        Df = length(Xs) + length(Ys) - 2,
        Nx = length(Xs),
        Ny = length(Ys),
        MuX = mean(Xs),
        MuY = mean(Ys),
        VarX = var(Xs),
        VarY = var(Ys),
        Sxy = ((Nx-1)*VarX+(Ny-1)*VarY)/Df,
        T = (MuX - MuY)/math:sqrt(Sxy*(1.0/Nx + 1.0/Ny)),
        Prob = 2*tdist:cdf(T,0.0,1.0,Df),
        {{statistic, T}, {p_value, Prob}}.

% ------------------------------------------------------------------------------
%  Welch's t test (Equal or unequal sample sizes, unequal variance)
% ------------------------------------------------------------------------------
wttest(Xs, Ys) ->
        Nx = length(Xs),
        Ny = length(Ys),
        MuX = mean(Xs),
        MuY = mean(Ys),
        VarX = var(Xs),
        VarY = var(Ys),
        T = (MuX - MuY)/math:sqrt(VarX/Nx + VarY/Ny),
        Df=?SQR(VarX/Nx + VarY/Ny)/(?SQR(VarX/Nx)/(Nx-1) + ?SQR(VarY/Ny)/(Ny-1)),
        Prob = 2*tdist:cdf(T,0.0,1.0,Df),
        {{statistic, T}, {p_value, Prob}}.


% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

mean_test() ->
        ?assertEqual(1.0, mean([1,1,1])),
        ?assertEqual(3.5, mean([1,2,3,4,5,6])),
        ?assertEqual(2.0, mean([1,2,3])).

var_test() ->
        ?assertEqual(0.0, var([1,1,1])),
        ?assertEqual(3.5, var([1,2,3,4,5,6])),
        ?assertEqual(1.0, var([1,2,3])).

sttest_test() ->
        Xs = [-1.20093942067179, 0.835770403625221, -0.773530756127177,
              -0.276594156215122, -0.271174108906289, 0.0974272885649251,
              0.404071994945201, 0.330078899108719, 0.499200718682945,
              -0.49703613199271],
        Ys = [-0.327520116434318, -0.195152014205557, 1.06504721546445,
              -0.658057764201224, -0.77539891485522, 0.695209191948099,
              1.30198002432061, 0.215325106726302, 0.347777749766435,
              -1.80090200342619],
        ?assertEqual({{statistic, -0.20174824628379243},
                      {p_value, 0.8423760415858366}},
                     sttest(Xs,Ys)).


wttest_test() ->
        Xs = [-1.20093942067179, 0.835770403625221, -0.773530756127177,
              -0.276594156215122, -0.271174108906289, 0.0974272885649251,
              0.404071994945201, 0.330078899108719, 0.499200718682945,
              -0.49703613199271],
        Ys = [-0.327520116434318, -0.195152014205557, 1.06504721546445,
              -0.658057764201224, -0.77539891485522, 0.695209191948099,
              1.30198002432061, 0.215325106726302, 0.347777749766435,
              -1.80090200342619],
        ?assertEqual({{statistic, -0.20174824628379243},
                      {p_value, 0.8426906495566168}},
                     wttest(Xs,Ys)).

-endif.
