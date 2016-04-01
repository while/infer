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
%  Student-t test
% ------------------------------------------------------------------------------
ttest(Xs, Ys) ->
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
        {{t, T}, {pval, Prob}}.


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

ttest_test() ->
        error.

-endif.
