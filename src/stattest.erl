% ==============================================================================
%  Statistical hypothesis tests
% ==============================================================================
-module(stattest).

-export([]).

% TODO: Remove export_all!
-compile([export_all]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

% ------------------------------------------------------------------------------
%  Simple helpers
% ------------------------------------------------------------------------------
mean(Xs) ->
        lists:sum(Xs)/length(Xs).

var(Xs) ->
        lists:foldl(fun(X,Acc) -> Acc + X*X end, 0, Xs)/(length(Xs) -1).

% ------------------------------------------------------------------------------
%  t-test
% ------------------------------------------------------------------------------
ttest(Xs, Ys) ->
                Df = length(Xs) + length(Ys) - 2,
                Nx = length(Xs),
                Ny = length(Ys),
                MuX = mean(Xs),
                MuY = mean(Ys),
                VarX = var(Xs),
                VarY = var(Ys),
                Svar=((Nx-1)*VarX+(Ny-1)*VarY)/Df,
                T = (MuX - MuY)/math:sqrt(Svar*(1.0/Nx + 1.0/Ny)),
                Prob = incgammabeta:betai(0.5*Df,0.5,Df/(Df+T*T)),
                {{t, T}, {pval, Prob}}.

