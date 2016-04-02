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
-define(assertEqualDigits(X,Y,D), ?assert(abs((X)-(Y)) < math:pow(10,-(D)))).
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

meanvar(Xs) ->
        N = length(Xs),
        {SumX,SumX2} = lists:foldl(fun(X,{AccX,AccX2}) -> {AccX+X, AccX2+?SQR(X)} end, {0,0}, Xs),
        {SumX/N, (SumX2 - ?SQR(SumX)/N)/(N-1)}.
        

% ------------------------------------------------------------------------------
%  Student's t test (Equal or unequal sample sizes, equal variance)
% ------------------------------------------------------------------------------
sttest(Xs, Ys) ->
        Df = length(Xs) + length(Ys) - 2,
        Nx = length(Xs),
        Ny = length(Ys),
        {MuX,VarX} = meanvar(Xs),
        {MuY,VarY} = meanvar(Ys),
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
        {MuX,VarX} = meanvar(Xs),
        {MuY,VarY} = meanvar(Ys),
        T = (MuX - MuY)/math:sqrt(VarX/Nx + VarY/Ny),
        Df=?SQR(VarX/Nx + VarY/Ny)/(?SQR(VarX/Nx)/(Nx-1) + ?SQR(VarY/Ny)/(Ny-1)),
        Prob = 2*tdist:cdf(T,0.0,1.0,Df),
        {{statistic, T}, {p_value, Prob}}.


% ------------------------------------------------------------------------------
%  F test
% ------------------------------------------------------------------------------
ftest(Xs,Ys) -> 
        Nx = length(Xs),
        Ny = length(Ys),
        VarX = var(Xs),
        VarY = var(Ys),
        {F,Df1,Df2} = if VarX >  VarY -> {VarX/VarY,Nx-1,Ny-1}
                       ; VarX =< VarY -> {VarY/VarX,Ny-1,Nx-1}
                      end,
        Prob = 2.0*fdist:cdf(F,Df1,Df2),
        {{statistic, F},
         {p_value, if Prob >  1.0 -> 2.0-Prob;
                      Prob =< 1.0 -> Prob end}}.


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

meanvar_test() ->
        ?assertEqual({1.0,0.0}, meanvar([1,1,1])),
        ?assertEqual({3.5,3.5}, meanvar([1,2,3,4,5,6])),
        ?assertEqual({2.0,1.0}, meanvar([1,2,3])).

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
        ?assertEqual({{statistic, -0.20174824628379245},
                      {p_value, 0.842690649556616}},
                     wttest(Xs,Ys)).

ftest_test() ->
        % Normal(0,1) variates
        Xs = [-0.205851093997324, -1.82828630625735, -0.180571751075289,
              -1.75489242033749, 0.344844472657044, 0.451987442089703,
              0.471725272234911, 0.7032831790146, -1.56967983963718,
              0.245924162151692, 0.153554768985877, -0.106792500229723,
              0.285680472688479, -0.0362077886137734, -0.912868712453373,
              -1.50022366711904, 0.742218361063315, -1.86029791668071,
              0.722963762247153, 0.396927361496996],
        % Normal(0,2^2) variates
        Ys = [-1.2740990151509, 0.723313469353654, -2.80484043451352,
              -1.90781044419281, 1.20992674834345, 2.61046965670307,
              1.07801428940928, -1.2794530489474, -0.457263911622332,
              0.214371515389186, 0.997417828931748, 2.68453565514012,
              3.07569976474604, -1.45615376944145, 1.3303228831093,
              -2.32969983774504, -0.795415018628534, 1.21099122567277,
              -2.48388806607184, -2.14442679550042],
        {{statistic,F},{p_value,Pval}} = ftest(Xs,Ys),
        ?assertEqual(3.924710587503716,F),
        ?assertEqual(0.004532196830164814,Pval).

ftest_equal_test() ->
        % Normal(0,1) variates
        Xs = [-0.205851093997324, -1.82828630625735, -0.180571751075289,
              -1.75489242033749, 0.344844472657044, 0.451987442089703,
              0.471725272234911, 0.7032831790146, -1.56967983963718,
              0.245924162151692, 0.153554768985877, -0.106792500229723,
              0.285680472688479, -0.0362077886137734, -0.912868712453373,
              -1.50022366711904, 0.742218361063315, -1.86029791668071,
              0.722963762247153, 0.396927361496996],
        {{statistic,F},{p_value,Pval}} = ftest(Xs,Xs),
        ?assertEqualDigits(1.0,F,14),
        ?assertEqualDigits(1.0,Pval,14).

-endif.
