% ==============================================================================
%  Gamma distribution
% ==============================================================================
-module(gamma).

-export([gampdf/3, gamcdf/3, gaminv/3, gammaln/1]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


% ------------------------------------------------------------------------------
%  gampdf - Uniform probability density function
% ------------------------------------------------------------------------------
gampdf(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
gampdf(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

gampdf(X,_,_) when X =< 0 -> 0.0;
gampdf(X,Alpha,Beta) ->
        Fac = Alpha*math:log(Beta)-gammaln(Alpha),
        math:exp(-Beta*X+(Alpha-1)*math:log(X)+Fac).

% ------------------------------------------------------------------------------
%  gamcdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
gamcdf(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
gamcdf(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

gamcdf(X,_,_) when X =< 0 -> 0.0;
gamcdf(X,Alpha,Beta) ->
        gammap(Alpha,Beta*X).

        
% ------------------------------------------------------------------------------
%  gaminv - Inverse gamma distribution function
% ------------------------------------------------------------------------------
gaminv(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
gaminv(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

gaminv(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
gaminv(P,Alpha,Beta) ->
        invgammap(P,Alpha)/Beta.


% ------------------------------------------------------------------------------
%  Gamma functions
% ------------------------------------------------------------------------------
gammaln(X) when X =< 0.0 -> {error, "Bad arg to gammln."};
gammaln(X) -> 
        Coef = [ 57.1562356658629235,    -59.5979603554754912, 
                 14.1360979747417471,     -0.491913816097620199,
                  0.339946499848118887e-4, 0.465236289270485756e-4,
                 -0.983744753048795646e-4, 0.158088703224912494e-3,
                 -0.210264441724104883e-3, 0.217439618115212643e-3,
                 -0.164318106536763890e-3, 0.844182239838527433e-4,
                 -0.261908384015814087e-4, 0.368991826595316234e-5 ],

        Tmp1 = X+5.24218750000000000,
        Tmp2 = (X+0.5)*math:log(Tmp1)-Tmp1,
        Ys = [X+I || I <- lists:seq(1,14)],
        Ser = lists:sum([0.999999999999997092|[C/Y || {C,Y} <- lists:zip(Coef,Ys)]]),
        Tmp2+math:log(2.5066282746310005*Ser/X).

gammap(_,_) -> not_implemented.
invgammap(_,_) -> not_implemented.

% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

gammaln_test() ->
        ?assertEqual(2.036327503417712, gammaln(0.123)),
        ?assertEqual(-0.09447840768115956, gammaln(1.234)),
        ?assertEqual(18.734347511936445, gammaln(12.5)).

gampdf_test() ->
        ?assertEqual(0.0, gampdf(0.0,2,2)),
        ?assertEqual(0.2651846416468157, gampdf(2.3,3,1)),
        ?assertEqual(0.5413411329464507, gampdf(1.0,2,2)).

gamcdf_test() ->
        ?assertEqual(0.0, gamcdf(0.0,1,1)).

gaminv_test() ->
        ?assertEqual(0.0, gaminv(0.0,1,1)).



-endif.



