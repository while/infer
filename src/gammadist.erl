% ==============================================================================
%  Gamma distribution
% ==============================================================================
-module(gammadist).

-export([gampdf/3, gamcdf/3, gaminv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(FPMIN, 1.00208418000449e-292).
-define(EPS, 2.22044604925031e-16).

% ------------------------------------------------------------------------------
%  gampdf - Uniform probability density function
% ------------------------------------------------------------------------------
gampdf(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
gampdf(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

gampdf(X,_,_) when X =< 0 -> 0.0;
gampdf(X,Alpha,Beta) ->
        Fac = Alpha*math:log(Beta)-incgammabeta:gammaln(Alpha),
        math:exp(-Beta*X+(Alpha-1)*math:log(X)+Fac).

% ------------------------------------------------------------------------------
%  gamcdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
gamcdf(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
gamcdf(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

gamcdf(X,_,_) when X =< 0 -> 0.0;
gamcdf(X,Alpha,Beta) ->
        incgammabeta:gammap(Alpha,Beta*X).

        
% ------------------------------------------------------------------------------
%  gaminv - Inverse gamma distribution function
% ------------------------------------------------------------------------------
gaminv(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
gaminv(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

gaminv(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
gaminv(P,Alpha,Beta) ->
        incgammabeta:invgammap(P,Alpha)/Beta.



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

gampdf_test() ->
        ?assertEqual(0.0, gampdf(0.0,2,2)),
        ?assertEqual(0.2651846416468157, gampdf(2.3,3,1)),
        ?assertEqual(1.9199765904692206e-04, gampdf(20,20,2)),
        ?assertEqual(0.5413411329464507, gampdf(1.0,2,2)).

gampdf_error_test() ->
        ?assertEqual({error,"Alpha parameter =< 0."}, gampdf(0.0,-1, 1)),
        ?assertEqual({error,"Beta parameter =< 0."},  gampdf(0.0, 1,-1)).

gamcdf_test() ->
        ?assertEqual(0.0, gamcdf(0.0,1,1)),
        ?assertEqual(0.5132987982791589, gamcdf(100,100,1)),
        ?assertEqual(0.5297427331607533, gamcdf(20,20,1)),
        ?assertEqual(0.9998236971022615, gamcdf(20,20,2)),
        ?assertEqual(0.12478121503252411, gamcdf(15,20,1)),
        ?assertEqual(0.9595723180054873, gamcdf(5,2,1)).

gamcdf_error_test() ->
        ?assertEqual({error,"Alpha parameter =< 0."}, gamcdf(0.0,-1, 1)),
        ?assertEqual({error,"Beta parameter =< 0."},  gamcdf(0.0, 1,-1)).

gaminv_test() ->
        ?assertEqual(0.0, gaminv(0.0,1,1)),
        ?assertEqual(0.6931471805599453, gaminv(0.5,1,1)),
        ?assertEqual(0.22746821155978653, gaminv(0.5,0.5,1)),
        ?assertEqual(4.670908882795983, gaminv(0.5,5,1)),
        ?assertEqual(2.432591025962664, gaminv(0.1,5,1)),
        ?assertEqual(4.865182051925328, gaminv(0.1,5,1/2)),
        ?assertEqual(1.216295512981332, gaminv(0.1,5,2)),
        ?assertEqual(3.9967947930263152, gaminv(0.9,5,2)).

-endif.



