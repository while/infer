% ==============================================================================
%  Beta distribution
% ==============================================================================
-module(beta).

-export([betapdf/3, betacdf/3, betainv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(SQR(X), ((X)*(X))).

% ------------------------------------------------------------------------------
%  betapdf - Probability density function
% ------------------------------------------------------------------------------
betapdf(_,Alpha,_) when Alpha =< 0 -> 
        {error,"Alpha parameter =< 0 in betapdf."};
betapdf(_,_,Beta) when Beta =< 0 -> 
        {error,"Beta parameter =< 0 in betapdf."};
betapdf(X,_,_) when X =< 0.0 orelse X >= 1.0 -> 
        {error,"Bad X parameter in betapdf."};

betapdf(X,Alpha,Beta) ->
        Fac = incgammabeta:gammaln(Alpha+Beta) - incgammabeta:gammaln(Alpha)
              - incgammabeta:gammaln(Beta),
        math:exp((Alpha-1.0)*math:log(X)+(Beta-1.0)*math:log(1.0-X)+Fac).

% ------------------------------------------------------------------------------
%  betacdf - Cumulative distribution function
% ------------------------------------------------------------------------------
betacdf(_,Alpha,_) when Alpha =< 0 -> 
        {error,"Alpha parameter =< 0 in betacdf."};
betacdf(_,_,Beta) when Beta =< 0 -> 
        {error,"Beta parameter =< 0 in betacdf."};
betacdf(X,_,_) when X =< 0.0 orelse X >= 1.0 -> 
        {error,"Bad X parameter in betacdf."};

betacdf(X,Alpha,Beta) ->
        incgammabeta:betai(Alpha,Beta,X).

        
% ------------------------------------------------------------------------------
%  betainv - Inverse cumulative distribution function
% ------------------------------------------------------------------------------
betainv(_,Alpha,_) when Alpha =< 0 -> 
        {error,"Alpha parameter =< 0 in betainv."};
betainv(_,_,Beta) when Beta =< 0 -> 
        {error,"Beta parameter =< 0 in betainv."};

betainv(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};

betainv(P,Alpha,Beta) -> 
        incgammabeta:invbetai(P,Alpha,Beta).




% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

betapdf_test() ->
        ?assertEqual(3.2804999999999933, betapdf(0.1,1,5)),
        ?assertEqual(3.8742048900000183, betapdf(0.1,1,10)),
        ?assertEqual(1.8749999999999938, betapdf(0.5,3,3)),
        ?assertEqual(1.0610329539459684, betapdf(0.9,0.5,0.5)).

betapdf_error_test() ->
        ?assertEqual({error,"Alpha parameter =< 0 in betapdf."},  betapdf(0.0,-1,1)),
        ?assertEqual({error,"Beta parameter =< 0 in betapdf."},  betapdf(0.0,1,-1)),
        ?assertEqual({error,"Bad X parameter in betapdf."},  betapdf(-1.0,1,1)),
        ?assertEqual({error,"Bad X parameter in betapdf."},  betapdf(2.0,1,1)).

betacdf_test() ->
        ?assertEqual(0.4095099999999993, betacdf(0.1,1,5)),
        ?assertEqual(0.6513215599000027, betacdf(0.1,1,10)),
        ?assertEqual(0.5000000000000016, betacdf(0.5,3,3)),
        ?assertEqual(0.7951672353008667, betacdf(0.9,0.5,0.5)).

betacdf_param_test() ->
        ?assertEqual({error,"Alpha parameter =< 0 in betacdf."},  betacdf(0.0,-1,1)),
        ?assertEqual({error,"Beta parameter =< 0 in betacdf."},  betacdf(0.0,1,-1)),
        ?assertEqual({error,"Bad X parameter in betacdf."},  betacdf(-1.0,1,1)),
        ?assertEqual({error,"Bad X parameter in betacdf."},  betacdf(2.0,1,1)).

betainv_test() ->
        ?assertEqual(0.020851637639023268, betainv(0.1,1,5)),
        ?assertEqual(0.010480741793785553, betainv(0.1,1,10)),
        ?assertEqual(0.499999999999999, betainv(0.5,3,3)),
        ?assertEqual(0.610181649469041, betainv(0.7,3,3)),
        ?assertEqual(0.9755282581475768, betainv(0.9,0.5,0.5)).

betainv_param_test() ->
        ?assertEqual({error,"Alpha parameter =< 0 in betainv."},  betainv(0.0,-1,1)),
        ?assertEqual({error,"Beta parameter =< 0 in betainv."},  betainv(0.0,1,-1)),
        ?assertEqual({error,"Invalid probability"},  betainv(-0.1,1,1)).

-endif.




