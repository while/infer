% ==============================================================================
%  F distribution
% ==============================================================================
-module(fdist).

-export([fdistpdf/3, fdistcdf/3, fdistinv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

% ------------------------------------------------------------------------------
%  fdistpdf - Probability density function
% ------------------------------------------------------------------------------
fdistpdf(_,Nu1,_) when Nu1 =< 0 -> {error,"Nu1 parameter =< 0 in F-dist."};
fdistpdf(_,_,Nu2) when Nu2 =< 0 -> {error,"Nu2 parameter =< 0 in F-dist."};

fdistpdf(X,_,_) when X < 0 -> 0.0;
fdistpdf(X,Nu1,Nu2) ->
        Fac = 0.5*(Nu1*math:log(Nu1) + Nu2*math:log(Nu2))
              + incgammabeta:gammaln(0.5*(Nu1+Nu2))
              - incgammabeta:gammaln(0.5*Nu1)
              - incgammabeta:gammaln(0.5*Nu2),
        math:exp((0.5*Nu1-1.0)*math:log(X)-0.5*(Nu1+Nu2)*math:log(Nu2+Nu1*X)+Fac).

% ------------------------------------------------------------------------------
%  fdistcdf - Cumulative distribution function
% ------------------------------------------------------------------------------
fdistcdf(_,Nu1,_) when Nu1 =< 0 -> {error,"Nu1 parameter =< 0 in F-dist."};
fdistcdf(_,_,Nu2) when Nu2 =< 0 -> {error,"Nu2 parameter =< 0 in F-dist."};

fdistcdf(X,_,_) when X < 0 -> 0.0;
fdistcdf(X,Nu1,Nu2) ->
        incgammabeta:betai(0.5*Nu1,0.5*Nu2,Nu1*X/(Nu2+Nu1*X)).

        
% ------------------------------------------------------------------------------
%  fdistinv - Inverse cumulative distribution function
% ------------------------------------------------------------------------------
fdistinv(_,Nu1,_) when Nu1 =< 0 -> {error,"Nu1 parameter =< 0 in F-dist."};
fdistinv(_,_,Nu2) when Nu2 =< 0 -> {error,"Nu2 parameter =< 0 in F-dist."};

fdistinv(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};

fdistinv(X,_,_) when X =< 0 -> 0.0;
fdistinv(P,Nu1,Nu2) ->
        InvB = incgammabeta:invbetai(P,0.5*Nu1,0.5*Nu2),
        Nu2*InvB/(Nu1*(1.0-InvB)).



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

fdistpdf_test() ->
        ?assertEqual(0.0, fdistpdf(-1.0,1,1)),
        ?assertEqual(0.1591549430918952, fdistpdf(1,1,1)),
        ?assertEqual(0.07502635967975876, fdistpdf(2,1,1)),
        ?assertEqual(0.19245008972987535, fdistpdf(1,2,1)),
        ?assertEqual(0.19245008972987535, fdistpdf(1,1,2)),
        ?assertEqual(0.010391328106475825, fdistpdf(10,2,1)),
        ?assertEqual(0.00614445752989748, fdistpdf(10,2,3)),
        ?assertEqual(0.10519242090293805, fdistpdf(2,10,1)),
        ?assertEqual(0.0027551055016847756, fdistpdf(2,100,100)).

fdistpdf_error_test() ->
        ?assertEqual({error,"Nu1 parameter =< 0 in F-dist."},  fdistpdf(0.0,-1,1)),
        ?assertEqual({error,"Nu2 parameter =< 0 in F-dist."},  fdistpdf(0.0,1,-1)).

fdistcdf_test() ->
        ?assertEqual(0.0, fdistcdf(-1.0,1,1)),
        ?assertEqual(0.5000000000000007, fdistcdf(1,1,1)),
        ?assertEqual(0.6081734479693931, fdistcdf(2,1,1)),
        ?assertEqual(0.42264973081037405, fdistcdf(1,2,1)),
        ?assertEqual(0.7817821097640074, fdistcdf(10,2,1)),
        ?assertEqual(0.4956475043831201, fdistcdf(2,10,1)).

fdistcdf_error_test() ->
        ?assertEqual({error,"Nu1 parameter =< 0 in F-dist."},  fdistcdf(0.0,-1,1)),
        ?assertEqual({error,"Nu2 parameter =< 0 in F-dist."},  fdistcdf(0.0,1,-1)).

fdistinv_test() ->
        ?assertEqual(0.0, fdistinv(0.0,1,1)),
        ?assertEqual(0.025085630936916625, fdistinv(0.1,1,1)),
        ?assertEqual(0.1172839506172838, fdistinv(0.1,2,1)),
        ?assertEqual(0.9999999999999966, fdistinv(0.5,1,1)),
        ?assertEqual(1.5000000000000018, fdistinv(0.5,2,1)),
        ?assertEqual(60.19498034404559, fdistinv(0.9,10,1)).

fdistinv_error_test() ->
        ?assertEqual({error,"Nu1 parameter =< 0 in F-dist."},  fdistinv(0.0,-1,1)),
        ?assertEqual({error,"Nu2 parameter =< 0 in F-dist."},  fdistinv(0.0,1,-1)),
        ?assertEqual({error,"Invalid probability"},  fdistinv(-0.1,1,1)).

-endif.




