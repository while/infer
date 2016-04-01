% ==============================================================================
%  F distribution
% ==============================================================================
-module(fdist).

-export([pdf/3, cdf/3, invcdf/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

% ------------------------------------------------------------------------------
%  pdf - Probability density function
% ------------------------------------------------------------------------------
pdf(_,Nu1,_) when Nu1 =< 0 -> {error,"Nu1 parameter =< 0 in F-dist."};
pdf(_,_,Nu2) when Nu2 =< 0 -> {error,"Nu2 parameter =< 0 in F-dist."};

pdf(X,_,_) when X < 0 -> 0.0;
pdf(X,Nu1,Nu2) ->
        Fac = 0.5*(Nu1*math:log(Nu1) + Nu2*math:log(Nu2))
              + incgammabeta:gammaln(0.5*(Nu1+Nu2))
              - incgammabeta:gammaln(0.5*Nu1)
              - incgammabeta:gammaln(0.5*Nu2),
        math:exp((0.5*Nu1-1.0)*math:log(X)-0.5*(Nu1+Nu2)*math:log(Nu2+Nu1*X)+Fac).

% ------------------------------------------------------------------------------
%  cdf - Cumulative distribution function
% ------------------------------------------------------------------------------
cdf(_,Nu1,_) when Nu1 =< 0 -> {error,"Nu1 parameter =< 0 in F-dist."};
cdf(_,_,Nu2) when Nu2 =< 0 -> {error,"Nu2 parameter =< 0 in F-dist."};

cdf(X,_,_) when X < 0 -> 0.0;
cdf(X,Nu1,Nu2) ->
        incgammabeta:betai(0.5*Nu1,0.5*Nu2,Nu1*X/(Nu2+Nu1*X)).

        
% ------------------------------------------------------------------------------
%  invcdf - Inverse cumulative distribution function
% ------------------------------------------------------------------------------
invcdf(_,Nu1,_) when Nu1 =< 0 -> {error,"Nu1 parameter =< 0 in F-dist."};
invcdf(_,_,Nu2) when Nu2 =< 0 -> {error,"Nu2 parameter =< 0 in F-dist."};

invcdf(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};

invcdf(X,_,_) when X =< 0 -> 0.0;
invcdf(P,Nu1,Nu2) ->
        InvB = incgammabeta:invbetai(P,0.5*Nu1,0.5*Nu2),
        Nu2*InvB/(Nu1*(1.0-InvB)).



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

pdf_test() ->
        ?assertEqual(0.0, pdf(-1.0,1,1)),
        ?assertEqual(0.1591549430918952, pdf(1,1,1)),
        ?assertEqual(0.07502635967975876, pdf(2,1,1)),
        ?assertEqual(0.19245008972987535, pdf(1,2,1)),
        ?assertEqual(0.19245008972987535, pdf(1,1,2)),
        ?assertEqual(0.010391328106475825, pdf(10,2,1)),
        ?assertEqual(0.00614445752989748, pdf(10,2,3)),
        ?assertEqual(0.10519242090293805, pdf(2,10,1)),
        ?assertEqual(0.0027551055016847756, pdf(2,100,100)).

pdf_error_test() ->
        ?assertEqual({error,"Nu1 parameter =< 0 in F-dist."},  pdf(0.0,-1,1)),
        ?assertEqual({error,"Nu2 parameter =< 0 in F-dist."},  pdf(0.0,1,-1)).

cdf_test() ->
        ?assertEqual(0.0, cdf(-1.0,1,1)),
        ?assertEqual(0.5000000000000007, cdf(1,1,1)),
        ?assertEqual(0.6081734479693931, cdf(2,1,1)),
        ?assertEqual(0.42264973081037405, cdf(1,2,1)),
        ?assertEqual(0.7817821097640074, cdf(10,2,1)),
        ?assertEqual(0.4956475043831201, cdf(2,10,1)).

cdf_error_test() ->
        ?assertEqual({error,"Nu1 parameter =< 0 in F-dist."},  cdf(0.0,-1,1)),
        ?assertEqual({error,"Nu2 parameter =< 0 in F-dist."},  cdf(0.0,1,-1)).

invcdf_test() ->
        ?assertEqual(0.0, invcdf(0.0,1,1)),
        ?assertEqual(0.025085630936916625, invcdf(0.1,1,1)),
        ?assertEqual(0.1172839506172838, invcdf(0.1,2,1)),
        ?assertEqual(0.9999999999999966, invcdf(0.5,1,1)),
        ?assertEqual(1.5000000000000018, invcdf(0.5,2,1)),
        ?assertEqual(60.19498034404559, invcdf(0.9,10,1)).

invcdf_error_test() ->
        ?assertEqual({error,"Nu1 parameter =< 0 in F-dist."},  invcdf(0.0,-1,1)),
        ?assertEqual({error,"Nu2 parameter =< 0 in F-dist."},  invcdf(0.0,1,-1)),
        ?assertEqual({error,"Invalid probability"},  invcdf(-0.1,1,1)).

-endif.




