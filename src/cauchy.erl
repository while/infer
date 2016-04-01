% ==============================================================================
%  Cauchy distribution
% ==============================================================================
-module(cauchy).

-export([pdf/3, cdf/3, invcdf/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(SQR(X), ((X)*(X))).

% ------------------------------------------------------------------------------
%  pdf - Cauchy probability density function
% ------------------------------------------------------------------------------
pdf(_,_,Sig) when Sig =< 0 -> 
        {error,"Sigma parameter =< 0 in Cauchy dist."};

pdf(X,Mu,Sig) ->
        0.318309886183790671/(Sig*(1+?SQR((X-Mu)/Sig))).

% ------------------------------------------------------------------------------
%  cdf - Cauchy cumulative distribution function
% ------------------------------------------------------------------------------
cdf(_,_,Sig) when Sig =< 0 -> 
        {error,"Sigma parameter =< 0 in Cauchy dist."};

cdf(X,Mu,Sig) ->
        0.5+0.318309886183790671*math:atan((X-Mu)/Sig).

        
% ------------------------------------------------------------------------------
%  invcdf - Inverse cauchy distribution function
% ------------------------------------------------------------------------------
invcdf(_,_,Sig) when Sig =< 0 ->
        {error,"Sigma parameter =< 0 in Cauchy dist."};

invcdf(P,_,_) when P < 0 orelse P > 1 ->
        {error,"Invalid probability in Cauchy dist"};

invcdf(P,Mu,Sig) ->
        Mu + Sig*math:tan(3.14159265358979324*(P-0.5)).



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).


pdf_test() ->
        ?assertEqual(0.07957747154594767, pdf(0.0,2,2)),
        ?assertEqual(0.2136307960965038, pdf(2.3,3,1)),
        ?assertEqual(0.15915494309189535, pdf(20,20,2)),
        ?assertEqual(0.12732395447351627, pdf(1.0,2,2)).

pdf_error_test() ->
        ?assertEqual({error,"Sigma parameter =< 0 in Cauchy dist."},  pdf(0.0, 1,-1)).

cdf_test() ->
        ?assertEqual(0.25, cdf(0.0,1,1)),
        ?assertEqual(0.5, cdf(100,100,1)),
        ?assertEqual(0.017665722888134616, cdf(2,20,1)),
        ?assertEqual(0.9823342771118654, cdf(20,2,1)),
        ?assertEqual(0.9371670418109989, cdf(-15,-20,1)),
        ?assertEqual(0.8975868006538049, cdf(5.0001,2,1)).

cdf_error_test() ->
        ?assertEqual({error,"Sigma parameter =< 0 in Cauchy dist."},  cdf(0.0, 1,-1)).

invcdf_test() ->
        ?assertEqual(-2.0776835371752527, invcdf(0.1,1,1)),
        ?assertEqual(19.675080303767093, invcdf(0.4,20,1)),
        ?assertEqual(34.49196962329062, invcdf(0.6,2,100)),
        ?assertEqual(11.155367074350504, invcdf(0.9,5,2)).

invcdf_error_test() ->
        ?assertEqual({error,"Sigma parameter =< 0 in Cauchy dist."},  invcdf(0.0, 1,-1)),
        ?assertEqual({error,"Invalid probability in Cauchy dist"},  invcdf(-0.1, 1,1)),
        ?assertEqual({error,"Invalid probability in Cauchy dist"},  invcdf(2.0, 1,1)).

-endif.



