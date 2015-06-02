% ==============================================================================
%  Cauchy distribution
% ==============================================================================
-module(cauchy).

-export([cauchypdf/3, cauchycdf/3, cauchyinv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(SQR(X), (X*X)).

% ------------------------------------------------------------------------------
%  cauchypdf - Cauchy probability density function
% ------------------------------------------------------------------------------
cauchypdf(_,_,Sig) when Sig =< 0 -> 
        {error,"Sigma parameter =< 0 in Cauchy dist."};

cauchypdf(X,Mu,Sig) ->
        0.318309886183790671/(Sig*(1+?SQR((X-Mu)/Sig))).

% ------------------------------------------------------------------------------
%  cauchycdf - Cauchy cumulative distribution function
% ------------------------------------------------------------------------------
cauchycdf(_,_,Sig) when Sig =< 0 -> 
        {error,"Sigma parameter =< 0 in Cauchy dist."};

cauchycdf(X,Mu,Sig) ->
        0.5+0.318309886183790671*math:atan((X-Mu)/Sig).

        
% ------------------------------------------------------------------------------
%  cauchyinv - Inverse cauchy distribution function
% ------------------------------------------------------------------------------
cauchyinv(_,_,Sig) when Sig =< 0 ->
        {error,"Sigma parameter =< 0 in Cauchy dist."};

cauchyinv(P,_,_) when P < 0 orelse P > 1 ->
        {error,"Invalid probability in Cauchy dist"};

cauchyinv(P,Mu,Sig) ->
        Mu + Sig*math:tan(3.14159265358979324*(P-0.5)).



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).


cauchypdf_test() ->
        ?assertEqual(0.07957747154594767, cauchypdf(0.0,2,2)),
        ?assertEqual(0.2136307960965038, cauchypdf(2.3,3,1)),
        ?assertEqual(0.15915494309189535, cauchypdf(20,20,2)),
        ?assertEqual(0.12732395447351627, cauchypdf(1.0,2,2)).

cauchypdf_error_test() ->
        ?assertEqual({error,"Sigma parameter =< 0 in Cauchy dist."},  cauchypdf(0.0, 1,-1)).

cauchycdf_test() ->
        ?assertEqual(0.25, cauchycdf(0.0,1,1)),
        ?assertEqual(0.5, cauchycdf(100,100,1)),
        ?assertEqual(0.017665722888134616, cauchycdf(2,20,1)),
        ?assertEqual(0.9823342771118654, cauchycdf(20,2,1)),
        ?assertEqual(0.9371670418109989, cauchycdf(-15,-20,1)),
        ?assertEqual(0.8975868006538049, cauchycdf(5.0001,2,1)).

cauchycdf_error_test() ->
        ?assertEqual({error,"Sigma parameter =< 0 in Cauchy dist."},  cauchycdf(0.0, 1,-1)).

cauchyinv_test() ->
        ?assertEqual(-2.0776835371752527, cauchyinv(0.1,1,1)),
        ?assertEqual(19.675080303767093, cauchyinv(0.4,20,1)),
        ?assertEqual(34.49196962329062, cauchyinv(0.6,2,100)),
        ?assertEqual(11.155367074350504, cauchyinv(0.9,5,2)).

-endif.



