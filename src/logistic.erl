% ==============================================================================
%  Logistic distribution
% ==============================================================================
-module(logistic).

-export([pdf/3, cdf/3, invcdf/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(SQR(X), ((X)*(X))).

% ------------------------------------------------------------------------------
%  pdf - Logistic probability density function
% ------------------------------------------------------------------------------
pdf(_,_,Sig) when Sig =< 0 -> 
        {error,"Sigma parameter =< 0 in Logistic dist."};

pdf(X,Mu,Sig) ->
        E = math:exp(-abs(1.81379936423421785*(X-Mu)/Sig)),
        1.81379936423421785*E/(Sig*?SQR(1+E)).

% ------------------------------------------------------------------------------
%  cdf - Logistic cumulative distribution function
% ------------------------------------------------------------------------------
cdf(_,_,Sig) when Sig =< 0 -> 
        {error,"Sigma parameter =< 0 in Logistic dist."};

cdf(X,Mu,Sig) ->
        E = math:exp(-1.81379936423421785*(X-Mu)/Sig),
        1/(1+E).

        
% ------------------------------------------------------------------------------
%  invcdf - Inverse logistic distribution function
% ------------------------------------------------------------------------------
invcdf(_,_,Sig) when Sig =< 0 ->
        {error,"Sigma parameter =< 0 in Logistic dist."};

invcdf(P,_,_) when P < 0 orelse P > 1 ->
        {error,"Invalid probability in Logistic dist"};

% TODO: Implement invcdf
invcdf(_,_,_) -> not_implemented.



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).


pdf_test() ->
        ?assertEqual(0.45344984105855446, pdf(0.0,0,1)),
        ?assertEqual(0.45344984105855446, pdf(100,100,1)),
        ?assertEqual(0.18582462414585452, pdf(0,1,2)),
        ?assertEqual(0.044974071984070316, pdf(1,0,10)),
        ?assertEqual(0.2186158850951135, pdf(1,0,1)).

pdf_error_test() ->
        ?assertEqual({error,"Sigma parameter =< 0 in Logistic dist."},
                     pdf(0.0, 1,-1)).

cdf_test() ->
        ?assertEqual(0.5, cdf(0.0,0,1)),
        ?assertEqual(0.5, cdf(100,100,1)),
        ?assertEqual(0.9999999999999933, cdf(20,2,1)),
        ?assertEqual(1.1516876350885883e-04, cdf(2,7,1)),
        ?assertEqual(0.9998848312364911, cdf(-15,-20,1)),
        ?assertEqual(0.673841274757147, cdf(5.0001,4.2,2)).

cdf_error_test() ->
        ?assertEqual({error,"Sigma parameter =< 0 in Logistic dist."},
                     cdf(0.0, 1,-1)).

invcdf_test() ->
        ?assertEqual(not_implemented, invcdf(0.1,1,1)).

invcdf_error_test() ->
        ?assertEqual({error,"Sigma parameter =< 0 in Logistic dist."},
                     invcdf(0.0, 1,-1)),
        ?assertEqual({error,"Invalid probability in Logistic dist"},
                     invcdf(-1.0, 1,1)).


-endif.



