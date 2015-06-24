% ==============================================================================
%  Logistic distribution
% ==============================================================================
-module(logistic).

-export([logisticpdf/3, logisticcdf/3, logisticinv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(SQR(X), ((X)*(X))).

% ------------------------------------------------------------------------------
%  logisticpdf - Logistic probability density function
% ------------------------------------------------------------------------------
logisticpdf(_,_,Sig) when Sig =< 0 -> 
        {error,"Sigma parameter =< 0 in Logistic dist."};

logisticpdf(X,Mu,Sig) ->
        E = math:exp(-abs(1.81379936423421785*(X-Mu)/Sig)),
        1.81379936423421785*E/(Sig*?SQR(1+E)).

% ------------------------------------------------------------------------------
%  logisticcdf - Logistic cumulative distribution function
% ------------------------------------------------------------------------------
logisticcdf(_,_,Sig) when Sig =< 0 -> 
        {error,"Sigma parameter =< 0 in Logistic dist."};

logisticcdf(X,Mu,Sig) ->
        E = math:exp(-1.81379936423421785*(X-Mu)/Sig),
        1/(1+E).

        
% ------------------------------------------------------------------------------
%  logisticinv - Inverse logistic distribution function
% ------------------------------------------------------------------------------
logisticinv(_,_,Sig) when Sig =< 0 ->
        {error,"Sigma parameter =< 0 in Logistic dist."};

logisticinv(P,_,_) when P < 0 orelse P > 1 ->
        {error,"Invalid probability in Logistic dist"};

logisticinv(_,_,_) -> not_implemented.



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).


logisticpdf_test() ->
        ?assertEqual(0.45344984105855446, logisticpdf(0.0,0,1)),
        ?assertEqual(0.45344984105855446, logisticpdf(100,100,1)),
        ?assertEqual(0.18582462414585452, logisticpdf(0,1,2)),
        ?assertEqual(0.044974071984070316, logisticpdf(1,0,10)),
        ?assertEqual(0.2186158850951135, logisticpdf(1,0,1)).

logisticpdf_error_test() ->
        ?assertEqual({error,"Sigma parameter =< 0 in Logistic dist."},  logisticpdf(0.0, 1,-1)).

logisticcdf_test() ->
        ?assertEqual(0.5, logisticcdf(0.0,0,1)),
        ?assertEqual(0.5, logisticcdf(100,100,1)),
        ?assertEqual(0.9999999999999933, logisticcdf(20,2,1)),
        ?assertEqual(1.1516876350885883e-04, logisticcdf(2,7,1)),
        ?assertEqual(0.9998848312364911, logisticcdf(-15,-20,1)),
        ?assertEqual(0.673841274757147, logisticcdf(5.0001,4.2,2)).

logisticcdf_error_test() ->
        ?assertEqual({error,"Sigma parameter =< 0 in Logistic dist."},  logisticcdf(0.0, 1,-1)).

logisticinv_test() ->
        ?assertEqual(not_implemented, logisticinv(0.1,1,1)).

-endif.



