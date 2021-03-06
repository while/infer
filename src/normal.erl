% ==============================================================================
%  Normal distribution
% ==============================================================================
-module(normal).

-export([pdf/3, cdf/3, invcdf/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


% ------------------------------------------------------------------------------
%  pdf - Normal probability density function
% ------------------------------------------------------------------------------
pdf(X,Mu,Sig) ->
        (0.398942280401432678/Sig)*math:exp(-0.5*math:pow((X - Mu)/Sig, 2)).

% ------------------------------------------------------------------------------
%  cdf - Normal cumulative distribution function
% ------------------------------------------------------------------------------
cdf(X,Mu,Sig) ->
        0.5*erfc(-0.707106781186547524*(X-Mu)/Sig).
        
% ------------------------------------------------------------------------------
%  invcdf - Inverse normal distribution function
% ------------------------------------------------------------------------------
invcdf(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
invcdf(P,Mu,Sig) ->
        -1.41421356237309505*Sig*inverfc(2*P)+Mu.

% ------------------------------------------------------------------------------
%  Complementary error function 
% ------------------------------------------------------------------------------
erfc(X) -> 
        1 - math:erf(X).

% ------------------------------------------------------------------------------
%  Inverse complementary error function
% ------------------------------------------------------------------------------
inverfc(P) when P >= 2 -> -100;
inverfc(P) when P =< 0 ->  100;
inverfc(P) ->
        PP = if P <  1 -> P
              ; P >= 1 -> 2 - P
             end,
        T = math:sqrt(-2*math:log(PP/2)),
        X0 = -0.70711*((2.30753+T*0.27061)/(1+T*(0.99229+T*0.04481)) - T),
        Err1 = erfc(X0) - PP,
        X1 = X0 + Err1/(1.12837916709551257*math:exp(-X0*X0)-X0*Err1),
        Err2 = erfc(X1) - PP,
        X2 = X1 + Err2/(1.12837916709551257*math:exp(-X1*X1)-X1*Err2),
        if P <  1 -> X2
         ; P >= 1 -> -X2
        end.


% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
%  Tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

pdf_test() ->
        ?assertEqual(0.0, pdf(-100,0,1)),
        ?assertEqual(0.3520653267642995, pdf(-0.5,0,1)),
        ?assertEqual(0.3989422804014327, pdf(0,0,1)),
        ?assertEqual(0.3520653267642995, pdf(0.5,0,1)),
        ?assertEqual(0.0, pdf(100,0,1)).


cdf_test() ->
        ?assertEqual(0.0, cdf(-100,0,1)),
        ?assertEqual(0.3085375387259869, cdf(-0.5,0,1)),
        ?assertEqual(0.5, cdf(0,0,1)),
        ?assertEqual(0.6914624612740131, cdf(0.5,0,1)),
        ?assertEqual(1.0, cdf(100,0,1)).

invcdf_test() ->
        ?assert(1.0e-16 >= invcdf(0.5,0,1)).

erfc_test() ->
        ?assertEqual(0.0, erfc(100)),
        ?assertEqual(1.0, erfc(0)),
        ?assertEqual(2.0, erfc(-100)).


inverfc_test() ->
        ?assertEqual(100, inverfc(0)),
        ?assertEqual(0.3708071585935579, inverfc(0.6)),
        ?assert(1.0e-16 >= inverfc(1)),
        ?assertEqual(-0.3708071585935579, inverfc(1.4)),
        ?assertEqual(-100, inverfc(2)).

invcdf_error_test() ->
        ?assertEqual({error,"Invalid probability"}, invcdf(-0.1,0,10)),
        ?assertEqual({error,"Invalid probability"}, invcdf(1.1,0,10)).

-endif.
