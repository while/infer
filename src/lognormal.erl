% ==============================================================================
%  Log-Normal distribution
% ==============================================================================
-module(lognormal).

-export([lognpdf/3, logncdf/3, logninv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


% ------------------------------------------------------------------------------
%  lognpdf - Log-Normal probability density function
% ------------------------------------------------------------------------------
lognpdf(X,_,_) when X =< 0  -> 0.0;
lognpdf(X,Mu,Sig) ->
        (0.398942280401432678/(Sig*X))*math:exp(-0.5*math:pow((math:log(X) - Mu)/Sig, 2)).

% ------------------------------------------------------------------------------
%  logncdf - Log-Normal cumulative distribution function
% ------------------------------------------------------------------------------
logncdf(X,_,_) when X =< 0  -> 0.0;
logncdf(X,Mu,Sig) ->
        0.5*erfc(-0.707106781186547524*(math:log(X)-Mu)/Sig).
        
% ------------------------------------------------------------------------------
%  logninv - Inverse loglognal distribution function
% ------------------------------------------------------------------------------
logninv(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
logninv(P,Mu,Sig) ->
        math:exp(-1.41421356237309505*Sig*inverfc(2.0*P)+Mu).

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


% ==============================================================================
%  Tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

lognpdf_test() ->
        ?assertEqual(0.0, lognpdf(-1.0,0,1)),
        ?assertEqual(0.0, lognpdf(0.0,0,1)),
        ?assertEqual(0.28159018901526833, lognpdf(0.1,0,1)),
        ?assertEqual(0.6274960771159244, lognpdf(0.5,0,1)),
        ?assertEqual(0.5347948320769198, lognpdf(0.7,0,1)),
        ?assertEqual(0.3989422804014327, lognpdf(1.0,0,1)),
        ?assertEqual(0.0028159018901526794, lognpdf(10.0,0,1)).


logncdf_test() ->
        ?assertEqual(0.0, logncdf(-1.0,0,1)),
        ?assertEqual(0.0, logncdf(0.0,0,1)),
        ?assertEqual(0.01065109934170011, logncdf(0.1,0,1)),
        ?assertEqual(0.24410859578558275, logncdf(0.5,0,1)),
        ?assertEqual(0.3606675826226491, logncdf(0.7,0,1)),
        ?assertEqual(0.5, logncdf(1.0,0,1)),
        ?assertEqual(0.9893489006582998, logncdf(10.0,0,1)).

logninv_test() ->
        ?assertEqual(0.27760624185200977, logninv(0.1,0,1)),
        ?assertEqual(1.0, logninv(0.5,0,1)),
        ?assertEqual(1.6894457434840042, logninv(0.7,0,1)).

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

logninv_error_test() ->
        ?assertEqual({error,"Invalid probability"}, logninv(-0.1,0,10)),
        ?assertEqual({error,"Invalid probability"}, logninv(1.1,0,10)).

-endif.
