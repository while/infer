% ==============================================================================
%  Log-Normal distribution
% ==============================================================================
-module(lognormal).

-export([lognpdf/3, logncdf/3, logninv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-include("src/erf.erl").

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


% ==============================================================================
%  EUnit tests
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

-endif.
