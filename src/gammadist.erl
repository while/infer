% ==============================================================================
%  Gamma distribution
% ==============================================================================
-module(gammadist).

-export([pdf/3, cdf/3, invcdf/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

% ------------------------------------------------------------------------------
%  pdf - Uniform probability density function
% ------------------------------------------------------------------------------
pdf(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
pdf(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

pdf(X,_,_) when X =< 0 -> 0.0;
pdf(X,Alpha,Beta) ->
        Fac = Alpha*math:log(Beta)-incgammabeta:gammaln(Alpha),
        math:exp(-Beta*X+(Alpha-1)*math:log(X)+Fac).

% ------------------------------------------------------------------------------
%  cdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
cdf(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
cdf(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

cdf(X,_,_) when X =< 0 -> 0.0;
cdf(X,Alpha,Beta) ->
        incgammabeta:gammap(Alpha,Beta*X).

        
% ------------------------------------------------------------------------------
%  invcdf - Inverse gamma distribution function
% ------------------------------------------------------------------------------
invcdf(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
invcdf(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

invcdf(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
invcdf(P,Alpha,Beta) ->
        incgammabeta:invgammap(P,Alpha)/Beta.



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

pdf_test() ->
        ?assertEqual(0.0, pdf(0.0,2,2)),
        ?assertEqual(0.2651846416468157, pdf(2.3,3,1)),
        ?assertEqual(1.9199765904692206e-04, pdf(20,20,2)),
        ?assertEqual(0.5413411329464507, pdf(1.0,2,2)).

pdf_error_test() ->
        ?assertEqual({error,"Alpha parameter =< 0."}, pdf(0.0,-1, 1)),
        ?assertEqual({error,"Beta parameter =< 0."},  pdf(0.0, 1,-1)).

cdf_test() ->
        ?assertEqual(0.0, cdf(0.0,1,1)),
        ?assertEqual(0.5132987982791589, cdf(100,100,1)),
        ?assertEqual(0.5297427331607533, cdf(20,20,1)),
        ?assertEqual(0.9998236971022615, cdf(20,20,2)),
        ?assertEqual(0.12478121503252411, cdf(15,20,1)),
        ?assertEqual(0.9595723180054873, cdf(5,2,1)).

cdf_error_test() ->
        ?assertEqual({error,"Alpha parameter =< 0."}, cdf(0.0,-1, 1)),
        ?assertEqual({error,"Beta parameter =< 0."},  cdf(0.0, 1,-1)).

invcdf_test() ->
        ?assertEqual(0.0, invcdf(0.0,1,1)),
        ?assertEqual(0.6931471805599453, invcdf(0.5,1,1)),
        ?assertEqual(0.22746821155978653, invcdf(0.5,0.5,1)),
        ?assertEqual(4.670908882795983, invcdf(0.5,5,1)),
        ?assertEqual(2.432591025962664, invcdf(0.1,5,1)),
        ?assertEqual(4.865182051925328, invcdf(0.1,5,1/2)),
        ?assertEqual(1.216295512981332, invcdf(0.1,5,2)),
        ?assertEqual(3.9967947930263152, invcdf(0.9,5,2)).

invcdf_error_test() ->
        ?assertEqual({error,"Alpha parameter =< 0."}, invcdf(0.0,-1, 1)),
        ?assertEqual({error,"Beta parameter =< 0."},  invcdf(0.0, 1,-1)),
        ?assertEqual({error,"Invalid probability"},  invcdf(-1.0,1,1)),
        ?assertEqual({error,"Invalid probability"},  invcdf(1.1,1,1)).

-endif.



