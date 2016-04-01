% ==============================================================================
%  Chi-square distribution
% ==============================================================================
-module(chisq).

-export([pdf/2, cdf/2, invcdf/2]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

% ------------------------------------------------------------------------------
%  pdf - Uniform probability density function
% ------------------------------------------------------------------------------
pdf(_,Nu) when Nu =< 0 -> {error,"Nu parameter =< 0 in Chi-square."};

pdf(X,_) when X < 0 -> 0.0;
pdf(X,Nu) ->
        Fac = 0.693147180559945309*(0.5*Nu) + incgammabeta:gammaln(0.5*Nu),
        math:exp(-0.5*(X-(Nu-2)*math:log(X))-Fac).

% ------------------------------------------------------------------------------
%  cdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
cdf(_,Nu) when Nu =< 0 -> {error,"Nu parameter =< 0 in Chi-square."};

cdf(X,_) when X < 0 -> 0.0;
cdf(X,Nu) ->
        incgammabeta:gammap(0.5*Nu,0.5*X).

        
% ------------------------------------------------------------------------------
%  invcdf - Inverse chisq distribution function
% ------------------------------------------------------------------------------
invcdf(_,Nu) when Nu =< 0 -> {error,"Nu parameter =< 0 in Chi-square."};

invcdf(P,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};

invcdf(X,_) when X =< 0 -> 0.0;
invcdf(P,Nu) ->
        2.0*incgammabeta:invgammap(P,0.5*Nu).


% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

pdf_test() ->
        ?assertEqual(0.0, pdf(-1.0,1)),
        ?assertEqual(0.24197072451914325, pdf(1,1)),
        ?assertEqual(0.10377687435514864, pdf(2,1)),
        ?assertEqual(0.3032653298563167, pdf(1,2)),
        ?assertEqual(0.003368973499542733, pdf(10,2)),
        ?assertEqual(0.00766415502440505, pdf(2,10)).

pdf_error_test() ->
        ?assertEqual({error,"Nu parameter =< 0 in Chi-square."},  pdf(0.0,-1)).

cdf_test() ->
        ?assertEqual(0.0, cdf(-1.0,1)),
        ?assertEqual(0.6826894921370856, cdf(1,1)),
        ?assertEqual(0.8427007929497147, cdf(2,1)),
        ?assertEqual(0.3934693402873665, cdf(1,2)),
        ?assertEqual(0.9932620530009145, cdf(10,2)),
        ?assertEqual(0.0036598468273437135, cdf(2,10)),
        ?assertEqual(0.5132987982791589, cdf(200,200)).

cdf_error_test() ->
        ?assertEqual({error,"Nu parameter =< 0 in Chi-square."},  cdf(0.0,-1)).

invcdf_test() ->
        ?assertEqual(0.0, invcdf(0.0,1)),
        ?assertEqual(0.015790774093431218, invcdf(0.1,1)),
        ?assertEqual(0.2107210313156527, invcdf(0.1,2)),
        ?assertEqual(0.45493642311957305, invcdf(0.5,1)),
        ?assertEqual(1.3862943611198906, invcdf(0.5,2)),
        ?assertEqual(15.987179172105261, invcdf(0.9,10)).

invcdf_error_test() ->
        ?assertEqual({error,"Nu parameter =< 0 in Chi-square."},  invcdf(0.5,-1)),
        ?assertEqual({error,"Invalid probability"},  invcdf(-0.5,1)),
        ?assertEqual({error,"Invalid probability"},  invcdf(2,1)).

-endif.




