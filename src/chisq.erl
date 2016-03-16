% ==============================================================================
%  Chi-square distribution
% ==============================================================================
-module(chisq).

-export([chisqpdf/2, chisqcdf/2, chisqinv/2]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

% ------------------------------------------------------------------------------
%  chisqpdf - Uniform probability density function
% ------------------------------------------------------------------------------
chisqpdf(_,Nu) when Nu =< 0 -> {error,"Nu parameter =< 0 in Chi-square."};

chisqpdf(X,_) when X < 0 -> 0.0;
chisqpdf(X,Nu) ->
        Fac = 0.693147180559945309*(0.5*Nu) + incgammabeta:gammaln(0.5*Nu),
        math:exp(-0.5*(X-(Nu-2)*math:log(X))-Fac).

% ------------------------------------------------------------------------------
%  chisqcdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
chisqcdf(_,Nu) when Nu =< 0 -> {error,"Nu parameter =< 0 in Chi-square."};

chisqcdf(X,_) when X < 0 -> 0.0;
chisqcdf(X,Nu) ->
        incgammabeta:gammap(0.5*Nu,0.5*X).

        
% ------------------------------------------------------------------------------
%  chisqinv - Inverse chisq distribution function
% ------------------------------------------------------------------------------
chisqinv(_,Nu) when Nu =< 0 -> {error,"Nu parameter =< 0 in Chi-square."};

chisqinv(P,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};

chisqinv(X,_) when X =< 0 -> 0.0;
chisqinv(P,Nu) ->
        2.0*incgammabeta:invgammap(P,0.5*Nu).


% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

chisqpdf_test() ->
        ?assertEqual(0.0, chisqpdf(-1.0,1)),
        ?assertEqual(0.24197072451914325, chisqpdf(1,1)),
        ?assertEqual(0.10377687435514864, chisqpdf(2,1)),
        ?assertEqual(0.3032653298563167, chisqpdf(1,2)),
        ?assertEqual(0.003368973499542733, chisqpdf(10,2)),
        ?assertEqual(0.00766415502440505, chisqpdf(2,10)).

chisqpdf_error_test() ->
        ?assertEqual({error,"Nu parameter =< 0 in Chi-square."},  chisqpdf(0.0,-1)).

chisqcdf_test() ->
        ?assertEqual(0.0, chisqcdf(-1.0,1)),
        ?assertEqual(0.6826894921370856, chisqcdf(1,1)),
        ?assertEqual(0.8427007929497147, chisqcdf(2,1)),
        ?assertEqual(0.3934693402873665, chisqcdf(1,2)),
        ?assertEqual(0.9932620530009145, chisqcdf(10,2)),
        ?assertEqual(0.0036598468273437135, chisqcdf(2,10)),
        ?assertEqual(0.5132987982791589, chisqcdf(200,200)).

chisqcdf_error_test() ->
        ?assertEqual({error,"Nu parameter =< 0 in Chi-square."},  chisqcdf(0.0,-1)).

chisqinv_test() ->
        ?assertEqual(0.0, chisqinv(0.0,1)),
        ?assertEqual(0.015790774093431218, chisqinv(0.1,1)),
        ?assertEqual(0.2107210313156527, chisqinv(0.1,2)),
        ?assertEqual(0.45493642311957305, chisqinv(0.5,1)),
        ?assertEqual(1.3862943611198906, chisqinv(0.5,2)),
        ?assertEqual(15.987179172105261, chisqinv(0.9,10)).

chisqinv_error_test() ->
        ?assertEqual({error,"Nu parameter =< 0 in Chi-square."},  chisqinv(0.5,-1)),
        ?assertEqual({error,"Invalid probability"},  chisqinv(-0.5,1)),
        ?assertEqual({error,"Invalid probability"},  chisqinv(2,1)).

-endif.




