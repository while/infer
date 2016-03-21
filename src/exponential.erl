% ==============================================================================
%  Exponential distribution
% ==============================================================================
-module(exponential).

-export([exppdf/2, expcdf/2, expinv/2, exprng/2]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


% ------------------------------------------------------------------------------
%  exppdf - Uniform probability density function
% ------------------------------------------------------------------------------
exppdf(_,Lambda) when Lambda =< 0 -> {error, "Lambda is smaller than zero."};
exppdf(X,_) when X < 0 -> 0.0;
exppdf(X,Lambda) ->
        Lambda*math:exp(-Lambda*X).

% ------------------------------------------------------------------------------
%  expcdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
expcdf(_,Lambda) when Lambda =< 0 -> {error, "Lambda is smaller than zero."};
expcdf(X,_) when X < 0 -> 0.0;
expcdf(X,Lambda) ->
        1-math:exp(-Lambda*X).
        
% ------------------------------------------------------------------------------
%  expinv - Inverse exponential distribution function
% ------------------------------------------------------------------------------
expinv(_,Lambda) when Lambda =< 0 -> {error, "Lambda is smaller than zero."};
expinv(P,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
expinv(P,Lambda) ->
        -math:log(1-P)/Lambda.


% ------------------------------------------------------------------------------
%  exprnd - RNG function
% ------------------------------------------------------------------------------
exprng(_,Lambda) when Lambda =< 0 -> {error, "Lambda is smaller than zero."};
exprng(N,Lambda) ->
        lists:map(fun(_) -> expinv(rand:uniform(),Lambda) end, lists:seq(1,N)).

% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).


% ------------------------------------------------------------------------------
%  pdf tests
% ------------------------------------------------------------------------------
exppdf_test() ->
        ?assertEqual(0.0, exppdf(-1.0,10)),
        ?assertEqual(10.0, exppdf(0.0,10)),
        ?assertEqual(1.353352832366127, exppdf(0.2,10)),
        ?assertEqual(0.06737946999085467, exppdf(0.5,10)),
        ?assertEqual(4.5399929762484856e-04, exppdf(1,10)),
        ?assertEqual(3.720075976020836e-43, exppdf(10,10)).

exppdf_error_test() ->
        ?assertEqual({error,"Lambda is smaller than zero."}, exppdf(1.0,-1)).


% ------------------------------------------------------------------------------
%  cdf tests
% ------------------------------------------------------------------------------
expcdf_test() ->
        ?assertEqual(0.0, expcdf(-1.0,10)),
        ?assertEqual(0.0, expcdf(0.0,10)),
        ?assertEqual(0.8646647167633873, expcdf(0.2,10)),
        ?assertEqual(0.9932620530009145, expcdf(0.5,10)),
        ?assertEqual(0.9990881180344455, expcdf(0.7,10)),
        ?assertEqual(0.9999546000702375, expcdf(1.0,10)),
        ?assertEqual(1.0, expcdf(10.0,10)).

expcdf_error_test() ->
        ?assertEqual({error,"Lambda is smaller than zero."}, expcdf(1.0,-1)).


% ------------------------------------------------------------------------------
%  inv tests
% ------------------------------------------------------------------------------
expinv_test() ->
        ?assertEqual(0.0, expinv(0.0,10)),
        ?assertEqual(0.010536051565782628, expinv(0.1,10)),
        ?assertEqual(0.06931471805599453, expinv(0.5,10)),
        ?assertEqual(0.1203972804325936, expinv(0.7,10)).

expinv_error_test() ->
        ?assertEqual({error,"Lambda is smaller than zero."}, expinv(1.0,-1)),
        ?assertEqual({error,"Invalid probability"}, expinv(-0.1,10)),
        ?assertEqual({error,"Invalid probability"}, expinv(1.1,10)).


% ------------------------------------------------------------------------------
%  rng tests
% ------------------------------------------------------------------------------
exprng_positive_test() ->
        [X] = exprng(1,1), 
        ?assert(X >= 0.0).

exprng_length_test() ->
        Xs = exprng(23,1), 
        ?assert(length(Xs) =:= 23).

exprng_error_test() ->
        ?assertEqual({error,"Lambda is smaller than zero."}, exprng(1,-1)).

-endif.

