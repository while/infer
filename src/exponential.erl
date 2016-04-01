% ==============================================================================
%  Exponential distribution
% ==============================================================================
-module(exponential).

-export([pdf/2, cdf/2, invcdf/2, rnd/2]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


% ------------------------------------------------------------------------------
%  pdf - Uniform probability density function
% ------------------------------------------------------------------------------
pdf(_,Lambda) when Lambda =< 0 -> {error, "Lambda is smaller than zero."};
pdf(X,_) when X < 0 -> 0.0;
pdf(X,Lambda) ->
        Lambda*math:exp(-Lambda*X).

% ------------------------------------------------------------------------------
%  cdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
cdf(_,Lambda) when Lambda =< 0 -> {error, "Lambda is smaller than zero."};
cdf(X,_) when X < 0 -> 0.0;
cdf(X,Lambda) ->
        1-math:exp(-Lambda*X).
        
% ------------------------------------------------------------------------------
%  invcdf - Inverse exponential distribution function
% ------------------------------------------------------------------------------
invcdf(_,Lambda) when Lambda =< 0 -> {error, "Lambda is smaller than zero."};
invcdf(P,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
invcdf(P,Lambda) ->
        -math:log(1-P)/Lambda.


% ------------------------------------------------------------------------------
%  exprnd - RNG function
% ------------------------------------------------------------------------------
rnd(_,Lambda) when Lambda =< 0 -> {error, "Lambda is smaller than zero."};
rnd(N,Lambda) ->
        lists:map(fun(_) -> invcdf(rand:uniform(),Lambda) end, lists:seq(1,N)).

% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).


% ------------------------------------------------------------------------------
%  pdf tests
% ------------------------------------------------------------------------------
pdf_test() ->
        ?assertEqual(0.0, pdf(-1.0,10)),
        ?assertEqual(10.0, pdf(0.0,10)),
        ?assertEqual(1.353352832366127, pdf(0.2,10)),
        ?assertEqual(0.06737946999085467, pdf(0.5,10)),
        ?assertEqual(4.5399929762484856e-04, pdf(1,10)),
        ?assertEqual(3.720075976020836e-43, pdf(10,10)).

pdf_error_test() ->
        ?assertEqual({error,"Lambda is smaller than zero."}, pdf(1.0,-1)).


% ------------------------------------------------------------------------------
%  cdf tests
% ------------------------------------------------------------------------------
cdf_test() ->
        ?assertEqual(0.0, cdf(-1.0,10)),
        ?assertEqual(0.0, cdf(0.0,10)),
        ?assertEqual(0.8646647167633873, cdf(0.2,10)),
        ?assertEqual(0.9932620530009145, cdf(0.5,10)),
        ?assertEqual(0.9990881180344455, cdf(0.7,10)),
        ?assertEqual(0.9999546000702375, cdf(1.0,10)),
        ?assertEqual(1.0, cdf(10.0,10)).

cdf_error_test() ->
        ?assertEqual({error,"Lambda is smaller than zero."}, cdf(1.0,-1)).


% ------------------------------------------------------------------------------
%  inv tests
% ------------------------------------------------------------------------------
invcdf_test() ->
        ?assertEqual(0.0, invcdf(0.0,10)),
        ?assertEqual(0.010536051565782628, invcdf(0.1,10)),
        ?assertEqual(0.06931471805599453, invcdf(0.5,10)),
        ?assertEqual(0.1203972804325936, invcdf(0.7,10)).

invcdf_error_test() ->
        ?assertEqual({error,"Lambda is smaller than zero."}, invcdf(1.0,-1)),
        ?assertEqual({error,"Invalid probability"}, invcdf(-0.1,10)),
        ?assertEqual({error,"Invalid probability"}, invcdf(1.1,10)).


% ------------------------------------------------------------------------------
%  rng tests
% ------------------------------------------------------------------------------
rnd_positive_test() ->
        [X] = rnd(1,1), 
        ?assert(X >= 0.0).

rnd_length_test() ->
        Xs = rnd(23,1), 
        ?assert(length(Xs) =:= 23).

rnd_error_test() ->
        ?assertEqual({error,"Lambda is smaller than zero."}, rnd(1,-1)).

-endif.

