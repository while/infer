% ==============================================================================
%  Uniform distribution
% ==============================================================================
-module(uniform).

-export([pdf/3, cdf/3, invcdf/3, rnd/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


% ------------------------------------------------------------------------------
%  pdf - Uniform probability density function
% ------------------------------------------------------------------------------
pdf(X,A,_) when X < A -> 0.0;
pdf(X,_,B) when X > B -> 0.0;
pdf(_,A,B) ->
        1/(B-A).

% ------------------------------------------------------------------------------
%  cdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
cdf(X,A,_) when X < A -> 0.0;
cdf(X,_,B) when X > B -> 1.0;
cdf(X,A,B) ->
        (X-A)/(B-A).
        
% ------------------------------------------------------------------------------
%  invcdf - Inverse uniform distribution function
% ------------------------------------------------------------------------------
invcdf(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
invcdf(P,A,B) ->
        A + P*(B-A).


% ------------------------------------------------------------------------------
%  rnd - RNG function
% ------------------------------------------------------------------------------
rnd(N,A,B) ->
        lists:map(fun(_) -> invcdf(rand:uniform(),A,B) end, lists:seq(1,N)).


% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

% ------------------------------------------------------------------------------
%  pdf tests
% ------------------------------------------------------------------------------
pdf_test() ->
        ?assertEqual(0.0, pdf(-0.999,0,10)),
        ?assertEqual(0.1, pdf(0,0,10)),
        ?assertEqual(0.1, pdf(2,0,10)),
        ?assertEqual(0.1, pdf(5,0,10)),
        ?assertEqual(0.1, pdf(7,0,10)),
        ?assertEqual(0.1, pdf(10,0,10)),
        ?assertEqual(0.0, pdf(10.001,0,10)).


% ------------------------------------------------------------------------------
%  cdf tests
% ------------------------------------------------------------------------------
cdf_test() ->
        ?assertEqual(0.0, cdf(-0.999,0,10)),
        ?assertEqual(0.0, cdf(0,0,10)),
        ?assertEqual(0.2, cdf(2,0,10)),
        ?assertEqual(0.5, cdf(5,0,10)),
        ?assertEqual(0.7, cdf(7,0,10)),
        ?assertEqual(1.0, cdf(10,0,10)),
        ?assertEqual(1.0, cdf(10.001,0,10)).


% ------------------------------------------------------------------------------
%  inv tests
% ------------------------------------------------------------------------------
invcdf_test() ->
        ?assertEqual(1.0, invcdf(0.1,0,10)),
        ?assertEqual(5.0, invcdf(0.5,0,10)),
        ?assertEqual(10.0, invcdf(1.0,0,10)).

invcdf_error_test() ->
        ?assertEqual({error,"Invalid probability"}, invcdf(-0.1,0,10)),
        ?assertEqual({error,"Invalid probability"}, invcdf(1.1,0,10)).


% ------------------------------------------------------------------------------
%  rng tests
% ------------------------------------------------------------------------------
rnd_positive_test() ->
        [X] = rnd(1,0,10), 
        ?assert(X >= 0.0 andalso X =< 10.0).

rnd_length_test() ->
        Xs = rnd(23,1,3), 
        ?assert(length(Xs) =:= 23).


-endif.
