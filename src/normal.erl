% ==============================================================================
%  Normal distribution
% ==============================================================================
-module(normal).

-export([normpdf/3, normcdf/3, norminv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

% ------------------------------------------------------------------------------
%  normpdf - Normal probability density function
% ------------------------------------------------------------------------------
normpdf(X,Mu,Sig) ->
        (0.398942280401432678/Sig)*math:exp(-0.5*math:pow((X - Mu)/Sig, 2)).

% ------------------------------------------------------------------------------
%  normcdf - Normal cumulative distribution function
% ------------------------------------------------------------------------------
normcdf(X,Mu,Sig) ->
        0.5*erfc(-0.707106781186547524*(X-Mu)/Sig).
        
% ------------------------------------------------------------------------------
%  norminv - Inverse normal distribution function
% ------------------------------------------------------------------------------
norminv(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
norminv(P,Mu,Sig) ->
        -1.41421356237309505*Sig*inverfc(2*P)+Mu.


% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

normpdf_test() ->
        ?assertEqual(0.0, normpdf(-100,0,1)),
        ?assertEqual(0.3520653267642995, normpdf(-0.5,0,1)),
        ?assertEqual(0.3989422804014327, normpdf(0,0,1)),
        ?assertEqual(0.3520653267642995, normpdf(0.5,0,1)),
        ?assertEqual(0.0, normpdf(100,0,1)).


normcdf_test() ->
        ?assertEqual(0.0, normcdf(-100,0,1)),
        ?assertEqual(0.3085375387259869, normcdf(-0.5,0,1)),
        ?assertEqual(0.5, normcdf(0,0,1)),
        ?assertEqual(0.6914624612740131, normcdf(0.5,0,1)),
        ?assertEqual(1.0, normcdf(100,0,1)).

norminv_test() ->
        ?assert(1.0e-16 >= norminv(0.5,0,1)).

-endif.
