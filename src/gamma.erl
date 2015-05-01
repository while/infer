% ==============================================================================
%  Gamma distribution
% ==============================================================================
-module(gamma).

-export([gammapdf/3, gammacdf/3, gammainv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


% ------------------------------------------------------------------------------
%  gammapdf - Uniform probability density function
% ------------------------------------------------------------------------------
gammapdf(X,_,_) when X < 0 -> 0.0;
gammapdf(X,Shape,Rate) ->
        not_implemented.

% ------------------------------------------------------------------------------
%  gammapcdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
gammacdf(X,_,_) when X < 0 -> 0.0;
gammacdf(X,Shape,Rate) ->
        not_implemented.

        
% ------------------------------------------------------------------------------
%  gammapinv - Inverse gamma distribution function
% ------------------------------------------------------------------------------
gammainv(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
gammainv(P,Shape,Rate) ->
        not_implemented.



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

gammapdf_test() ->
        ?assertEqual(0.0, gammapdf(0.0,1,1)).


gammacdf_test() ->
        ?assertEqual(0.0, gammacdf(0.0,1,1)),

gammainv_test() ->
        ?assertEqual(0.0, gammainv(0.0,1,1)).

-endif.



