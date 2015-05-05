% ==============================================================================
%  Uniform distribution
% ==============================================================================
-module(uniform).

-export([unifpdf/3, unifcdf/3, unifinv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


% ------------------------------------------------------------------------------
%  unifpdf - Uniform probability density function
% ------------------------------------------------------------------------------
unifpdf(X,A,_) when X < A -> 0.0;
unifpdf(X,_,B) when X > B -> 0.0;
unifpdf(_,A,B) ->
        1/(B-A).

% ------------------------------------------------------------------------------
%  unifcdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
unifcdf(X,A,_) when X < A -> 0.0;
unifcdf(X,_,B) when X > B -> 1.0;
unifcdf(X,A,B) ->
        (X-A)/(B-A).
        
% ------------------------------------------------------------------------------
%  unifinv - Inverse uniform distribution function
% ------------------------------------------------------------------------------
unifinv(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
unifinv(P,A,B) ->
        A + P*(B-A).



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

unifpdf_test() ->
        ?assertEqual(0.0, unifpdf(-0.999,0,10)),
        ?assertEqual(0.1, unifpdf(0,0,10)),
        ?assertEqual(0.1, unifpdf(2,0,10)),
        ?assertEqual(0.1, unifpdf(5,0,10)),
        ?assertEqual(0.1, unifpdf(7,0,10)),
        ?assertEqual(0.1, unifpdf(10,0,10)),
        ?assertEqual(0.0, unifpdf(10.001,0,10)).

unifcdf_test() ->
        ?assertEqual(0.0, unifcdf(-0.999,0,10)),
        ?assertEqual(0.0, unifcdf(0,0,10)),
        ?assertEqual(0.2, unifcdf(2,0,10)),
        ?assertEqual(0.5, unifcdf(5,0,10)),
        ?assertEqual(0.7, unifcdf(7,0,10)),
        ?assertEqual(1.0, unifcdf(10,0,10)),
        ?assertEqual(1.0, unifcdf(10.001,0,10)).

unifinv_test() ->
        ?assertEqual(1.0, unifinv(0.1,0,10)),
        ?assertEqual(5.0, unifinv(0.5,0,10)),
        ?assertEqual(10.0, unifinv(1.0,0,10)).

unifinv_error_test() ->
        ?assertEqual({error,"Invalid probability"}, unifinv(-0.1,0,10)),
        ?assertEqual({error,"Invalid probability"}, unifinv(1.1,0,10)).

-endif.
