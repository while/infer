% ==============================================================================
%  Empirical distribution functions
% ==============================================================================
-module(empirical).

-export([ecdf/1, ecdf/2]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


%-------------------------------------------------------------------------------
% Empirical cumulative distribution function
%-------------------------------------------------------------------------------
ecdf(Xs) ->
        fun(X) -> ecdf(X,Xs) end.

ecdf(X,Xs) ->
       ecdf_(X,lists:sort(Xs),0.0)/length(Xs). 

ecdf_(_,[],Idx) -> Idx;
ecdf_(X,[H|_],Idx) when X < H -> Idx;
ecdf_(X,[_|T],Idx) -> ecdf_(X,T,Idx+1).


% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

ecdf_fun_test() ->
        Xs = [7,5,6,4,8,9,3,2,10,1],
        F = ecdf(Xs),
        ?assertEqual(0.0, F(-100)),
        ?assertEqual(0.0, F(0)),
        ?assertEqual(0.2, F(2)),
        ?assertEqual(0.4, F(4)),
        ?assertEqual(0.6, F(6)),
        ?assertEqual(0.8, F(8)),
        ?assertEqual(1.0, F(10)),
        ?assertEqual(1.0, F(100)).

ecdf_test() ->
        ?assertEqual(0.0, ecdf(-100, [1,2,3,5,7,6,4,8,9,10])),
        ?assertEqual(0.0, ecdf(0, [1,2,3,5,7,6,4,8,9,10])),
        ?assertEqual(0.2, ecdf(2, [1,2,3,4,5,6,7,8,9,10])),
        ?assertEqual(0.4, ecdf(4, [1,2,3,4,5,6,7,8,9,10])),
        ?assertEqual(0.6, ecdf(6, [1,2,3,4,5,6,7,8,9,10])),
        ?assertEqual(0.8, ecdf(8, [10,2,3,4,5,6,7,8,9,1])),
        ?assertEqual(1.0, ecdf(10, [1,2,3,4,5,6,7,8,9,10])),
        ?assertEqual(1.0, ecdf(100, [1,2,3,4,5,6,7,8,9,10])).




-endif.
