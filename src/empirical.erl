% ==============================================================================
%  Empirical distribution functions
% ==============================================================================
-module(empirical).

-export([ecdf/1, ecdf/2, ecdf_ci/3]).

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


%-------------------------------------------------------------------------------
% Comfidence set for empirical cumulative distribution function
%-------------------------------------------------------------------------------
ecdf_ci(X,Xs,Alpha) ->
       Fn = ecdf(Xs),
       Eps = math:sqrt(0.5/length(Xs)*math:log(2/Alpha)),
       {max(Fn(X)-Eps, 0), min(Fn(X)+Eps, 1)}.


% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

ecdf_fun_test() ->
        Xs = [7,5,6,4,8,9,3,2,10,1],
        Fn = ecdf(Xs),
        ?assertEqual(0.0, Fn(-100)),
        ?assertEqual(0.0, Fn(0)),
        ?assertEqual(0.2, Fn(2)),
        ?assertEqual(0.4, Fn(4)),
        ?assertEqual(0.6, Fn(6)),
        ?assertEqual(0.8, Fn(8)),
        ?assertEqual(1.0, Fn(10)),
        ?assertEqual(1.0, Fn(100)).


ecdf_test() ->
        ?assertEqual(0.0, ecdf(-100, [1,2,3,5,7,6,4,8,9,10])),
        ?assertEqual(0.0, ecdf(0, [1,2,3,5,7,6,4,8,9,10])),
        ?assertEqual(0.2, ecdf(2, [1,2,3,4,5,6,7,8,9,10])),
        ?assertEqual(0.4, ecdf(4, [1,2,3,4,5,6,7,8,9,10])),
        ?assertEqual(0.6, ecdf(6, [1,2,3,4,5,6,7,8,9,10])),
        ?assertEqual(0.8, ecdf(8, [10,2,3,4,5,6,7,8,9,1])),
        ?assertEqual(1.0, ecdf(10, [1,2,3,4,5,6,7,8,9,10])),
        ?assertEqual(1.0, ecdf(100, [1,2,3,4,5,6,7,8,9,10])).


ecdf_ci_test() ->
        ?assertEqual({0,0.4294694083467376},
                     ecdf_ci(0.0, [1,2,3,4,5,6,7,8,9,10], 0.05)),
        ?assertEqual({0.07053059165326242,0.9294694083467376}, 
                     ecdf_ci(5.0, [1,2,3,4,5,6,7,8,9,10], 0.05)),
        ?assertEqual({0,0.5294694083467376},
                     ecdf_ci(1.0, [1,2,3,4,5,6,7,8,9,10], 0.05)).

-endif.
