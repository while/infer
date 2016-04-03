% ==============================================================================
%  Smoothing Kernels
% ==============================================================================
-module(kernels).

-export([boxcar/1, gaussian/1, epanechnikov/1, tricube/1, triangular/1]).

-define(SQR(X), ((X)*(X))).
-define(CUBE(X), ((X)*(X)*(X))).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-define(assertEqualDigits(X,Y,D), ?assert(abs((X)-(Y)) < math:pow(10,-(D)))).
-endif.

% ------------------------------------------------------------------------------
%  Indicator function
% ------------------------------------------------------------------------------
i(X) when abs(X) =< 1 -> 1.0;
i(_) -> 0.0.

% ------------------------------------------------------------------------------
%  Boxcar Kernel
% ------------------------------------------------------------------------------
boxcar(X) ->
        0.5*i(X).


% ------------------------------------------------------------------------------
%  Gaussian Kernel
% ------------------------------------------------------------------------------
gaussian(X) ->
        0.398942280401432678*math:exp(-0.5*?SQR(X)).


% ------------------------------------------------------------------------------
%  Epanechnikov Kernel
% ------------------------------------------------------------------------------
epanechnikov(X) ->
        0.75*(1-?SQR(X))*i(X).


% ------------------------------------------------------------------------------
%  Tricube Kernel
% ------------------------------------------------------------------------------
tricube(X) ->
        0.8641975308641975*?CUBE(1-?CUBE(abs(X)))*i(X).


% ------------------------------------------------------------------------------
%  Triangular Kernel
% ------------------------------------------------------------------------------
triangular(X) ->
        (1-abs(X))*i(X).

% ------------------------------------------------------------------------------
%  Cosine Kernel
% ------------------------------------------------------------------------------
cosine(X) ->
       0.7853981633974483*math:cos(1.570796326794897*X)*i(X).


% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

indicator_test() ->
        ?assertEqual(1.0, i(1.0)),
        ?assertEqual(1.0, i(0.0)),
        ?assertEqual(1.0, i(0.5)),
        ?assertEqual(1.0, i(-0.5)),
        ?assertEqual(0.0, i(10.0)).

boxcar_test() ->
        ?assertEqual(0.5, boxcar(1.0)),
        ?assertEqual(0.5, boxcar(0.0)),
        ?assertEqual(0.5, boxcar(0.5)),
        ?assertEqual(0.5, boxcar(-0.5)),
        ?assertEqual(0.0, boxcar(10.0)),
        ?assertEqual(0.0, boxcar(-1.01)).

gaussian_test() ->
        ?assertEqual(0.24197072451914337, gaussian(1.0)),
        ?assertEqual(0.3989422804014327, gaussian(0.0)),
        ?assertEqual(0.3520653267642995, gaussian(0.5)),
        ?assertEqual(0.3520653267642995, gaussian(-0.5)),
        ?assertEqual(7.69459862670642e-23, gaussian(10.0)),
        ?assertEqual(0.23955109772801336, gaussian(-1.01)).

epanechnikov_test() ->
        ?assertEqual(0.0, epanechnikov(1.0)),
        ?assertEqual(0.75, epanechnikov(0.0)),
        ?assertEqual(0.5625, epanechnikov(0.5)),
        ?assertEqual(0.5625, epanechnikov(-0.5)),
        ?assertEqual(0.0, epanechnikov(10.0)),
        ?assertEqual(0.0, epanechnikov(-1.01)).

tricube_test() ->
        ?assertEqual(0.0, tricube(1.0)),
        ?assertEqual(70/81, tricube(0.0)),
        ?assertEqual(0.5789448302469136, tricube(0.5)),
        ?assertEqual(0.5789448302469136, tricube(-0.5)),
        ?assertEqual(0.0, tricube(10.0)),
        ?assertEqual(0.0, tricube(-1.01)).

triangular_test() ->
        ?assertEqual(0.0, triangular(1.0)),
        ?assertEqual(1.0, triangular(0.0)),
        ?assertEqual(0.5, triangular(0.5)),
        ?assertEqual(0.5, triangular(-0.5)),
        ?assertEqual(0.0, triangular(10.0)),
        ?assertEqual(0.0, triangular(-1.01)).

cosine_test() ->
        ?assertEqualDigits(0.0, cosine(1.0), 14),
        ?assertEqual(0.7853981633974483, cosine(0.0)),
        ?assertEqual(0.5553603672697957, cosine(0.5)),
        ?assertEqual(0.5553603672697957, cosine(-0.5)),
        ?assertEqual(0.0, cosine(10.0)),
        ?assertEqual(0.0, cosine(-1.01)).

-endif.
