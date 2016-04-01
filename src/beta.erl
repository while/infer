% ==============================================================================
%  Beta distribution
% ==============================================================================
-module(beta).

-export([pdf/3, cdf/3, invcdf/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(SQR(X), ((X)*(X))).

% ------------------------------------------------------------------------------
%  pdf - Probability density function
% ------------------------------------------------------------------------------
pdf(_,Alpha,_) when Alpha =< 0 -> 
        {error,"Alpha parameter =< 0 in pdf."};
pdf(_,_,Beta) when Beta =< 0 -> 
        {error,"Beta parameter =< 0 in pdf."};
pdf(X,_,_) when X =< 0.0 orelse X >= 1.0 -> 
        {error,"Bad X parameter in pdf."};

pdf(X,Alpha,Beta) ->
        Fac = incgammabeta:gammaln(Alpha+Beta) - incgammabeta:gammaln(Alpha)
              - incgammabeta:gammaln(Beta),
        math:exp((Alpha-1.0)*math:log(X)+(Beta-1.0)*math:log(1.0-X)+Fac).

% ------------------------------------------------------------------------------
%  cdf - Cumulative distribution function
% ------------------------------------------------------------------------------
cdf(_,Alpha,_) when Alpha =< 0 -> 
        {error,"Alpha parameter =< 0 in cdf."};
cdf(_,_,Beta) when Beta =< 0 -> 
        {error,"Beta parameter =< 0 in cdf."};
cdf(X,_,_) when X =< 0.0 orelse X >= 1.0 -> 
        {error,"Bad X parameter in cdf."};

cdf(X,Alpha,Beta) ->
        incgammabeta:betai(Alpha,Beta,X).

        
% ------------------------------------------------------------------------------
%  invcdf - Inverse cumulative distribution function
% ------------------------------------------------------------------------------
invcdf(_,Alpha,_) when Alpha =< 0 -> 
        {error,"Alpha parameter =< 0 in invcdf."};
invcdf(_,_,Beta) when Beta =< 0 -> 
        {error,"Beta parameter =< 0 in invcdf."};

invcdf(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};

invcdf(P,Alpha,Beta) -> 
        incgammabeta:invbetai(P,Alpha,Beta).




% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

pdf_test() ->
        ?assertEqual(3.2804999999999933, pdf(0.1,1,5)),
        ?assertEqual(3.8742048900000183, pdf(0.1,1,10)),
        ?assertEqual(1.8749999999999938, pdf(0.5,3,3)),
        ?assertEqual(1.0610329539459684, pdf(0.9,0.5,0.5)).

pdf_error_test() ->
        ?assertEqual({error,"Alpha parameter =< 0 in pdf."},  pdf(0.0,-1,1)),
        ?assertEqual({error,"Beta parameter =< 0 in pdf."},  pdf(0.0,1,-1)),
        ?assertEqual({error,"Bad X parameter in pdf."},  pdf(-1.0,1,1)),
        ?assertEqual({error,"Bad X parameter in pdf."},  pdf(2.0,1,1)).

cdf_test() ->
        ?assertEqual(0.4095099999999993, cdf(0.1,1,5)),
        ?assertEqual(0.6513215599000027, cdf(0.1,1,10)),
        ?assertEqual(0.5000000000000016, cdf(0.5,3,3)),
        ?assertEqual(0.7951672353008667, cdf(0.9,0.5,0.5)).

cdf_param_test() ->
        ?assertEqual({error,"Alpha parameter =< 0 in cdf."},  cdf(0.0,-1,1)),
        ?assertEqual({error,"Beta parameter =< 0 in cdf."},  cdf(0.0,1,-1)),
        ?assertEqual({error,"Bad X parameter in cdf."},  cdf(-1.0,1,1)),
        ?assertEqual({error,"Bad X parameter in cdf."},  cdf(2.0,1,1)).

invcdf_test() ->
        ?assertEqual(0.020851637639023268, invcdf(0.1,1,5)),
        ?assertEqual(0.010480741793785553, invcdf(0.1,1,10)),
        ?assertEqual(0.499999999999999, invcdf(0.5,3,3)),
        ?assertEqual(0.610181649469041, invcdf(0.7,3,3)),
        ?assertEqual(0.9755282581475768, invcdf(0.9,0.5,0.5)).

invcdf_param_test() ->
        ?assertEqual({error,"Alpha parameter =< 0 in invcdf."},  invcdf(0.0,-1,1)),
        ?assertEqual({error,"Beta parameter =< 0 in invcdf."},  invcdf(0.0,1,-1)),
        ?assertEqual({error,"Invalid probability"},  invcdf(-0.1,1,1)).

-endif.




