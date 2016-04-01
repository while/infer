% ==============================================================================
%  Weibull distribution
% ==============================================================================
-module(weibull).

-export([pdf/3, cdf/3, invcdf/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


% ------------------------------------------------------------------------------
%  pdf - Uniform probability density function
% ------------------------------------------------------------------------------
pdf(X,_,_) when X < 0 -> 0.0;
pdf(X,A,B) ->
        (A/B)*math:pow(X/B,A-1)*math:exp(-math:pow(X/B,A)).

% ------------------------------------------------------------------------------
%  weibullpcdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
cdf(X,_,_) when X < 0 -> 0.0;
cdf(X,A,B) ->
        1-math:exp(-math:pow(X/B,A)).

        
% ------------------------------------------------------------------------------
%  weibullpinv - Inverse weibull distribution function
% ------------------------------------------------------------------------------
invcdf(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
invcdf(P,A,B) ->
        B*math:pow(-math:log(1-P),1/A).



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

pdf_test() ->
        ?assertEqual(0.0, pdf(-1.0,2,1)),
        ?assertEqual(0.0, pdf(0.0,2,1)),
        ?assertEqual(0.7788007830714049, pdf(0.5,2,1)),
        ?assertEqual(0.7357588823428847, pdf(1.0,2,1)),
        ?assertEqual(0.3161976736855930, pdf(1.5,2,1)),
        ?assertEqual(0.07326255555493671, pdf(2.0,2,1)),
        ?assertEqual(0.009652270681138546, pdf(2.5,2,1)),
        ?assertEqual(0.0007404588245200774, pdf(3.0,2,1)).

cdf_test() ->
        ?assertEqual(0.0, cdf(-1.0,2,1)),
        ?assertEqual(0.0, cdf(0.0,2,1)),
        ?assertEqual(0.22119921692859512, cdf(0.5,2,1)),
        ?assertEqual(0.6321205588285577, cdf(1.0,2,1)),
        ?assertEqual(0.8946007754381357, cdf(1.5,2,1)),
        ?assertEqual(0.9816843611112658, cdf(2.0,2,1)),
        ?assertEqual(0.9980695458637723, cdf(2.5,2,1)),
        ?assertEqual(0.9998765901959134, cdf(3.0,2,1)).

invcdf_test() ->
        ?assertEqual(0.0, invcdf(0.0,2,1)),
        ?assertEqual(0.4723807270774388, invcdf(0.2,2,1)),
        ?assertEqual(0.8325546111576977, invcdf(0.5,2,1)),
        ?assertEqual(1.0972569454443821, invcdf(0.7,2,1)),
        ?assertEqual(1.5174271293851465, invcdf(0.9,2,1)).

invcdf_error_test() ->
        ?assertEqual({error,"Invalid probability"}, invcdf(-0.1,0,10)),
        ?assertEqual({error,"Invalid probability"}, invcdf(1.1,0,10)).
-endif.


