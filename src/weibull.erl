% ==============================================================================
%  Weibull distribution
% ==============================================================================
-module(weibull).

-export([weibullpdf/3, weibullcdf/3, weibullinv/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.


% ------------------------------------------------------------------------------
%  weibullpdf - Uniform probability density function
% ------------------------------------------------------------------------------
weibullpdf(X,_,_) when X < 0 -> 0.0;
weibullpdf(X,A,B) ->
        (A/B)*math:pow(X/B,A-1)*math:exp(-math:pow(X/B,A)).

% ------------------------------------------------------------------------------
%  weibullpcdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
weibullcdf(X,_,_) when X < 0 -> 0.0;
weibullcdf(X,A,B) ->
        1-math:exp(-math:pow(X/B,A)).

        
% ------------------------------------------------------------------------------
%  weibullpinv - Inverse weibull distribution function
% ------------------------------------------------------------------------------
weibullinv(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
weibullinv(P,A,B) ->
        B*math:pow(-math:log(1-P),1/A).



% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

weibullpdf_test() ->
        ?assertEqual(0.0, weibullpdf(-1.0,2,1)),
        ?assertEqual(0.0, weibullpdf(0.0,2,1)),
        ?assertEqual(0.7788007830714049, weibullpdf(0.5,2,1)),
        ?assertEqual(0.7357588823428847, weibullpdf(1.0,2,1)),
        ?assertEqual(0.3161976736855930, weibullpdf(1.5,2,1)),
        ?assertEqual(0.07326255555493671, weibullpdf(2.0,2,1)),
        ?assertEqual(0.009652270681138546, weibullpdf(2.5,2,1)),
        ?assertEqual(0.0007404588245200774, weibullpdf(3.0,2,1)).


weibullcdf_test() ->
        ?assertEqual(0.0, weibullcdf(0.0,2,1)),
        ?assertEqual(0.22119921692859512, weibullcdf(0.5,2,1)),
        ?assertEqual(0.6321205588285577, weibullcdf(1.0,2,1)),
        ?assertEqual(0.8946007754381357, weibullcdf(1.5,2,1)),
        ?assertEqual(0.9816843611112658, weibullcdf(2.0,2,1)),
        ?assertEqual(0.9980695458637723, weibullcdf(2.5,2,1)),
        ?assertEqual(0.9998765901959134, weibullcdf(3.0,2,1)).

weibullinv_test() ->
        ?assertEqual(0.0, weibullinv(0.0,2,1)),
        ?assertEqual(0.4723807270774388, weibullinv(0.2,2,1)),
        ?assertEqual(0.8325546111576977, weibullinv(0.5,2,1)),
        ?assertEqual(1.0972569454443821, weibullinv(0.7,2,1)),
        ?assertEqual(1.5174271293851465, weibullinv(0.9,2,1)).

-endif.


