% ==============================================================================
%  Student-t distribution
% ==============================================================================
-module(tdist).

-export([pdf/4, cdf/4, invcdf/4]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(SQR(X), ((X)*(X))).

% ------------------------------------------------------------------------------
%  pdf - Probability density function
% ------------------------------------------------------------------------------
pdf(_,_,Sig,_) when Sig =< 0 -> 
        {error,"Sig parameter =< 0 in pdf."};
pdf(_,_,_,Nu) when Nu =< 0 -> 
        {error,"Nu parameter =< 0 in pdf."};

pdf(T,Mu,Sig,Nu) ->
        Np = 0.5*(Nu + 1.0),
        Fac = incgammabeta:gammaln(Np)-incgammabeta:gammaln(0.5*Nu),
        math:exp(-Np*math:log(1.0+?SQR((T-Mu)/Sig)/Nu)+Fac)
         / (math:sqrt(3.14159265358979324*Nu)*Sig).

% ------------------------------------------------------------------------------
%  cdf - Cumulative distribution function
% ------------------------------------------------------------------------------
cdf(_,_,Sig,_) when Sig =< 0 -> 
        {error,"Sig parameter =< 0 in cdf."};
cdf(_,_,_,Nu) when Nu =< 0 -> 
        {error,"Nu parameter =< 0 in cdf."};

cdf(T,Mu,Sig,Nu) ->
        P = 0.5*incgammabeta:betai(0.5*Nu, 0.5, Nu/(Nu+?SQR((T-Mu)/Sig))),
        if T >= Mu -> 1.0 - P
         ; T <  Mu -> P
        end.

        
% ------------------------------------------------------------------------------
%  invcdf - Inverse cumulative distribution function
% ------------------------------------------------------------------------------
invcdf(_,_,Sig,_) when Sig =< 0 -> 
        {error,"Sig parameter =< 0 in invcdf."};
invcdf(_,_,_,Nu) when Nu =< 0 -> 
        {error,"Nu parameter =< 0 in invcdf."};

invcdf(P,_,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};

invcdf(P,Mu,Sig,Nu) -> 
        X = incgammabeta:invbetai(2.0*min(P,1.0-P), 0.5*Nu, 0.5), 
        if P >= 0.5 -> Mu + Sig*math:sqrt(Nu*(1.0-X)/X)
         ; P <  0.5 -> Mu - Sig*math:sqrt(Nu*(1.0-X)/X)
        end.




% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

pdf_test() ->
        ?assertEqual(0.37960668982249457, pdf(0,0,1,5)),
        ?assertEqual(0.3891083839660307, pdf(0,0,1,10)),
        ?assertEqual(0.23036198922913836, pdf(-1,0,1,10)),
        ?assertEqual(0.23036198922913836, pdf(1,0,1,10)).

pdf_error_test() ->
        ?assertEqual({error,"Sig parameter =< 0 in pdf."}, pdf(0.0,0,-1,1)),
        ?assertEqual({error,"Nu parameter =< 0 in pdf."}, pdf(0.0,0,1,-1)).

cdf_test() ->
        ?assertEqual(0.5, cdf(0,0,1,5)),
        ?assertEqual(0.5, cdf(0,0,1,10)),
        ?assertEqual(0.8295534338489696, cdf(1,0,1,10)),
        ?assertEqual(0.17044656615103038, cdf(-1,0,1,10)).

cdf_param_test() ->
        ?assertEqual({error,"Sig parameter =< 0 in cdf."}, cdf(0.0,0,-1,1)),
        ?assertEqual({error,"Nu parameter =< 0 in cdf."}, cdf(0.0,0,1,-1)).

invcdf_test() ->
        ?assertEqual(-1.4758840488244813, invcdf(0.1,0,1,5)),
        ?assertEqual(-1.372183641110338, invcdf(0.1,0,1,10)),
        ?assertEqual(0.0, invcdf(0.5,0,1,10)),
        ?assertEqual(0.5415280387550168, invcdf(0.7,0,1,10)),
        ?assertEqual(2.7637694581126953, invcdf(0.99,0,1,10)).

invcdf_param_test() ->
        ?assertEqual({error,"Sig parameter =< 0 in invcdf."}, invcdf(0.0,0,-1,1)),
        ?assertEqual({error,"Nu parameter =< 0 in invcdf."}, invcdf(0.0,0,1,-1)),
        ?assertEqual({error,"Invalid probability"}, invcdf(-0.1,0,1,1)).

-endif.




