% ==============================================================================
%  Student-t distribution
% ==============================================================================
-module(tdist).

-export([tdistpdf/4, tdistcdf/4, tdistinv/4]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(SQR(X), ((X)*(X))).

% ------------------------------------------------------------------------------
%  tdistpdf - Probability density function
% ------------------------------------------------------------------------------
tdistpdf(_,_,Sig,_) when Sig =< 0 -> 
        {error,"Sig parameter =< 0 in tdistpdf."};
tdistpdf(_,_,_,Nu) when Nu =< 0 -> 
        {error,"Nu parameter =< 0 in tdistpdf."};

tdistpdf(T,Mu,Sig,Nu) ->
        Np = 0.5*(Nu + 1.0),
        Fac = incgammabeta:gammaln(Np)-incgammabeta:gammaln(0.5*Nu),
        math:exp(-Np*math:log(1.0+?SQR((T-Mu)/Sig)/Nu)+Fac)
         / (math:sqrt(3.14159265358979324*Nu)*Sig).

% ------------------------------------------------------------------------------
%  tdistcdf - Cumulative distribution function
% ------------------------------------------------------------------------------
tdistcdf(_,_,Sig,_) when Sig =< 0 -> 
        {error,"Sig parameter =< 0 in tdistcdf."};
tdistcdf(_,_,_,Nu) when Nu =< 0 -> 
        {error,"Nu parameter =< 0 in tdistcdf."};

tdistcdf(T,Mu,Sig,Nu) ->
        P = 0.5*incgammabeta:betai(0.5*Nu, 0.5, Nu/(Nu+?SQR((T-Mu)/Sig))),
        if T >= Mu -> 1.0 - P
         ; T <  Mu -> P
        end.

        
% ------------------------------------------------------------------------------
%  tdistinv - Inverse cumulative distribution function
% ------------------------------------------------------------------------------
tdistinv(_,_,Sig,_) when Sig =< 0 -> 
        {error,"Sig parameter =< 0 in tdistinv."};
tdistinv(_,_,_,Nu) when Nu =< 0 -> 
        {error,"Nu parameter =< 0 in tdistinv."};

tdistinv(P,_,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};

tdistinv(P,Mu,Sig,Nu) -> 
        X = incgammabeta:invbetai(2.0*min(P,1.0-P), 0.5*Nu, 0.5), 
        if P >= 0.5 -> Mu + Sig*math:sqrt(Nu*(1.0-X)/X)
         ; P <  0.5 -> Mu - Sig*math:sqrt(Nu*(1.0-X)/X)
        end.




% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

tdistpdf_test() ->
        ?assertEqual(0.37960668982249457, tdistpdf(0,0,1,5)),
        ?assertEqual(0.3891083839660307, tdistpdf(0,0,1,10)),
        ?assertEqual(0.23036198922913836, tdistpdf(-1,0,1,10)),
        ?assertEqual(0.23036198922913836, tdistpdf(1,0,1,10)).

tdistpdf_error_test() ->
        ?assertEqual({error,"Sig parameter =< 0 in tdistpdf."},  tdistpdf(0.0,0,-1,1)),
        ?assertEqual({error,"Nu parameter =< 0 in tdistpdf."},  tdistpdf(0.0,0,1,-1)).

tdistcdf_test() ->
        ?assertEqual(0.5, tdistcdf(0,0,1,5)),
        ?assertEqual(0.5, tdistcdf(0,0,1,10)),
        ?assertEqual(0.8295534338489696, tdistcdf(1,0,1,10)),
        ?assertEqual(0.17044656615103038, tdistcdf(-1,0,1,10)).

tdistcdf_param_test() ->
        ?assertEqual({error,"Sig parameter =< 0 in tdistcdf."},  tdistcdf(0.0,0,-1,1)),
        ?assertEqual({error,"Nu parameter =< 0 in tdistcdf."},  tdistcdf(0.0,0,1,-1)).

tdistinv_test() ->
        ?assertEqual(-1.4758840488244813, tdistinv(0.1,0,1,5)),
        ?assertEqual(-1.372183641110338, tdistinv(0.1,0,1,10)),
        ?assertEqual(0.0, tdistinv(0.5,0,1,10)),
        ?assertEqual(0.5415280387550168, tdistinv(0.7,0,1,10)),
        ?assertEqual(2.7637694581126953, tdistinv(0.99,0,1,10)).

tdistinv_param_test() ->
        ?assertEqual({error,"Sig parameter =< 0 in tdistinv."},  tdistinv(0.0,0,-1,1)),
        ?assertEqual({error,"Nu parameter =< 0 in tdistinv."},  tdistinv(0.0,0,1,-1)),
        ?assertEqual({error,"Invalid probability"},  tdistinv(-0.1,0,1,1)).

-endif.




