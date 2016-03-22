% ==============================================================================
%  Incomplete gamma and beta functions with friends
% ==============================================================================
-module(incgammabeta).

-export([gammaln/1, gammap/2, invgammap/2, betai/3, invbetai/3]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(FPMIN, 1.00208418000449e-292).
-define(EPS, 2.22044604925031e-16).
-define(SWITCH, 3000).
%
% ------------------------------------------------------------------------------
%  Gamma functions
% ------------------------------------------------------------------------------

%% Natural logarithm of Gamma function
gammaln(X) when X =< 0.0 -> {error, "Bad arg to gammln."};
gammaln(X) -> 
        Coef = [ 57.1562356658629235,    -59.5979603554754912, 
                 14.1360979747417471,     -0.491913816097620199,
                  0.339946499848118887e-4, 0.465236289270485756e-4,
                 -0.983744753048795646e-4, 0.158088703224912494e-3,
                 -0.210264441724104883e-3, 0.217439618115212643e-3,
                 -0.164318106536763890e-3, 0.844182239838527433e-4,
                 -0.261908384015814087e-4, 0.368991826595316234e-5 ],

        Tmp1 = X+5.24218750000000000,
        Tmp2 = (X+0.5)*math:log(Tmp1)-Tmp1,
        Ys = [X+I || I <- lists:seq(1,14)],
        Ser = lists:sum([0.999999999999997092|[C/Y || {C,Y} <- lists:zip(Coef,Ys)]]),
        Tmp2+math:log(2.5066282746310005*Ser/X).

%% Incomplete Gamma function
gammap(A,X) when X < 0 orelse A =< 0 -> {error, "Bad args in gammap."};
gammap(_,0.0) -> 0.0;
gammap(A,X) when A >= 100 -> gammapapprox(A,X,1);
gammap(A,X) when X =< A + 1.0 -> gser(A,X);
gammap(A,X) -> 1.0 - gcf(A,X).


%% Abscissas for Gauss-Legendre quadrature
gauleg18_y() ->
        [ 0.0021695375159141994, 0.011413521097787704, 0.027972308950302116,
          0.051727015600492421,  0.082502225484340941, 0.12007019910960293,
          0.16415283300752470,   0.21442376986779355,  0.27051082840644336,
          0.33199876341447887,   0.39843234186401943,  0.46931971407375483,
          0.54413605556657973,   0.62232745288031077,  0.70331500465597174,
          0.78649910768313447,   0.87126389619061517,  0.95698180152629142 ].

%% Abscissas for Gauss-Legendre quadrature
gauleg18_w() -> 
        [ 0.0055657196642445571, 0.012915947284065419, 0.020181515297735382,
          0.027298621498568734,  0.034213810770299537, 0.040875750923643261,
          0.047235083490265582,  0.053244713977759692, 0.058860144245324798,
          0.064039797355015485,  0.068745323835736408, 0.072941885005653087,
          0.076598410645870640,  0.079687828912071670, 0.082187266704339706,
          0.084078218979661945,  0.085346685739338721, 0.085983275670394821 ].

%% Approximate incomplete gamma function by Gauss-Legendre quadrature
gammapapprox(A,X,Psig) ->
        A1 = A-1.0,
        LnA1 = math:log(A1),
        SqrtA1 = math:sqrt(A1),
        Gln = gammaln(A),
        Xu = if X >  A1 -> max(A1 + 11.5*SqrtA1, X + 6.0*SqrtA1)
              ; X =< A1 -> max(0.0, min(A1 - 7.5*SqrtA1, X - 5.0*SqrtA1))
             end,
        Ts = [X + (Xu-X)*Y || Y <- gauleg18_y()],
        Sum = lists:sum([W*math:exp(-(T-A1) + A1*(math:log(T)-LnA1)) || {W,T} <- lists:zip(gauleg18_w(),Ts)]),
        Ans = Sum*(Xu-X)*math:exp(A1*(LnA1-1.0) - Gln),
        case Psig of
                1 -> if Ans >  0.0 -> 1.0 - Ans
                      ; Ans =< 0.0 -> -Ans
                     end;
                0 -> if Ans >= 0.0 -> Ans
                      ; Ans <  0.0 -> 1.0 + Ans
                     end
        end.


%% Approximate incomplete gamma function by seriies expansion
gser(A,X) -> gser_(A, X, A, 1.0/A, 1.0/A).

%% Tail recusrion for gser
gser_(A,X,_,Del,Sum) when abs(Del) < abs(Sum)*?EPS -> 
        Sum*math:exp(-X+A*math:log(X)-gammaln(A));
gser_(A,X,Ap,Del,Sum) ->
        Ap2 = Ap + 1,
        Del2 = Del*X/Ap2,
        Sum2 = Sum+Del2,
        gser_(A, X, Ap2, Del2, Sum2).


%%
gcf(A,X) -> 
        B = X + 1.0 - A,
        D = 1.0/B,
        gcf_(A,X,B,1.0/?FPMIN,D,D,1).

gcf_(A,X,_,C,D,H,_) when abs(D*C-1) =< ?EPS ->
        math:exp(-X+A*math:log(X)-gammaln(A))*H;

gcf_(A,X,B,C,D,H,I) ->
        An = -I*(I-A),
        B2 = B + 2,
        D2 = if abs(An*D+B2) < ?FPMIN -> 1.0/?FPMIN
              ; abs(An*D+B2) >=?FPMIN -> 1.0/(An*D+B2)
             end,
        C2 = if abs(B2+An/C) < ?FPMIN -> ?FPMIN
              ; abs(B2+An/C) >=?FPMIN -> B2+An/C
             end,
        gcf_(A,X,B2,C2,D2,H*D2*C2,I+1).


%% Inverse of gammap
invgammap(_,A) when A =< 0 -> {error, "Parameter a must be larger then zero in invgammap."};

invgammap(P,A) when P >= 1 -> max(100.0, A + 100.0 * math:sqrt(A));
invgammap(P,_) when P =< 0 -> 0.0;

invgammap(P,A) when A > 1 ->
        PP = min(P,1-P),
        T = math:sqrt(-2*math:log(PP)),
        Sign = if P <  0.5 -> -1
                ; P >= 0.5 -> 1
               end,
        X0 = Sign*(2.30753+T*0.27061)/(1+T*(0.99229+T*0.04481)) - T,
        X = max(1.0e-3,A*math:pow(1-1/(9*A)-X0/(3*math:sqrt(A)),3)),
        invgammap_(P,A,X,12);

invgammap(P,A) ->
        T = 1 - A*(0.253+A*0.12),
        X = if P <  T -> math:pow(P/T,1.0/A)
             ; P >= T -> 1 - math:log(1-(P-T)/(1-T))
            end,
        invgammap_(P,A,X,12).


invgammap_(_,_,X,_) when X =< 0 -> 0.0; 
invgammap_(_,_,X,0) -> X; 
invgammap_(P,A,X,J) -> 
        A1 = A - 1,
        Gln = gammaln(A),
        Err = gammap(A,X) - P,
        T0 = if A >  1 -> 
                        Lna1 = math:log(A1),
                        Afac = math:exp(A1*(Lna1-1)-Gln),
                        Afac*math:exp(-(X-A1)+A1*(math:log(X)-Lna1))
             ; A =< 1 -> 
                        math:exp(-X+A1*math:log(X)-Gln)
            end,
        U = Err/T0,
        T = U/(1-0.5*min(1.0, U*((A-1)/X - 1))),
        Xi = if X-T =< 0 -> 0.5*X
              ; X-T >  0 -> X-T
             end,
        Stop = abs(T) < ?EPS*X,
        if Stop -> invgammap_(P,A,Xi,0)
         ; true -> invgammap_(P,A,Xi,J-1)
        end.


% ------------------------------------------------------------------------------
%  Beta functions
% ------------------------------------------------------------------------------

%% Parameter checks 
betai(A,_,_) when A =< 0.0 -> {error, "Bad arg A in function betai."};
betai(_,B,_) when B =< 0.0 -> {error, "Bad arg B in function betai."};
betai(_,_,X) when X < 0.0 orelse X > 1.0 -> {error, "Bad arg X in function betai."};
%% Border cases
betai(_,_,0) -> 0.0;
betai(_,_,1) -> 1.0;
% betai(A,B,1) -> math:exp(gammaln(A) + gammaln(B) - gammaln(A+B));
%% Use approximations for extreme params
betai(A,B,X) when A > ?SWITCH andalso B > ?SWITCH -> betaiapprox(A,B,X);
%% The rest use the normal definition
betai(A,B,X) -> 
        Bt = math:exp(gammaln(A+B) - gammaln(A) - gammaln(B) + A*math:log(X)
                      + B*math:log(1.0-X)),
        if X <  ((A+1.0)/(A+B+2.0)) -> Bt*betacf(A,B,X)/A 
         ; X >= ((A+1.0)/(A+B+2.0)) -> 1.0 - Bt*betacf(B,A,1.0-X)/B
        end.


%% Use a minimum value if absollute value is less than the cap.
absmin(X,Cap) ->
	if abs(X) <  Cap -> Cap
         ; abs(X) >= Cap -> X
        end.


%% 
betacf(A,B,X) ->
	Qab = A + B,
	Qap = A + 1.0,
	C   = 1.0,
	D  = 1.0/absmin(1.0 - Qab*X/Qap, ?FPMIN),
	H = D,
	betacf_(A,B,X,C,D,H,1).


%% betacf helper function
betacf_(_,_,_,_,_,H,10000) -> H;
betacf_(A,B,X,C,D,H,M) ->
	Qab = A + B,
	Qap = A + 1.0,
	Qam = A - 1.0,
	M2 = 2*M,

	AA = M*(B-M)*X/((Qam+M2)*(A+M2)),
	D0 = 1.0/absmin(1.0+AA*D, ?FPMIN),
	C0 = absmin(1.0+AA/C, ?FPMIN),
	H0 = H*D0*C0,

	AB = -(A+M)*(Qab+M)*X/((A+M2)*(Qap+M2)),
	D1 = 1.0/absmin(1.0+AB*D0, ?FPMIN),
	C1 = absmin(1.0+AB/C0, ?FPMIN),
	H1 = H0*D1*C1,

	N  = case abs(D1*C1 - 1.0) =< ?EPS of
               true  -> 10000;
	       false -> M + 1
	     end,
	betacf_(A,B,X,C1,D1,H1,N).

%%
betaiapprox(_,_,_) -> {error, "Not implemented"}.


%%
invbetai(P,_,_) when P =< 0 -> 0.0;
invbetai(P,_,_) when P >= 1 -> 1.0;
invbetai(P,A,B) when A >= 1 andalso B >= 1 -> 
        PP = min(P, 1-P),
        T = math:sqrt(-2.0*math:log(PP)),
        X0 = (2.30753+T*0.27061)/(1.0+T*(0.99229+T*0.04481)) - T,
        X = if P =< 0.5 -> -X0
             ; P >  0.5 ->  X0
            end,
        Al = (X*X-3.0)/6.0,
        H = 2.0/(1.0/(2.0*A-1.0)+1.0/(2.0*B-1.0)),
        W = (X*math:sqrt(Al+H)/H) 
                - (1.0/(2.0*B-1)-1.0/(2.0*A-1.0))
                * (Al+5.0/6.0-2.0/(3.0*H)),
        Xn = A/(A+B*math:exp(2.0*W)),
        invbetai_(P,A,B,Xn,0);

invbetai(P,A,B) -> 
	Lna = math:log(A/(A+B)),
	Lnb = math:log(B/(A+B)),
	T = math:exp(A*Lna)/A,
	U = math:exp(B*Lnb)/B,
	W = T + U,
	X = if P < T/W -> math:pow(A*W*P,1.0/A)
	     ; true    -> 1.0 - math:pow(B*W*(1.0-P),1.0/B)
	    end,
        invbetai_(P,A,B,X,0).


%% 
invbetai_(_,_,_,0.0,_) -> 0.0;
invbetai_(_,_,_,1.0,_) -> 1.0;
invbetai_(_,_,_,X,10) -> X;

invbetai_(P,A,B,X,N) -> 
        Afac = -gammaln(A)-gammaln(B)+gammaln(A+B),
        Err = betai(A,B,X) - P,
        T = math:exp((A-1)*math:log(X)+(B-1)*math:log(1.0-X) + Afac),
        U = Err/T,
        Tn = U/(1.0 - 0.5*min(1.0,U*((A-1)/X - (B-1)/(1.0-X)))),
        X0 = X - Tn,
        Xn = if X0 =< 0.0 -> 0.5*(X0 + T)
              ; X0 >= 1.0 -> 0.5*(X0 + T + 1.0)
              ; true      -> X0
             end,

        invbetai_(P,A,B,Xn,N+1).


% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

gammaln_test() ->
        ?assertEqual(2.036327503417712, gammaln(0.123)),
        ?assertEqual(-0.09447840768115956, gammaln(1.234)),
        ?assertEqual(18.734347511936445, gammaln(12.5)).

gammaln_param_test() ->
        ?assertEqual({error,"Bad arg to gammln."}, gammaln(-0.123)).

gammap_param_test() ->
        ?assertEqual({error, "Bad args in gammap."}, gammap(-1.0, 1.0)),
        ?assertEqual({error, "Bad args in gammap."}, gammap( 1.0,-1.0)).

gammap_test() ->
        ?assertEqual(0.0, gammap(1.0, 0.0)).

gammapapprox_test() ->
        ?assertEqual(0.147137348442258092, gammap(111, 100)).

invgammap_param_test() ->
        ?assertEqual({error, "Parameter a must be larger then zero in invgammap."},
                    invgammap(1.0,-1)).

invgammap_large_p_test() ->
        ?assertEqual(100.0, invgammap(2.0,0.1)),
        ?assertEqual(100.0, invgammap(200.0,0.1)),
        ?assertEqual(143.421356237309510107, invgammap(2.0,2.0)).

invgammap_small_A_test() ->
        ?assertEqual(1.119032936887381702, invgammap(0.98,0.1)).

betai_param_test() ->
        ?assertEqual({error, "Bad arg A in function betai."}, betai(-1,1,1)),
        ?assertEqual({error, "Bad arg B in function betai."}, betai(1,-1,1)),
        ?assertEqual({error, "Bad arg X in function betai."}, betai(1,1,-1)).

betai_test() ->
        ?assertEqual(1.0, betai(1,1,1)),
        ?assertEqual(0.0, betai(1,1,0)),
        ?assertEqual(9/10, betai(1,1,0.9)),
        ?assertEqual(0.7499999999999998, betai(1,2,0.5)),
        ?assertEqual(0.99, betai(1,2,0.9)),
        ?assertEqual(0.11086997774500017, betai(3,10,0.1)),
        ?assertEqual(5.455000000000031e-9, betai(10,3,0.1)),
        ?assertEqual(1.0, betai(2,1,1)).


invbetai_test() ->
        ?assertEqual(0.0, invbetai(-1,1,1)),
        ?assertEqual(0.0, invbetai(0,1,1)),
        ?assertEqual(1.0, invbetai(1,1,1)),
        ?assertEqual(1.0, invbetai(10,1,1)),
        ?assertEqual(0.4999999999999997, invbetai(0.5,2,2)),
        ?assertEqual(0.950548152829731596, invbetai(0.9,10,2)),
        ?assertEqual(0.06424517809580489, invbetai(0.001,7,20)).

-endif.
