% ==============================================================================
%  Gamma distribution
% ==============================================================================
-module(gammadist).

-export([gampdf/3, gamcdf/3, gaminv/3, gammaln/1]).

-ifdef(TEST).
-include_lib("eunit/include/eunit.hrl").
-endif.

-define(FPMIN, 1.00208418000449e-292).
-define(EPS, 2.22044604925031e-16).

% ------------------------------------------------------------------------------
%  gampdf - Uniform probability density function
% ------------------------------------------------------------------------------
gampdf(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
gampdf(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

gampdf(X,_,_) when X =< 0 -> 0.0;
gampdf(X,Alpha,Beta) ->
        Fac = Alpha*math:log(Beta)-gammaln(Alpha),
        math:exp(-Beta*X+(Alpha-1)*math:log(X)+Fac).

% ------------------------------------------------------------------------------
%  gamcdf - Uniform cumulative distribution function
% ------------------------------------------------------------------------------
gamcdf(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
gamcdf(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

gamcdf(X,_,_) when X =< 0 -> 0.0;
gamcdf(X,Alpha,Beta) ->
        gammap(Alpha,Beta*X).

        
% ------------------------------------------------------------------------------
%  gaminv - Inverse gamma distribution function
% ------------------------------------------------------------------------------
gaminv(_,Alpha,_) when Alpha =< 0 -> {error,"Alpha parameter =< 0."};
gaminv(_,_,Beta) when Beta =< 0 -> {error,"Beta parameter =< 0."};

gaminv(P,_,_) when P < 0 orelse P > 1 -> {error,"Invalid probability"};
gaminv(P,Alpha,Beta) ->
        invgammap(P,Alpha)/Beta.


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

invgammap(P,A) when P >= 1 -> math:max(100.0, A + 100.0 * math:sqrt(A));
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

% ==============================================================================
%  EUnit tests
% ------------------------------------------------------------------------------
-ifdef(TEST).

gammaln_test() ->
        ?assertEqual(2.036327503417712, gammaln(0.123)),
        ?assertEqual(-0.09447840768115956, gammaln(1.234)),
        ?assertEqual(18.734347511936445, gammaln(12.5)).

gammap_test() ->
        ?assertEqual({error, "Bad args in gammap."}, gammap(-1.0, 1.0)),
        ?assertEqual({error, "Bad args in gammap."}, gammap( 1.0,-1.0)),
        ?assertEqual(0.0, gammap(1.0, 0.0)).

gampdf_test() ->
        ?assertEqual(0.0, gampdf(0.0,2,2)),
        ?assertEqual(0.2651846416468157, gampdf(2.3,3,1)),
        ?assertEqual(1.9199765904692206e-04, gampdf(20,20,2)),
        ?assertEqual(0.5413411329464507, gampdf(1.0,2,2)).

gampdf_error_test() ->
        ?assertEqual({error,"Alpha parameter =< 0."}, gampdf(0.0,-1, 1)),
        ?assertEqual({error,"Beta parameter =< 0."},  gampdf(0.0, 1,-1)).

gamcdf_test() ->
        ?assertEqual(0.0, gamcdf(0.0,1,1)),
        ?assertEqual(0.5132987982791589, gamcdf(100,100,1)),
        ?assertEqual(0.5297427331607533, gamcdf(20,20,1)),
        ?assertEqual(0.9998236971022615, gamcdf(20,20,2)),
        ?assertEqual(0.12478121503252411, gamcdf(15,20,1)),
        ?assertEqual(0.9595723180054873, gamcdf(5,2,1)).

gamcdf_error_test() ->
        ?assertEqual({error,"Alpha parameter =< 0."}, gamcdf(0.0,-1, 1)),
        ?assertEqual({error,"Beta parameter =< 0."},  gamcdf(0.0, 1,-1)).

gaminv_test() ->
        ?assertEqual(0.0, gaminv(0.0,1,1)),
        ?assertEqual(0.6931471805599453, gaminv(0.5,1,1)),
        ?assertEqual(0.22746821155978653, gaminv(0.5,0.5,1)),
        ?assertEqual(4.670908882795983, gaminv(0.5,5,1)),
        ?assertEqual(2.432591025962664, gaminv(0.1,5,1)),
        ?assertEqual(4.865182051925328, gaminv(0.1,5,1/2)),
        ?assertEqual(1.216295512981332, gaminv(0.1,5,2)),
        ?assertEqual(3.9967947930263152, gaminv(0.9,5,2)).

-endif.


