% -------------------------------------------------------------
% functions defining the BVP for limitcycles
% -------------------------------------------------------------

function f = BVP_LC_f(odefile,t,xp,p,T, time)
f = t-T*feval(odefile, time, xp, p{:});
