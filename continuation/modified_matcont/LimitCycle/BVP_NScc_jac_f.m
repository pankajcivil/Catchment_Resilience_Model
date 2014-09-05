% -------------------------------------------------------------
% functions defining the jacobian of the BVP for flip (PD) bifurcations
% -------------------------------------------------------------

% function
function f = BVP_PD_jac_f(odefile,xp,p,T,theta,tp)
global lds
wploc = lds.wp/lds.dt(tp);
range = lds.phases;
for c=lds.cols
    
  % TJP get mesh time and add to all ode calls below
  meshindex= lds.ncol*(tp-1)+c;    
  mesh_time = lds.finemsh(meshindex);
  
  sysjac(range,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',odejac(mesh_time,xp(:,c),p));
  range = range + lds.nphase;
end
b=lds.wt([1 1 1 2 2 2 3 3 3 4 4 4 5 5 5],[1 1 1 2 2 2 3 3 3 4 4 4 ])';

f = [wploc-T*sysjac+2*1i*theta*b];
