% -------------------------------------------------------------
% functions defining the jacobian of the BVP for BPC bifurcations
% -------------------------------------------------------------

% function
function f = BVP_BPC_jac_f(odefile,xp,p,T,tp)
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
f = [wploc-T*sysjac    lds.BP_psi((tp-1)*lds.ncol_coord + lds.col_coords)'];
