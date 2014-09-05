% -------------------------------------------------------------
% functions defining the jacobian of the BVP for limitcycles
% -------------------------------------------------------------
function f = BVP_LC_jac_f(odefile,xp,p,T,tp)
global lds
sysjacp = zeros(lds.ncol_coord,length(lds.ActiveParams));
wploc = lds.wp/lds.dt(tp);
range = lds.phases;
% xp:value of polynomial on each collocation point
for c=lds.cols
    
  % TJP get mesh time and add to all ode calls below
  meshindex= lds.ncol*(tp-1)+c;    
  mesh_time = lds.finemsh(meshindex);
  
  xt = xp(:,c);
  frhs(range) = feval(odefile, mesh_time, xt, p{:});
  jac=odejac(mesh_time, xt,p);
  sysjac(range,:) = fastkron(lds.ncol,lds.nphase,lds.wt(:,c)',jac);
  jacp = odejacp(mesh_time, xt,p);
  sysjacp(range,:) = jacp(:,lds.ActiveParams);
  range = range + lds.nphase;
end
f = [wploc-T*sysjac    -frhs'    -T*sysjacp];
