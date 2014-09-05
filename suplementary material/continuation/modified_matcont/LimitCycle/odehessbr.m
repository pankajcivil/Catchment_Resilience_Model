function h=odehessbr(t,x,p)
global lds cds

if cds.options.SymDerivativeP >= 2
  h = feval(lds.HessiansP, t, x, p{:});
  h = h(:,:,lds.BranchParam);
else
  for i=lds.BranchParam
    p1 = p; p1{i} = p1{i}-cds.options.Increment;
    p2 = p; p2{i} = p2{i}+cds.options.Increment;
    h(:,:,i) = odejac(t,x,p2)-odejac(t,x,p1);
  end
  h = h(:,:,lds.BranchParam)/(2*cds.options.Increment);
end
