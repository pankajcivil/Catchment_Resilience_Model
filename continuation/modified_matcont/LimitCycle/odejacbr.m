function j=odejacbr(t,x,p)
global lds cds

if cds.options.SymDerivativeP >= 1
  j = feval(lds.JacobianP, t, x, p{:});
  j = j(:,lds.BranchParam);
else
  for i=lds.BranchParam
    p1 = p; p1{i} = p1{i}-cds.options.Increment;
    p2 = p; p2{i} = p2{i}+cds.options.Increment;
    j(:,i) = feval(lds.func, t, x, p2{:})-feval(lds.func, t, x, p1{:});
  end
    j = j(:,lds.BranchParam)/(2*cds.options.Increment);
end
