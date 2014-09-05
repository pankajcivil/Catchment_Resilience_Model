function h=odehess(t,x,p)
global lds cds

if cds.options.SymDerivative >= 2
  h = feval(lds.Hessians, t, x, p{:});
else
  for i=lds.phases
    x1 = x; x1(i) = x1(i)-cds.options.Increment;
    x2 = x; x2(i) = x2(i)+cds.options.Increment;
    h(:,:,i) = odejac(t,x2,p)-odejac(t,x1,p);
  end
  h = h/(2*cds.options.Increment);
end
