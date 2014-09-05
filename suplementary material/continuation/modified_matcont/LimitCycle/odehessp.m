function h=odehessp(t,x,p)%dxdp
global lds cds

if cds.options.SymDerivativeP >= 2
  h = feval(lds.HessiansP, t, x, p{:});
  h = h(:,:,lds.ActiveParams);
else
  for i=lds.ActiveParams
    p1 = p; p1{i} = p1{i}-cds.options.Increment;
    p2 = p; p2{i} = p2{i}+cds.options.Increment;
    h(:,:,i) = odejac(t, x,p2)-odejac(t, x,p1);
  end
  h = h(:,:,lds.ActiveParams)/(2*cds.options.Increment);
end
