function h=odehesspbr(t,x,p)%dbrdp
global lds cds

for i=lds.ActiveParams
    p1 = p; p1{i} = p1{i}-cds.options.Increment;
    p2 = p; p2{i} = p2{i}+cds.options.Increment;
    h(:,:,i) = odejacbr(t,x,p2)-odejacbr(t,x,p1);
end
h = h(:,:,lds.ActiveParams)/(2*cds.options.Increment);

