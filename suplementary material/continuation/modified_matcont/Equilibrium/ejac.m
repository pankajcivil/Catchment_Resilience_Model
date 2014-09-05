function j=ejac(x,p)
global eds cds

if cds.options.SymDerivative >= 1
    % This feature for efficicent calculation of the Jacobian was implemented for 
    % Peterson et al. 2012 PDE model for the hill-slope Bousinesq equation.
    if isempty(cds.options.odeset.Jacobian_Options);        
        error(' Calculation of the model Jacobian must be vectorised in order to achieve acceptable performance!');
    else
        j = full( feval(eds.Jacobian,  [] , x, p{:}, cds.options.odeset.Jacobian_Options.g, cds.options.odeset.Jacobian_Options.pattern ));
    end    
else
  for i=lds.phases
    x1 = x; x1(i) = x1(i)-cds.options.Increment;
    x2 = x; x2(i) = x2(i)+cds.options.Increment;
    j(:,i) = feval(eds.func, 0, x2, p{:})-feval(eds.func, 0, x1, p{:});
  end
  j = j/(2*cds.options.Increment);
end

