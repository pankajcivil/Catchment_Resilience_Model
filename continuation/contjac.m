 function jac = contjac(x,varargin)
% function jac = contjac(x)
% cds.pJac

%
%  jacobian(x)
%
%    Calculates jacobian matrix of F(x), which is to be found in curve file, which is global
global cds

%TJP: check if the jacopbian needs to be stored
storeJac=true;
if ~isempty(varargin)
    if islogical(varargin{1});
        storeJac=varargin{1};
    end
end

if nargin <= 1 
  error('contjac needs a point');
end
if ~isempty(cds.pJacX)
  if x == cds.pJacX
    jac = cds.pJac;    
    return;
  end
end
try
    symjac  = cds.symjac;
catch symjac=0;end

if ~storeJac;
    cds.pJac = [];
    cds.pJacX =[];
end

if symjac
%if cds.options.SymDerivative >=1
  jac =  feval(cds.curve_jacobian, x);  
else
  x1 = x;
  x2 = x;
  % WM: mmm... is this really needed?
  jac = ones(cds.ndim-1,cds.ndim)*NaN;

  for j=1:cds.ndim  %cols
    x1(j) = x(j) - cds.options.Increment;
    x2(j) = x(j) + cds.options.Increment;

    % WM: Removed temporary result variables
    jac(:,j) = feval(cds.curve_func, x2)-feval(cds.curve_func, x1);
    
    x1(j) = x(j);
    x2(j) = x(j);
  end

  % WM: Moved out of loop
  % TJP: removed setting of jac
  cds.pJac = jac/(2*cds.options.Increment);
  
end
% %TJP: TO conserve memory do not store jac and x
if storeJac;
    cds.pJac = jac;
    cds.pJacX = x;
end

% % newjac = jac;
% % save newjac newjac
% goodjac = jac;
% save goodjac goodjac
% pause

%cds.pJacX = x;

%SD:calculates jacobian
