function [x,v,i] = newtcorr(x0, v0)
%
% Newton corrections, internal routine
%
global cds 
%global doJacRichardExtrapolation;
x = real(x0);
v = real(v0);

%display(['any complex x0 terms: ',num2str(any(~real(x0) & x0~=0 ))]);
%display(['any complex v0 terms: ',num2str(any(~real(v0) & v0~=0))]);

% TJP
clear x0 v0;
lasterror('reset');

RT = [];
RT(cds.ndim) = 1;R=RT';
clear RT;
STPMX=100;

%doJacRichardExtrapolation=false;

% set some iteration params and adjust for conversion of initial est
maxCorrIters = cds.options.MaxCorrIters;
maxints = cds.options.MaxSecantIters; 
if sum(v)==1;
    maxints=1;
    maxCorrIters = 2*maxCorrIters;
%    doJacRichardExtrapolation=false;
end
%doJacRichardExtrapolation=true;
%try
    for i = 1:maxCorrIters

      % Calc. Jacobian
      B = [contjac(x,false); v'];       
          
      Q = getQ(cds.curve_func, x);
      
      % repeat twice with same jacobian, calculating the jacobian is usually 
      % a lot more expensive than solving a system
      its=0;
      prev_normQ = 0; normQ=1e29;
      prev_normdx = 0; normdx=1e29;
      
      doNextIteration=true;      
      while doNextIteration;
        if isnan(norm(Q)) || isinf(norm(Q))
            x = [];
            v = [];
            return;
        end
        
        if cds.options.MoorePenrose
          lastwarn('');
          
          D = B\[Q R];          
          if ~isempty(lastwarn)
            x = [];
            v = [];
            return;
          end
          v = D(:,2);
          v = v/norm(v);                    
          
          dx = D(:,1);
        else
          dx = B\Q;          
        end
        
        % Do a line search step
        stpmax = STPMX * max(norm(x),size(x,1));
        [xnew, Qnew, lambdaReduced, err] = lnsrch( cds.curve_func, x, Q , B, -dx, stpmax);
        
        its=its+1;
        prev_normdx = normdx;
        prev_normQ = normQ;
        normdx = norm(x-xnew);
        normQ = norm(Qnew);   

        if its==1 && ~err && normdx < cds.options.VarTolerance && normQ < cds.options.FunTolerance    
            % Convergence achieved.
            % To avoid errors in calc of v, a secant approx to v cannot be
            % used to calc final v.
            
            x=xnew;
            %display(['   ...t=',num2str(its),'(FINAL), norm(dx) =',num2str(normdx,6),', norm(Q) =',num2str(normQ,6) ...
            %   ', max(abs(dx)) =',num2str(max(abs(dx)),6),', max(abs(Q)) =',num2str(max(abs(Q)),6) ]);
            
            clear Q B D dx;
            B = contjac(x,false);
            v = ([B ; v']\R);
            v = v/norm(v);
            
            return;
            
        elseif ~err && prev_normQ>normQ && prev_normdx > normdx;
            % Only update x if solution converging                        
            
            % Update x
            x=xnew;            
            
            % If converged to below both criteria, then set the next
            % iteration to a non-secant step.
            if normdx < cds.options.VarTolerance && normQ < cds.options.FunTolerance
                doNextIteration=false;   
            else
                doNextIteration=true;                                    
            end           
            
            %display(['   ...t=',num2str(its),', norm(dx) =',num2str(normdx,6),', norm(Q) =',num2str(normQ,6) ...
            %   ', max(abs(dx)) =',num2str(max(abs(dx)),6),', max(abs(Q)) =',num2str(max(abs(Q)),6) ]);
                       
        else
            %if err || prev_normQ <= normQ || prev_normdx <= normdx;
            % if solution is NOT converging do not update x and replace v
            % with v_prev which is for the non updates x
            
            if ~err ;
                x=xnew;            
            else
                display('...Reducing param step size as error detected (possible line search error)');
                x = [];
                v = [];
                return;
            end
            maxints=1;
            doNextIteration=false;
            %display(['   ...Recalc B as non-convergence at t=',num2str(its),', norm(dx) =',num2str(normdx,6),', norm(Q) =',num2str(normQ,6) ...
            %   ', max(abs(dx)) =',num2str(max(abs(dx)),6),', max(abs(Q)) =',num2str(max(abs(Q)),6) ]);
        end
        
        if its>=maxints;
            doNextIteration=false;
        end
        if lambdaReduced;
            maxints=1;
            doNextIteration=false;
            %doJacRichardExtrapolation=true;
        end
        clear D dx;
      end  
      % NOTE: this is vital so that the memory requirements do not exceed RAM  
      clear B;
    end
    
    display('... Newton convergence not achieved. Step size to be reduced and convergence re-attempted');    
    x = [];
    v = [];
    
%   catch  
%        x = [];
%        v = [];    
%        display('...Reducing param step size as model error detected');
%        return;     
%   end    
%SD:Newton corrections pal/mp
end

% Get fucntion for which a root is required
%------------------------------------------------------
function Q = getQ(func, x)
    Q = [feval(func, x); 0];
end
%------------------------------------------------------    

% Line search and bracketing 
% from Press et al 2007 Numerical Recipes section 9.7.1 pages478-780
% Note:
% g = jacobian matrix B
% p = detla x (ie newton step)
%------------------------------------------------------
function [x,f, lambdaReduced, err] = lnsrch(func, xold, fold, B, p, stpmax)
    ALF = 1e-5; % constant ensuring sufficient decrease in func value
    TOLX = eps;
    
    % some declarations
    n = size(xold,1);    
    err = false;
    lambdaReduced = false;
    p_norm = norm(p);    
    f_min_old = 0.5*dot(fold,fold);
    
    % scale p is attempt too big
    if p_norm > stpmax;        
        p = p .* stpmax/p_norm;
    end
    
    % Slope is reculated from B when F is updated
    slope = sum( (B'*fold) .*p);    
    if slope >0        
        % Roundoff problem in lnsrch
        display(['Warning: Rounding error: lnsrch slope= ', num2str(slope,6), ' which is >0']);
        err = true;
        return;
    end
    
    % compute lambda_min
    test = max(abs(p)./max(abs(xold), ones(n,1)) );
    alamin = TOLX/test;
    
    % get exact alamin
    temp = alamin + 1;
    alamin = temp - 1;
    
    % initialise lambda
    alam=1.0;                           % always try full newton step first
%    int = 0;
    
%     best_f = fold;
%     best_f_min = f_min_old;
%     best_x = xold;
%     best_int = 0;
%     best_lambda = 0;
    
    % start main loop
    while(1)
       
       % check if alam is below machine precision. If so return prior est of x and f 
       if alam<=alamin;

%           if best_int==0;
               display(['Error: lambda constrained to lambda_min, ie ', num2str(alam,6), ' < ', num2str(alamin,6),'. No reduction in f achieved.' ]);
               err=true;
%           else            
%            display(['Warning: lambda constrained to lambda_min, ie ', num2str(alam,6), ' < ', num2str(alamin,6),'. Returning x and f with lambda=1.' ]);   
%            display(['....     lambda=1, best interation: ', num2str(best_lambda,6), ' , ', num2str(best_int,6) ]);            
%            x = x_initial;
%            f = f_initial;                        
%           end
           return
       end
        
       x = xold + alam.*p;
       f = getQ(func,x);
       f_min = 0.5*dot(f,f);
       %int =int+1;
       
       % Slope is reculated from B when F is updated
       slope = sum( (B'*f) .*p);    
       
       if alam==1;
            f_initial = f;
            %f_min = f_min;
            x_initial = x;           
            %best_int = int;
            %best_lambda = alam;           
       end
       
       % convergence of delta x
       % sufficient function decrease    
       if f_min <= f_min_old + ALF * alam * slope;
           return;
       %backtrack    
       else
           if alam==1.0;
               % first backtrack
               tmplam = -slope/(2 * (f_min - f_min_old - slope) );
           else
               % subsequant backtracks
               %display('...... reducing newton step (lambda<1)');
               lambdaReduced=true;               
               rhs1 = f_min - f_min_old - alam * slope;
               rhs2 = f_min2 - f_min_old - alam2 * slope;
               a = ( rhs1/alam^2 - rhs2/alam2^2)/(alam - alam2);
               b = ( -alam2 * rhs1 / alam^2 + alam * rhs2/alam2^2)/( alam - alam2);
               
               if a==0;
                   tmplam = -slope/(2*b);
               else
                   disc = b^2 - 3*a*slope;
                   if disc<0;
                       tmplam = 0.5*alam;
                   elseif b<=0;
                       tmplam = (-b + sqrt(disc))/(3*a);
                   else
                       tmplam = -slope/(b + sqrt(disc));
                   end
               end
               if tmplam>0.5*alam;
                   tmplam=0.5*alam;         % try lambda = 0.5 lambda_1
               end
           end
       end
       
       alam2=alam;
       f_min2=f_min;       
       alam = max(tmplam, 0.1*alam);        % try lambda = 0.1 lambda_1
       
       % try again
    end
end
%------------------------------------------------------