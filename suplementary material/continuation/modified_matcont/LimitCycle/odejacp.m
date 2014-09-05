function jac=odejacp(t,u,p)

global lds cds  

u_typical = 1000; 
jac=repmat(0,size(u,1),size(p,1));
doJacRichardExtrapolation=true;
if doJacRichardExtrapolation
   for i=1:length(p) 
      
        % get del into an extact number       
        del = cds.options.Increment;         

        if isreal(cds.options.Increment);
            del = max(eps,del);
            temp = p{i} + del;
            del = temp - p{i};          
        else
            del = max(1i*eps(u_typical),del);
            temp = 1i*u_typical + del;
            del = temp - 1i*u_typical;      
        end

       jac(:,i) = Ridders(lds, u,  t, p, i, del);   
       
   end
else
    for i=1:length(p)

        % get del into an extact number       
        del = cds.options.Increment;         

        if isreal(cds.options.Increment);
            del = max(eps,del);
            temp = p{i} + del;
            del = temp - p{i};          
        else
            del = max(1i*eps(u_typical),del);   
            temp = 1i*u_typical + del;
            del = temp - 1i*u_typical;      
        end

        if del==0;
          warning('Parameter Jacobian del==0. May be due to machine precision rounding');
        end

        % calc centre weighted finite diff
        p1 = p; p1{i} = p1{i}-del;
        p2 = p; p2{i} = p2{i}+del;   
        jac(:,i) = 1./(2.*del) * ( feval(lds.func, t, u, p2{:}) - feval(lds.func, t, u, p1{:}) );       
    end
    
end
end

function [dFdu, err] = Ridders(lds, u,  t, params, param_index, h)
% Ridder's method for polynomial extrapolation of centre weighted finite
% diff derivatives.
% From Press et al , 2007, Numerical Recipes, p 230 secion 5.7, editon 3.
    ntab = 10;
    con = 1.4;
    con2 = con^2;
    safe=2;
    big = 1e29;
    a = cell(ntab, ntab);

    if h==0;
        display('h must be >0')
        error(' ');
    end

    hh=h;

    % get est for h
    p1 = params; p1{param_index} = p1{param_index} - hh;
    p2 = params; p2{param_index} = p2{param_index} + hh;   
    a{1,1} = 1./(2*hh).*(feval(lds.func, t, u, p2{:})-feval(lds.func, t, u, p1{:}));          

    err=big;
    ii=2;
    
    while ii<=ntab
       
        hh = hh/con;
        
        % try smaller step of h
        p1 = params; p1{param_index} = p1{param_index} - hh;
        p2 = params; p2{param_index} = p2{param_index} + hh;   
        a{1,ii} = 1./(2*hh).*(feval(lds.func, t, u, p2{:})-feval(lds.func, t, u, p1{:}));          
               
        if norm(a{1,ii})==0;
           dFdu = a{1,ii};
           return
        end
        fac=con2;
        
        % compute extrapolations
        jj=2;
        while jj<=ii
           a{jj,ii}  = (a{jj-1,ii}.*fac - a{jj-1,ii-1}) ./ (fac - 1);
           fac = con2*fac;
           errt = norm( max( abs(a{jj,ii} - a{jj-1,ii}), abs(a{jj,ii} - a{jj-1,ii-1})) );
           
           
           if errt < err;
               err=errt;
               dFdu = a{jj,ii};
               
               if any(any(imag(dFdu)~=0));
                   display('dFdu term is imag within Ridders Algorithm');
                   error(' ');
               end
           end
           jj=jj+1;
        end
        
        if norm(abs(a{ii,ii} - a{ii-1,ii-1})) > safe*err;
           return;
        end
        ii=ii+1;
    end    
end
