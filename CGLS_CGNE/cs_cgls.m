function [x, dnrm, snrm, enrm, obj2] = cs_cgls(A, b, P, C, F, ...
                                    m_outer, tol, m_iter, eps, lambda, x_mod)

[m,n] = size(A);
g = size(F,1);

   x0 = ones(n,1); 
   dif = norm( x0 );
   it = 0;

   dnrm = [];
   snrm = [];
   enrm = [];
   obj2 = [];

   while ( dif > tol & it < m_outer )

      it = it + 1;
%       w = (abs(F*x0)+eps).^(p-2);
%       PP = spdiags(1./sqrt(w),0,n,n);
      if it == 1
         W = speye(g,g); 
      else
         W = spdiags(1./(abs(F*x0)+1e-15),0,g,g);
      end
      R = F'*W*F;
    
      [x, res1, res2, err, obj] = cgls( A, b, x0, P, C, R, ...
                                m_iter, eps, lambda, x_mod,1 );
        
%         [x, res1, res2, err, obj] = cgls( A, b, zeros(n,1), P, C, R, ...
%                                 m_iter, eps, lambda, x_mod,1 );
                            
%       dif = norm( x - x0 );
%       sol = sum( abs(x) );
%       err = norm( x - x_mod );
      
      x0 = x(:,end);

%       dnrm = [dnrm, dif];
%       snrm = [snrm, sol];
%       enrm = [enrm, err];

%       dnrm = [dnrm, res1];
%       snrm = [snrm, res2];
%       enrm = [enrm, err];
%       obj2 = [obj2, obj];

      dnrm = [dnrm, [res1 NaN*zeros(1,m_iter-length(err))]];
      snrm = [snrm, [res2 NaN*zeros(1,m_iter-length(err))]];
      enrm = [enrm, [err NaN*zeros(1,m_iter-length(err))]];
      
      for k = 1:size(x,2)
        obj(k) = 1/2*abs((b-A*x(:,k))'*(C\(b-A*x(:,k))))+1/2*lambda*sum(abs(F*x(:,k)));
      end

      
      obj2 = [obj2, [obj NaN*zeros(1,m_iter-length(err))]];

   end

end
