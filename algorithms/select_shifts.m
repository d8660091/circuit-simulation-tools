function [w, q]= select_shifts(H,which);
%
%     Compute the eigenvalues of H and sort according to
%
%     which = 'LR', 'SR', 'LM','LA','SA'
%
%     Input:  H     -- a real square upper Hessenberg matrix
%
%             which -- a string indicating which eigenvalues are
%                      wanted
%
%     Output: w     -- eigenvalues of H sorted according to which
%                      most wanted to least wanted order
%             q     -- absolute value of last component of eigenvectors
%                      of H in same order as the eigenvalues in w
%
%
      [q,w] = eig(H);
      w = diag(w); 

      m = length(w);
      q = q(m,:);
%
%     select filter mechanism by activating appropriate choice below
%
      switch which


      case {'LR', 'LA'} %   sort for largest real part 
%
         [s,ir] = sort(-real(w));

%        shifts are smallest real part

      case {'SR', 'SA'} %   sort for smallest real part 
%
         [s,ir] = sort(real(w));

%        shifts are largest real part

      case {'LM'}       %   sort for largest magnitude 
%
         [s,ir] = sort(-abs(w));

      end
      w = w(ir);
      q = abs(q(ir));
 
