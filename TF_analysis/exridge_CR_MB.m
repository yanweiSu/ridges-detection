function [c,e] = exridge_CR_MB(Tx,lambda,beta,jump,C)
% exridge : extracts  the ridge curve by maximising some energy. 
%   The result is an approximation computd by a greedy algorithm.
%   The algorithm uses several random initializations and then
%   forward/backward search.
%
% INPUTS:   
%   Tx : synchrosqueezed transform
%   lambda : the parameter associated with first derivative penalization (usually equal to zero)
%   beta   : the parameter associated with the second derivative penalization (usually equal to zero)
%   jump   : the parameter associated with the chirp rate (estimated via
%   synchrosqueezing)
%    C     : real constant
% OUTPUTS:  
%   c : vector containing the indexes of the ridge. length :size(Tx,2).
%   e : energy of the returned ridge, scalar.


%Et = log(abs(Tx)+eps^0.25);
Et     = abs(Tx).^2;
[na,N] = size(Tx);

% Parameters of the algorithm.
da = jump; 
ng = min(60,floor(N/8)); % Number of random initializations.

ccur = zeros(1,N);
c = ccur;    
e = -Inf;


for k = floor(linspace(N/(ng+1),N-N/(ng+1),ng))
    [ecur,idx] = max(Et(:,k));        
    ccur(k) = idx;
%     Iq = max(1,idx-da(idx,k)):min(na,idx+da(idx,k));
%     [ecur,idx] = max(Et(Iq,k-1));        
%     ccur(k-1) = idx;
   
%     idx = ccur(k);
    % forward step
    for b=k+1:N
     I = intersect(ccur(b-1)+da(ccur(b-1),b-1)-C:ccur(b-1)+da(ccur(b-1),b-1)+C,1:na);
     if isempty(I)
         if ccur(b-1)+da(ccur(b-1),b-1)>na
             I = na;
         else
             I = 1;
         end;
     end;
     [etmp,idx] = max(Et(I,b));
     ccur(b)=idx+I(1)-1;
     ecur = ecur + etmp;
    end

    % backward step
    idx = ccur(k);
    for b=k-1:-1:1
        I = intersect(ccur(b+1)-da(ccur(b+1),b+1)-C:ccur(b+1)-da(ccur(b+1),b+1)+C,1:na);
         if isempty(I)
         if ccur(b+1)-da(ccur(b+1),b+1)>na
             I = na;
         else
             I = 1;
         end;
        end;
        [etmp,idx] = max(Et(I,b));
        ccur(b)=idx+I(1)-1;
        ecur = ecur + etmp;
    end
    if ecur> e
        e = ecur;
        c = ccur;
    end
end

end
