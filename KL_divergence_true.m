function D = KL_divergence_true(p,q)
%KL_divergence_true  computes Kullback-Leiblier divergence D(p||q) in nats 
%for the input discrete distribution pair.

% This function returns a scalar divergence when the input distribution p and q are  
% vectors of probability masses, or returns in a row vector the columnwise 
% entropies of the input probability matrix p and q.

% Error-check of the input distribution pair
p0 = p(:);
q0 = q(:);

if (any(size(p) ~= size(q)))
    error('p and q must be equal sizes.');
end

if any(imag(p0)) || any(isinf(p0)) || any(isnan(p0)) || any(p0<0) || any(p0>1)
    error('The probability elements of p must be real numbers between 0 and 1.');
elseif any(abs(sum(p)-1) > sqrt(eps))
    error('Sum of the probability elements of p must equal 1.');
end

if any(imag(q0)) || any(isinf(q0)) || any(isnan(q0)) || any(q0<=0) || any(q0>1)
    error('The probability elements of q must be real numbers between 0 and 1.');
elseif any(abs(sum(p)-1) > sqrt(eps))
    error('Sum of the probability elements of q must equal 1.');
end

D = sum(p.*log((p./q)));  
end