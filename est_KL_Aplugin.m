function est = est_KL_Aplugin(samp_p,samp_q)
% est = est_KL_Aplugin(samp_p,samp_q)
%
% Augmented plug-in estimation of KL divergence (in nats) of the input samples
%
% This function returns a scalar estimation of KL divergence D(p||q) when samp_p 
% and samp_q are vectors, or returns a (row-) vector consisting of the estimation 
% of KL divergence of each column of samples when samp_p and samp_q are matrix.
%
% Input:
% ----- samp_p: a vector or matrix which can only contain integers. The input
%             data type can be any interger classes such as uint8/int8/
%             uint16/int16/uint32/int32/uint64/int64, or floating-point 
%             such as single/double. 
% ----- samp_q: same conditions as sampP. Must have the same number of
%              columns if a matrix.
% Output:
% ----- est: the KL divergence (in nats) of the input vector or that of each 
%            column of the input matrix. The output data type is double. 


% Error-check
    if ~(isequal(samp_p, fix(samp_p)) && isequal(samp_q, fix(samp_q)))
        error('Input sample must only contain integers.');
    end

    if isrow(samp_p)
        samp_p = samp_p.';
    end
    [m, wid_p] = size(samp_p);

    if isrow(samp_q)
        samp_q = samp_q.';
    end
    [n, wid_q] = size(samp_q);
    
    if ~(isequal(wid_p, wid_q))
        error('Number of trials are not equal.');
    end  
    
    est = zeros(1,wid_q);
    for i=1:wid_q
        temp_p = samp_p(:,i);
        temp_q = samp_q(:,i);
        edge_p = unique(temp_p);
        hist_p = histc(temp_p,edge_p);
        hist_q = histc(temp_q(ismember(temp_q,edge_p)),edge_p);
        est(i) = sum((hist_p/m).*log((hist_p/m)./((hist_q+1)/n)));
    end

end