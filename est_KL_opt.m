function est = est_KL_opt(samp_p,samp_q,k)
% est = est_KL_opt(samp_p,samp_q,k,f)
%
% Minimax optimal estimator proposed by Bu, Zou, Liang and Veeravalli.
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
%------ k: The true alphabet size of the distribution pair p and q.
%------ f: The true desity ratio upper bound of the distribution pair p and q.
%
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
    
    if ~(isequal(wid_p, wid_q))
        error('Number of trials are not equal.');
    end     
    c_0=1.2;
    
    order = floor(c_0*log(k));  %order of the polynomial approximation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    persistent coeffs;
    if isempty(coeffs)
        load poly_coeff_KL.mat coeffs; %load polynomial coefficients
    end
    coeff = coeffs(order,:);  
    coeff(isnan(coeff))=[];
    
    coeff(1)=[];
    coeff=-coeff;
    
    
    est = zeros(1,wid_q);
    for i=1:wid_p
        temp_p = samp_p(:,i);
        temp_q = samp_q(:,i);
        edge_p = unique(temp_p);
        hist_p = histc(temp_p,edge_p);
        hist_q = histc(temp_q(ismember(temp_q,edge_p)),edge_p);

        c_1 = 2*max(set_c1(hist_p,k), 1/(1.9*log(k)));
        c_2 = c_1/2 ;
         
        thres = c_2 * log(k);
        
        hist_p_poly = hist_p(hist_p<thres);
        hist_p_plugin = hist_p(hist_p>=thres);
%         temp = entro_poly(hist_p_poly, k, m, coeff, c_1);
        entrop_est = sum(entro_poly(hist_p_poly, k, m, coeff, c_1));
        entrop_est = entrop_est + sum((hist_p_plugin/m).*log(hist_p_plugin/m)-1/(2*m));
        
%         c_1 = 2*max(set_c1(hist_q,k), 1/(1.9*log(k)));
%         c_2 = c_1/2;         
%         thres = c_2 * log(k);
        
        hist_q_poly = hist_q(hist_q<thres);
        hist_q_plugin = hist_q(hist_q>=thres);
        plogq_est = sum((hist_p(hist_q<thres)/m).*plogq_poly(hist_q_poly, k, n, coeff, c_1));
        plogq_est = plogq_est + sum((hist_p(hist_q>=thres)/m).*(log((hist_q_plugin+1)/n)-1./(2*(hist_q_plugin+1))));

        est(i) = (entrop_est - plogq_est);
        est(i) = max(est(i),0);
    end

end

function output = entro_poly(x, k,m, g_coeff, c_1)
% Compute the polynomial for entropy part p*log(p)
    order = length(g_coeff);
    index = 0:order-1;
    moment = repmat(gamma(x+1),1,order)./gamma(repmat(x+1,1,order)-repmat(index+1,length(x),1));
    output = 1/m*(moment*(g_coeff./((c_1*log(k)).^index))'-log(m/(c_1*log(k)))*x);
end


function output = plogq_poly(x, k, n, g_coeff, c_1)
% Compute the polynomial for p*log(q)
    order = length(g_coeff);
    index = 0:order-1;
    moment = repmat(gamma(x+1),1,order)./gamma(repmat(x+1,1,order)-repmat(index,length(x),1));
    output = moment*(g_coeff./((c_1*log(k)).^index))'-log(n/(c_1*log(k)));
end
 
function output = set_c1( hist, n)
% Set the optimal constant
    V1 = [0.3303 0.4679];     
    V2 = [-0.530556484842359,1.09787328176926,0.184831781602259];   
    f = histc(hist, unique(hist));
    f1nonzero = f(1,:) > 0;
    if  any(f1nonzero)
        if n < 200
            c_1(f1nonzero) = polyval(V1, log(n./f(1,f1nonzero)));   
        else
            n2f1_small = f1nonzero & log(n./f(1,:)) <= 1.5;
            n2f1_large = f1nonzero & log(n./f(1,:)) > 1.5;
            c_1(n2f1_small) = polyval(V2, log(n./f(1,n2f1_small)));  
            c_1(n2f1_large) = polyval(V1, log(n./f(1,n2f1_large)));  
        end
    end
    output = c_1(f1nonzero); 
end
