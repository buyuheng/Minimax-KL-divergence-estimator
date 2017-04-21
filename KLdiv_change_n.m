clc;
close all; 
clear;  


C = 3;   % Constant before the number of samples from q
num = 20;   % Number of points sampled in the plot
mc_times = 10;  % Total number of Monte-Carlo trials for each alphabet size 

record_k = 10^4;    % Alphabet size

record_m = ceil(logspace(3, 5, num)); % Number of samples from p



for iter = num:-1:1 
    k = record_k
    m = record_m(iter)

%%  First distribution pair
%     f=5;    % Density ratio
%     dist_p = ones(k,1)/k;    
%     dist_q = dist_p/f;
%     dist_q(1) = 1-sum(dist_q(2:end));
    
%%  Zipf distribution pair

    alpha = 1;      
    beta = 0.6;
    
    dist_p = 1./[1: k].^alpha;
    dist_p = dist_p/(sum(dist_p));
    
    dist_q = 1./[1: k].^beta;
    dist_q = dist_q/(sum(dist_q));
    f = max(dist_p./dist_q);     % Compute density ratio

%%    
     
    n = ceil(C*f*record_m(iter))    % Number of samples from q
    
    true_KL(iter) = KL_divergence_true(dist_p,dist_q);     
    samp_p = randsmpl(dist_p, m, mc_times, 'int32');    
    samp_q = randsmpl(dist_q, n, mc_times, 'int32');
   
    tic
    record_Aplugin = est_KL_Aplugin(samp_p,samp_q);  
    toc
    tic
    record_opt = est_KL_opt(samp_p,samp_q,k);  
    toc

    Aplugin_err(iter) = sqrt(mean((record_Aplugin - true_KL(iter)).^2));  
    opt_err(iter) = sqrt(mean((record_opt - true_KL(iter)).^2));
end 

figure(1) 
semilogx(record_m, opt_err,'b-s','LineWidth',2,'MarkerFaceColor','b'); 
hold on;
semilogx(record_m, Aplugin_err,'m-.o','LineWidth',2,'MarkerFaceColor','r'); 
hold on;


legend('BZLV A-plugin','BZLV opt');
xlabel('m');
ylabel('RMSE/nats');
