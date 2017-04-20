clc;
close all; 
clear;  


C_1 = 2;  % Constant before the number of samples from p
C_2 = 1;  % Constant before the number of samples from q
num = 20; % Number of points sampled in the plot
mc_times = 10;  % Total number of Monte-Carlo trials for each alphabet size 

record_k = ceil(logspace(3,6, num)); % Alphabet size

record_m = ceil(C_1*record_k./log(record_k)); % Number of samples from p
record_n = ceil(C_2*record_k./log(record_k)); 


 
for iter = num:-1:1 
    k = record_k(iter)
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
   
    n = ceil(f*record_n(iter))   % Number of samples from q
    
    true_KL(iter) = KL_divergence_true(dist_p,dist_q);     
    samp_p = randsmpl(dist_p, m, mc_times, 'int32');    
    samp_q = randsmpl(dist_q, n, mc_times, 'int32');
    
    tic
    record_Aplugin = est_KL_Aplugin(samp_p,samp_q);  
    toc
    tic
    record_Mplugin = est_KL_Mplugin(samp_p,samp_q);  
    toc
    tic
    record_opt = est_KL_opt(samp_p,samp_q,k);  
    toc

    Aplugin_err(iter) = sqrt(mean((record_Aplugin - true_KL(iter)).^2));  
    Mplugin_err(iter) = sqrt(mean((record_Mplugin - true_KL(iter)).^2)); 
    opt_err(iter) = sqrt(mean((record_opt - true_KL(iter)).^2));
end 

figure(1) 
semilogx(record_k, opt_err,'b-s','LineWidth',2,'MarkerFaceColor','b'); 
hold on;
semilogx(record_k, Aplugin_err,'m-.o','LineWidth',2,'MarkerFaceColor','r'); 
hold on;
semilogx(record_k, Mplugin_err,'r-*','LineWidth',2,'MarkerFaceColor','k'); 
hold on;

legend('BZLV opt','BZLV A-plugin','HJW M-plugin');
xlabel('k');
ylabel('RMSE/nats');
