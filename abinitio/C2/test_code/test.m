SNRs = [1/2,1/4,1/8,1/16];
nTests             = numel(SNRs);
res_cl_detec_rate  = zeros(4,nTests);
res_rank1_rate     = zeros(4,nTests);
res_err_in_degrees = zeros(4,nTests);

nImages = 100;

for k=1:nTests
    snr = SNRs(k);
    
    [cl_detec_rate,rank1_rate,err_in_degrees] = do_test(nImages,snr,false,false);
    res_cl_detec_rate(1,k) = cl_detec_rate;
    res_rank1_rate(1,k) = rank1_rate;
    res_err_in_degrees(1,k) = median(err_in_degrees);
    
    
    [cl_detec_rate,rank1_rate,err_in_degrees] = do_test(nImages,snr,true,false);
    res_cl_detec_rate(2,k) = cl_detec_rate;
    res_rank1_rate(2,k) = rank1_rate;
    res_err_in_degrees(2,k) = median(err_in_degrees);
    
    
    [cl_detec_rate,rank1_rate,err_in_degrees] = do_test(nImages,snr,false,true);
    res_cl_detec_rate(3,k) = cl_detec_rate;
    res_rank1_rate(3,k) = rank1_rate;
    res_err_in_degrees(3,k) = median(err_in_degrees);
    
    
    [cl_detec_rate,rank1_rate,err_in_degrees] = do_test(nImages,snr,true,true);
    res_cl_detec_rate(4,k) = cl_detec_rate;
    res_rank1_rate(4,k) = rank1_rate;
    res_err_in_degrees(4,k) = median(err_in_degrees); 
end

figure; 
plot(res_err_in_degrees(1,:)); hold on; 
plot(res_err_in_degrees(2,:)); 
plot(res_err_in_degrees(3,:)); 
plot(res_err_in_degrees(4,:)); 

hold off; 

legend('false,false','true,false','false,true','true,true')
