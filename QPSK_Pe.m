% Initializations
snr_db= -4:1:10; 
num_sim= 10e5;
num_snr_db=length(snr_db);
A=1;
Pe_analytical = zeros(1,num_snr_db);
%probability of error using analytical method 
for i= 1:num_snr_db
    snr_lin = 10^(0.1*snr_db(i));
    Pe_analytical(i) = 2*qfunc(sqrt(snr_lin))*(1-0.5*qfunc(sqrt(snr_lin)));
end
%probability of error using Monte Carlo Simulation
for i=1:num_snr_db
    snr_lin = 10^(0.1*snr_db(i));
    count = 0;
    for b = 1:1:num_sim
        n= (sqrt(1/(2*snr_lin))*randn(1)+1j*sqrt(1/(2*snr_lin))*randn(1)); 
        a=rand;
        if a<=0.25
            s=A;
        end
        if (0.25 < a)&&(a <= 0.5)
            s=1j*A;
        end
        if (0.5 < a)&&(a <= 0.75)
            s=-A;
        end
        if (0.75 < a)&&(a <= 1)
            s = -1j*A;
        end
        r= s+n;
        msg = [A,0+1j*A,-A,-1j*A];
        for j=1:1:4
            msgl = msg(j);
            p(j)= (norm(r-msgl))^2;
        end
        [min_value,index] = min(p);
        if(msg(index)~=s)
            count = count+1;
        end
    end
    Pe_simulation(i) = count/num_sim;
end
semilogy (snr_db,Pe_analytical,'r');
hold on;

semilogy (snr_db,Pe_simulation,'b');
xlabel('SNR(dB)');
ylabel('Pe');
title('Pe as a function of snr');
legend('Analytical','Monte Carlo Simulation');


