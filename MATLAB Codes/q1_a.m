clc; clear; close all;

Fs   = 1000;           
Ts   = 1/Fs;
tmax = 1;              
t    = 0:Ts:tmax;   
Nideal = length(t);

signals = {
    sin(2*pi*10*t) + 0.7*sin(2*pi*25*t);           
    sin(2*pi*5*t) + 0.6*sin(2*pi*15*t);            
    sin(2*pi*2*t) + 0.5*sin(2*pi*30*t) + 0.4*sin(2*pi*50*t);  
};

K_values = 1:4;
Final_MAE = zeros(3, length(K_values));  

for sig_idx = 1:3
    xtrue = signals{sig_idx};   

    L  = 10;
    Tfine = Ts/L;
    tfine = 0:Tfine:tmax;
    Nfine = length(tfine);

    switch sig_idx
        case 1
            xtrue_fine = sin(2*pi*10*tfine) + 0.7*sin(2*pi*25*tfine);
        case 2
            xtrue_fine = sin(2*pi*5*tfine) + 0.6*sin(2*pi*15*tfine);
        case 3
            xtrue_fine = sin(2*pi*2*tfine) + 0.5*sin(2*pi*30*tfine) + 0.4*sin(2*pi*50*tfine);
    end

    for kk = 1:length(K_values)
        K = K_values(kk);
        Delta = Ts/10;          
        n  = 0:floor(tmax/Ts);
        kn = randi([-K K], size(n));     
        tn = n*Ts + kn*Delta;            
        valid = (tn >= 0) & (tn <= tmax);
        tn = tn(valid);
        n  = n(valid);

        switch sig_idx
            case 1
                x_hat = sin(2*pi*10*tn) + 0.7*sin(2*pi*25*tn);
            case 2
                x_hat = sin(2*pi*5*tn) + 0.6*sin(2*pi*15*tn);
            case 3
                x_hat = sin(2*pi*2*tn) + 0.5*sin(2*pi*30*tn) + 0.4*sin(2*pi*50*tn);
        end

        x_rec = interp1(tn, x_hat, tfine, 'pchip', 'extrap'); % piecewise cubic hermite polynomial

        Nit = 100;
        cutoff = 120; 

        for it = 1:Nit
            X = fft(x_rec);
            X = fftshift(X);

            fshift = (-Nfine/2 : Nfine/2-1)*(1/(Nfine*Tfine));
            mask = abs(fshift) <= cutoff;

            X(~mask) = 0;
            X = ifftshift(X);

            x_bl = real(ifft(X));

            idx = round(tn/Tfine) + 1;
            idx = max(1, min(idx, Nfine));

            x_bl(idx) = x_hat;
            x_rec = x_bl;
        end

        idx_ideal = round((0:Nideal-1)*Ts/Tfine) + 1;
        idx_ideal = max(1, min(idx_ideal, Nfine));

        x_rec_samples = x_rec(idx_ideal);
        x_true_samples = xtrue;

        abs_err_samples = abs(x_rec_samples - x_true_samples);
        Final_MAE(sig_idx, kk) = mean(abs_err_samples);

        figure;
        sgtitle(sprintf("Signal %d - Results for K = %d", sig_idx, K));

        subplot(3,1,1);
        plot(tfine, xtrue_fine,'b','LineWidth',1.1); hold on;
        stem(tn, x_hat,'r','Marker','none');
        title("True signal and jittered samples"); xlabel("t");
        legend("True", "Jittered");

        subplot(3,1,2);
        plot(t, abs_err_samples,'b','LineWidth',1.2);
        title("Absolute Error at Ideal Sampling Instants");
        xlabel("Time"); ylabel("|error|");

        subplot(3,1,3);
        stem(t, x_rec_samples,'b'); hold on;
        stem(t, x_true_samples,'r--');
        title("Reconstructed vs True (Ideal Sample Instants)");
        xlabel("t");
        legend("Reconstructed","True");
    end
end

figure;
hold on;
colors = ['r', 'g', 'b'];
for sig_idx = 1:3
    plot(K_values, Final_MAE(sig_idx, :), '-o', 'MarkerFaceColor', colors(sig_idx));
end
grid on;
title("Final Reconstruction MAE vs Jitter Level K");
xlabel("Jitter K");
ylabel("Final MAE");
legend('Signal 1', 'Signal 2', 'Signal 3');