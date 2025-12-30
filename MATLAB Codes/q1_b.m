fs = 1000;                  
Ts = 1/fs;
t = 0:Ts:1;
N = length(t);

BW = 40;                
N_iter = 150;             
p_values = 0.01:0.01:0.10;

signals = {
    sin(2*pi*5*t);                                  
    sin(2*pi*5*t) + 0.6*sin(2*pi*15*t);            
    sin(2*pi*5*t) + (1/3)*sin(2*pi*15*t) + (1/5)*sin(2*pi*25*t)
};

for s = 1:3
    x = signals{s};
    errors = zeros(size(p_values));
    x_plot = [];
    x_gp_plot = [];
    mask_plot = [];

    for pp = 1:length(p_values)
        p = p_values(pp);
        mask = rand(size(x)) > p;
        x_tilde = x .* mask;
        known_idx = find(mask==1);
        x_init = interp1(known_idx, x(known_idx), 1:N, 'pchip', 'extrap');
        x_gp = x_init;
        L = N;
        f = (-L/2:L/2-1)*(fs/L);
        H = abs(f) <= BW;

        for k = 1:N_iter
            Xf = fftshift(fft(x_gp));
            Xf = Xf .* H;
            x_band = real(ifft(ifftshift(Xf)));
            x_band(mask==1) = x(mask==1);
            x_gp = x_band;
        end

        errors(pp) = mean(abs(x - x_gp));

        if abs(p - 0.05) < 1e-6
            x_plot = x;
            x_gp_plot = x_gp;
            mask_plot = mask;
        end
    end

    figure(s*2-1);
    subplot(3,1,1);
    stem(x_plot,'b');
    grid on;
    title('Original Signal');

    subplot(3,1,2);
    stem(find(mask_plot==1), x_plot(mask_plot==1),'g'); hold on;
    stem(find(mask_plot==0), zeros(sum(mask_plot==0),1),'r');
    title('Available and Missing Samples');

    subplot(3,1,3);
    stem(x_plot,'k--'); hold on;
    stem(x_gp_plot,'m');
    grid on;
    title('Reconstructed Signal');
    legend('Original','Reconstructed');

    figure(s*2);
    plot(p_values, errors,'-o','LineWidth',2);
    grid on;
    xlabel('p'); ylabel('MAE');
    title('Error vs Missing Probability');
end