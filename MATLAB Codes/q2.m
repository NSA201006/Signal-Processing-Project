fileList = {'File1.wav', 'File2.wav', 'File3.wav','File4.wav'}; 

for fileIdx = 1:numel(fileList)
    [audio, fs] = audioread(fileList{fileIdx});
    audio = mean(audio,2); 
    %Performing Stft and then finding the spectrogram analysis
    window_size = 512; 
    hop_size    = 128;
    window      = hann(window_size, 'periodic');
    
    [S, f, t] = stft(audio, fs, 'Window', window, ...
        'OverlapLength', window_size-hop_size, 'FFTLength', window_size);
    
    magnitude_spectrum  = abs(S);
    Spectral_flux = sum(diff(magnitude_spectrum, 1, 2).^2, 1);
    Spectral_flux = [Spectral_flux(1), Spectral_flux];            
    
    fluxNorm = (Spectral_flux - min(Spectral_flux)) / (max(Spectral_flux) - min(Spectral_flux));
    
    mu    = mean(fluxNorm);
    sigma = std(fluxNorm);
    kThr  = 1.0;                                   
    threshold = mu + kThr * sigma;
    
    useful_frames = find(fluxNorm > threshold);
    used_frames_times  = t(useful_frames);
    
    fprintf('\nFile: %s | # Hits detected: %d\n', fileList{fileIdx}, numel(used_frames_times));
    
    hit_duration = zeros(size(used_frames_times));
    num_frames = size(magnitude_spectrum,2);
    
    for k = 1:numel(useful_frames)
        idx = useful_frames(k);
        frame_of_energy = magnitude_spectrum(:, idx); 
        decay_threshold = 0.2 * max(frame_of_energy);
        
        last_idx = idx;
        for fj = idx+1:num_frames
            if max(magnitude_spectrum(:,fj)) < decay_threshold
                last_idx = fj;
                break;
            end
        end
        
        hit_duration(k) = t(last_idx) - t(idx);
    end
    
    
    parameter_array = zeros(numel(useful_frames), 3); 
    
    for k = 1:numel(useful_frames)
        fr  = magnitude_spectrum(:, useful_frames(k));
        cent = sum(f .* fr) / sum(fr);
        bw   = sqrt(sum(((f - cent).^2) .* fr) / sum(fr));
        [~, mxidx] = max(fr);
        
        parameter_array(k, :) = [cent, bw, f(mxidx)];
    end
    
    
    cluster_max = min([5, numel(useful_frames)]);
    
    if numel(useful_frames) > 1
        Y = pdist(parameter_array, 'euclidean');
        Z = linkage(Y, 'ward');
        clusters = cluster(Z, 'maxclust', cluster_max);
    else
        clusters = 1;
    end
    
    fprintf('  Hit #   Time (s)   Duration (s)   Instrument Group\n');
    for k = 1:numel(useful_frames)
        fprintf('  %3d    %7.3f    %7.3f      %d\n', ...
            k, used_frames_times(k), hit_duration(k), clusters(k));
    end
    fprintf('Estimated # of distinct instruments in this file: %d\n', ...
        numel(unique(clusters)));
    
    figure;
    uniq_clusters = unique(clusters);
    colors = lines(numel(uniq_clusters)); 
    
    hold on;
    for c = 1:numel(uniq_clusters)
        idx = clusters == uniq_clusters(c);
        plot(used_frames_times(idx), hit_duration(idx), ...
             'LineStyle', 'none', ...               
             'Marker', 'o', ...                     
             'MarkerFaceColor', colors(c,:), ...   
             'MarkerEdgeColor', colors(c,:), ...    
             'MarkerSize', 8);                      
    end
    hold off;
    
    xlabel('Hit Time (s)');
    ylabel('Hit Duration (s)');
    title(sprintf('Drum Hit Clustering - Time Domain (%s)', fileList{fileIdx}));
    grid on;
    legend(arrayfun(@(x) sprintf('Group %d', x), uniq_clusters, 'UniformOutput', false), ...
           'Location', 'best');
end