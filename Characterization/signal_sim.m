function [signals] = signal_sim(n_mics,SS,X_S,fs,X_mics,N,n_sources,c,defase)
    
    signals = zeros(n_mics,length(SS(1,:)));
    %time delay source-microphone
    w = [0 1:N/2 (-N/2+1):-1]/N*fs;                 %frequencies to be analized
    
    y = zeros(n_sources,N);
    for k = 1:n_mics
        for ns = 1:n_sources
            x_f = fft(SS(ns,:));                                 %Fourier Transform
            y_f = zeros(1,length(x_f));                     %output signal
            time = sqrt((X_S(ns,1) - X_mics(k,1)).^2 + (X_S(ns,2) - X_mics(k,2)).^2)/c;
            for f = 1:N
                e=exp(-1i*2*pi*w(f)*time);    % steering vector for this frequency
                y_f(f) = x_f(f)*e;
            end
            y(ns,:) = real(ifft(y_f));
            signals(k,:) = signals(k,:) + y(ns,:);
            
        end
    % Synchro analysis
    end
%     for i = 1:n_mics
%     [pks,locs] = findpeaks(signals(i,:));
%     ind_max(i) = locs(1);
%      %sync_index(i,1) = sync_index(i,1)+locs(1)-1;
%     end
    %Fourier Transform
    for m=1:n_mics
        signals(m,:) = fft(signals(m,:));
    end
end
