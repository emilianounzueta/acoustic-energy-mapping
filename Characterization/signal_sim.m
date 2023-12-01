function [signals,times] = signal_sim(n_mics,SS,X_S,X_mics,N,n_sources,c,t)
fs = 48000; 
%defase = 0;
w = [0 1:N/2+1 (-N/2+1):-1]/N*fs;
signals = zeros(n_mics,length(SS(1,:)));
y = zeros(n_sources,N);
times = zeros(n_mics,N);
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
    signals(k,:) = [signals(k,ceil(time*fs):end) signals(k,1:floor(time*fs))]; 
    signals(k,:) = signals(k,:).*hann(length(signals(k,:)))';
    times(k,:) = [t(ceil(time*fs):end) t(1:floor(time*fs))]; 
end
% signal_temp = zeros(1,length(signals(1,:)));
% signal_temp(1:defase) = signals(1,end-defase+1:end);
% signal_temp(defase+1:end) = signals(1,1:end-defase);
% signals(1,:) = [zeros(1,defase) signal_temp(defase+1:end)];

%Fourier Transform
for m=1:n_mics
    signals(m,:) = fft(signals(m,:));
end
end
