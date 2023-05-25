function [E_MVDR,E_FREQ,E_DAS] = mapping(signals,N,n_mics,i,w_c,phase_diff_threshold)
%% Beamformers
for f = 2:N
    %% DELAY AND SUM
    o_f_DAS(f) = (w_c(:,f)'*signals(:,f))/n_mics;
    %% MVDR
    R = signals(:,f)*signals(:,f)';
    for m =1:n_mics
        R(m,m) = 1.1*R(m,m);
    end
    inv_R = pinv(R);
    w_a = w_c(:,f);
    w_o = (inv_R*w_a)/(w_a'*inv_R*w_a);
    o_f_MVDR(f) = w_o'*signals(:,f);
    %% FREQ_MASK
    aligned_f(1) = angle(signals(1,f));
    for m = 2:n_mics
        aligned_f(m) = angle(w_c(m,f)'*signals(m,f));
    end
    phase_diffs_sum = 0;
    phase_diffs_num = 0;
    for m1 = 1:n_mics-1
        for m2 = m1+1:n_mics
            this_diff = abs(aligned_f(m1)-aligned_f(m2));
            if this_diff > pi
                this_diff = 2*pi - this_diff;
            end
            phase_diffs_sum = phase_diffs_sum + this_diff;
            phase_diffs_num = phase_diffs_num + 1;
        end
    end
    phase_diff = phase_diffs_sum/phase_diffs_num;
    if(phase_diff < phase_diff_threshold)
        freq_mask(f) = 1;
    else
        freq_mask(f) = 0;
    end
    o_f_FREQ(f) =freq_mask(f)*signals(1,f);
end

o_MVDR= real(ifft(o_f_MVDR));
o_FREQ = real(ifft(o_f_FREQ));
o_DAS= real(ifft(o_f_DAS));
%E_MVDR = o_MVDR.^2;
%E_FREQ = o_FREQ.^2;
%E_DAS = o_DAS.^2;
E_MVDR = sqrt(norm(o_MVDR));
E_FREQ = sqrt(norm(o_FREQ));
E_DAS = sqrt(norm(o_DAS));
%E_MVDR = trapz(times(i,:),E_MVDR);
%E_FREQ = trapz(times(i,:),E_FREQ);
%E_DAS = trapz(times(i,:),E_DAS);
end