function [T_p,w_c] = steering_vector(X_p,X_mics,c,n_mics,i,N,w)
    w_c = zeros(n_mics,N);
    w_c(1,:) = ones(1,N);
    t_ref = sqrt((X_p(i,1) - X_mics(1,1))^2 + (X_p(i,2) - X_mics(1,2))^2)/c;
    for k = 2:n_mics
        T_p(i,k) = sqrt((X_p(i,1) - X_mics(k,1))^2 + (X_p(i,2) - X_mics(k,2))^2)/c;
        T_p(i,k) = T_p(i,k) - t_ref;
    end
    for m = 1:n_mics-1
        for f = 1:N
            w_c(m+1,f)=exp(-1i*2*pi*w(f)*T_p(i,m+1));    % steering vector for this frequency
        end
    end
    
end