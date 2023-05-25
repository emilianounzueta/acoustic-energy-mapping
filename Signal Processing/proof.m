function [X_mics,X_p,x_p,y_p] = proof(x_mic1,y_mic1,x,y,n_mics,n_proof_x,n_proof_y)
    X_mics = zeros(n_mics,2);
    for i = 1:n_mics
        if i <= n_mics/2
            X_mics(i,1) = -x_mic1;
        else
            X_mics(i,1) = x_mic1;
        end
        X_mics(i,2) = y_mic1*(-1)^i;
    end
    n_rows_p = n_proof_x-1;                     % n_rows = n_mics_x + 1
    n_columns_p= n_proof_y-1;                   % n_columns = n_mics_y + 1    
    delta_x_p = (x(end)-x(1))/n_rows_p;         % separation of proof points in x
    delta_y_p = (y(end)-y(1))/n_columns_p;      % separation of proof points in y
    x_p = zeros(1,n_proof_x);
    y_p = zeros(1,n_proof_y);  
    for i = 1:n_proof_x
        for j = 1:n_proof_y
            x_p(i) = x(1)+((i-1)*delta_x_p);
            y_p(j) = y(1)+((j-1)*delta_y_p);
        end
        X_p((i-1)*n_proof_y+1:i*n_proof_y,1) = x_p(i);
        for j = 1:n_proof_y
            X_p((i-1)*n_proof_y+j,2) = y_p(j);
        end
    end
end 
