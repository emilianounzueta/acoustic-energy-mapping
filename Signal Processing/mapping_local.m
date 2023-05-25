%% Program to generate an acoustic energy map of a WASN
close all
clear
clc
nframes = 1024;
c = 343;
fs = 48000;
ts = 1/fs;
%% Proof points, mics position & steering vector
x_mic1 = 0.32;    
y_mic1 = 0.325;   
n_mics = 4;
x_proof_limit = 0.4;
y_proof_limit = 0.4;
grid_size = 200;
x = -x_proof_limit:2*x_proof_limit/grid_size:x_proof_limit;
y = -y_proof_limit:2*x_proof_limit/grid_size:y_proof_limit;
n_proof_x = 15;
n_proof_y = 15;
n_proof = n_proof_x*n_proof_y;
[X_mics,X_p,x_p,y_p] = proof(x_mic1,y_mic1,x,y,n_mics,n_proof_x,n_proof_y);
%%
t_synchrorecording_folder = input('Name of recording folder: ','s');
n = input('Number of recordings: ');
folders = char(zeros(n, 8));

for i = 1+0:n+0
    if i < 10
        folders(i-0,:) = ['prueba0', num2str(i)];
    else
        folders(i-0,:) = ['prueba', num2str(i)];
    end
end 

full_times = cell(n,1);
full_datas = cell(n,1);
begin_end = cell(n,1);
for i = 1:n
    [full_times{i},full_datas{i},begin_end{i}] = signal_adquisicion(folders(i,:),recording_folder);
end 

mics = [1 2 3 4];
sample_size = zeros(n,4);
beg_end = zeros(n,2);
begin = zeros(n,4);

%% Getting matrix for all signals
for i = 1:n
     beg_end = cell2mat(begin_end(i,1));
     begin(i,:) = beg_end(:,1);
end
%% Statistic analizis
% time difference in recording
for i = 1:n
    for j = 1:4
        time_diff_beg(i,j) = begin(i,j) - begin(i,4); %ref = rasp5 
    end
end

samples_diff_beg = floor(time_diff_beg*48000);
[max_diff,index_max] = max(samples_diff_beg,[],2);
epoch_start_synchro = zeros(n,4);
for i = 1:n
    for j = 1:4
        if j ~= index_max
            epoch_start_synchro(i,j) = time_diff_beg(i,index_max(i))-time_diff_beg(i,j);
        else
            epoch_start_synchro(i,j) = 0;
        end
    end
end
epoch_samples = floor(epoch_start_synchro*48000);

%% Signal adquisition & processing
E_MVDR = zeros(n,n_proof);
E_DAS = zeros(n,n_proof);
E_FREQ = zeros(n,n_proof);
E_SPR = zeros(n,n_proof);

synchro_error = zeros(n,4);
%wait = waitbar(0,'Barra de progreso');
sync_index = zeros(n,4);
syncro_times = zeros(n,4);
count = 0;
for i = 1:n
    %waitbar(i/n);
    count = count+1;
    times_full = cell2mat(full_times(i,1));
    datas_full = cell2mat(full_datas(i,1));
    [signals,times,sync_index] = processing(times_full,datas_full,i,begin,sync_index);

    N = length(signals);
    w = [0 1:N/2+1 (-N/2+1):-1]/N*fs;
    syncro_times(i,:) = [times(2,sync_index(i,1)) times(3,sync_index(i,2)) times(1,sync_index(i,3)) times(4,sync_index(i,4))]
    %% Acoustic map generation
    amp_out = 10;                           %post-amplification for beamformer output (original: 1)
    phase_diff_threshold = 5*pi/180;      %mask threshold in degrees (original:0.01)

    figure(5000+i)
    plot(real(ifft(signals(1,:))),'k')
    hold on
    plot(real(ifft(signals(2,:))),'g')
    plot(real(ifft(signals(3,:))),'b')
    plot(real(ifft(signals(4,:))),'r')
    ylabel('Amplitud')
    xlabel('Time [s]')
    title('Input signals')
    legend('data1','data3','data2','data4')  
end
    Pow = zeros(1,n_proof);
    Pow_n = 0;
    for j = 1:n_proof 
        [T_p,w_c] = steering_vector(X_p,X_mics,c,n_mics,j,N,w);
        o_f_MVDR = zeros(1,N);
        o_f_FREQ = zeros(1,N);
        o_f_FREQ(1) = signals(1,1);
        Pow_n = 0;
        [E_MVDR(i,j),E_FREQ(i,j),E_DAS(i,j)] = mapping(signals,N,n_mics,i,w_c,phase_diff_threshold);
        for h = 1:n_mics
            for jj = h+1:n_mics
                q_rm = sqrt((X_p(j,1)-X_mics(h,1))^2+(X_p(j,2)-X_mics(h,2))^2);
                q_rn = sqrt((X_p(j,1)-X_mics(jj,1))^2+(X_p(j,2)-X_mics(jj,2))^2);
                tau = (q_rm-q_rn)/c;
                R_tmp = signals(h,:).*conj(signals(jj,:));
                R_tmpabs=abs(R_tmp);
                R_tmpabs(R_tmpabs == 0) = eps;
                R_nm = (R_tmp./R_tmpabs).*exp(1i*w(1:end-1)*tau);
                Pow_n = Pow_n + sum(R_nm);
            end
        end
        Pow(j) = Pow_n;
        %[val_srp(count),val_srp_ind(count)]  = max(Pow);
        E_SPR(i,j) = norm(Pow(j));
    end
%     figure(6000+20+i)
%     subplot(1,2,1)
%     scatter3(X_p(:,1),X_p(:,2),E_MVDR(i,:),'*')
%     subtitle('MVDR beamformer','FontSize',12)
%     subplot(1,2,2)
%     scatter3(X_p(:,1),X_p(:,2),E_DAS(i,:),'*')
%     subtitle('Frequency mask (beamformer)','FontSize',12)
%     sgtitle(['Acoustic energy mapping for the ',num2str(i+0),' recording'])
end 
%close(wait);
for kk = 1:n
    count = kk;
    E_DAS_slides = zeros(n_proof_x,n_proof_y);
    E_MVDR_slides = zeros(n_proof_x,n_proof_y);
    E_SPR_slides = zeros(n_proof_x,n_proof_y);
    E_FREQ_slides = zeros(n_proof_x,n_proof_y);
    for mm = 1:n_proof_x
        E_DAS_slides(mm,:) = E_DAS(kk,(mm-1)*n_proof_x+1:(mm)*n_proof_x);
        E_MVDR_slides(mm,:) = E_MVDR(kk,(mm-1)*n_proof_x+1:(mm)*n_proof_x);
        E_SPR_slides(mm,:) = E_SPR(kk,(mm-1)*n_proof_x+1:(mm)*n_proof_x);
        E_FREQ_slides(mm,:) = E_FREQ(kk,(mm-1)*n_proof_x+1:(mm)*n_proof_x);
    end
    x_label = string(round(-x_proof_limit:2*x_proof_limit/(n_proof_x-1):x_proof_limit,2));
    y_label = string(round(-x_proof_limit:2*x_proof_limit/(n_proof_x-1):x_proof_limit,2));
    %cdl = k.XDisplayLabels;
    %labels = [sprintf('%0.2f',-x_proof_limit); repmat('     ',(size(cdl,1)-1)/2-1, (size(cdl,2)));'0.00 ' ;repmat('     ',(size(cdl,1)-1)/2-1, (size(cdl,2))); sprintf('%0.2f ',x_proof_limit)];
%     figure(300+kk)
%     set(gcf,'Position',[200 200 700 600])
%     k = heatmap(x_label,y_label,E_DAS_slides,'Colormap',parula,'FontSize',14);
%     s = struct(k);
%     s.YAxis.TickLabelRotation = 90;
%     s.XAxis.TickLabelRotation = 0;
%     k.XLabel = 'X coordinate [m]'; k.YLabel = 'Y coordinate [m]';
%     k.Title = 'Heatmap for DAS bemaformer';
%     cdl = k.XDisplayLabels;
%     k.XDisplayLabels = [sprintf('%0.2f',x_proof_limit); repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2)));'0.00' ;repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2))); sprintf('%0.2f',x_proof_limit)];
%     k.YDisplayLabels = [sprintf('%0.2f',x_proof_limit); repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2)));'0.00' ;repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2))); sprintf('%0.2f',x_proof_limit)];
%     k.ColorbarVisible = 'off';
%     figure(400+kk)
%     set(gcf,'Position',[200 200 700 600])
%     k = heatmap(x_label,y_label,E_SPR_slides,'Colormap',parula,'FontSize',14);
%     s = struct(k);
%     s.YAxis.TickLabelRotation = 90;
%     s.XAxis.TickLabelRotation = 0;
%     k.XLabel = 'X coordinate [m]'; k.YLabel = 'Y coordinate [m]';
%     k.Title = 'Heatmap for SRP bemaformer';
%     cdl = k.XDisplayLabels;
%     k.XDisplayLabels = [sprintf('%0.2f',x_proof_limit); repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2)));'0.00' ;repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2))); sprintf('%0.2f',x_proof_limit)];
%     k.YDisplayLabels = [sprintf('%0.2f',x_proof_limit); repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2)));'0.00' ;repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2))); sprintf('%0.2f',x_proof_limit)];
%     k.ColorbarVisible = 'off';
%     figure(500+kk)
%     set(gcf,'Position',[200 200 700 600])
%     k = heatmap(x_label,y_label,E_FREQ_slides,'Colormap',parula,'FontSize',14);
%     s = struct(k);
%     s.YAxis.TickLabelRotation = 90;
%     s.XAxis.TickLabelRotation = 0;
%     k.XLabel = 'X coordinate [m]'; k.YLabel = 'Y coordinate [m]';
%     k.Title = 'Heatmap for PBM bemaformer';
%     cdl = k.XDisplayLabels;                                    % Current Display Labels
%     k.XDisplayLabels = [sprintf('%0.2f',x_proof_limit); repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2)));'0.00' ;repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2))); sprintf('%0.2f',x_proof_limit)];
%     k.YDisplayLabels = [sprintf('%0.2f',x_proof_limit); repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2)));'0.00' ;repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2))); sprintf('%0.2f',x_proof_limit)];
%     k.ColorbarVisible = 'off';
%     figure(600+kk)
%     set(gcf,'Position',[200 200 700 600])
%     k = heatmap(x_label,y_label,E_MVDR_slides,'Colormap',parula,'FontSize',14);
%     s = struct(k);
%     s.YAxis.TickLabelRotation = 90;
%     s.XAxis.TickLabelRotation = 0;
%     k.XLabel = 'X coordinate [m]'; k.YLabel = 'Y coordinate [m]';
%     k.Title = 'Heatmap for MVDR bemaformer';
%     cdl = k.XDisplayLabels;                                    % Current Display Labels
%     k.XDisplayLabels = [sprintf('%0.2f',x_proof_limit); repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2)));'0.00' ;repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2))); sprintf('%0.2f',x_proof_limit)];
%     k.YDisplayLabels = [sprintf('%0.2f',x_proof_limit); repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2)));'0.00' ;repmat('    ',(size(cdl,1)-1)/2-1, (size(cdl,2))); sprintf('%0.2f',x_proof_limit)];
%     k.ColorbarVisible = 'off';
    for i = 1:n_proof_x
        Max_das(1,i) = max(E_DAS_slides(i,:));
        Max_das(2,i) = max(E_DAS_slides(:,i));
        Max_mvdr(1,i) = max(E_MVDR_slides(i,:));
        Max_mvdr(2,i) = max(E_MVDR_slides(:,i));
        Max_spr(1,i) = max(E_SPR_slides(i,:));
        Max_spr(2,i) = max(E_SPR_slides(:,i));
        Max_freq(1,i) = max(E_FREQ_slides(i,:));
        Max_freq(2,i) = max(E_FREQ_slides(:,i));
    end
    [max_das_x(count),ind_max_dasx] = max(Max_das(1,:),[],'all');
    [max_das_Y(count),ind_max_dasy] = max(Max_das(2,:),[],'all');
    [max_mvdr_x(count),ind_max_mvdrx] = max(Max_mvdr(1,:),[],'all');
    [max_mvdr_Y(count),ind_max_mvdry] = max(Max_mvdr(2,:),[],'all');
    [max_spr_x(count),ind_max_sprx] = max(Max_spr(1,:),[],'all');
    [max_spr_Y(count),ind_max_spry] = max(Max_spr(2,:),[],'all');
    [max_freq_x(count),ind_max_freqx] = max(Max_freq(1,:),[],'all');
    [max_freq_Y(count),ind_max_freqy] = max(Max_freq(2,:),[],'all');
    das_x(count)= x_p(ind_max_dasx); das_y(count) = y_p(ind_max_dasy);
    mvdr_x(count)= x_p(ind_max_mvdrx); mvdr_y(count) = y_p(ind_max_mvdry);
    spr_x(count)= x_p(ind_max_sprx); spr_y(count) = y_p(ind_max_spry);
    freq_x(count)= x_p(ind_max_freqx); freq_y(count) = y_p(ind_max_freqx);
    
%     figure(700+count)
%     set(gcf,'Position',[0 0 700 650])
%     plot(X_mics(:,1),X_mics(:,2),'r*','MarkerSize',15)
%     hold on 
%     plot(0,0,'bo','MarkerSize',15)
%     plot(das_x(count),-das_y(count),'^','color',[0.8500 0.3250 0.0980],'MarkerSize',20,'LineWidth',2)
%     plot(mvdr_x(count),-mvdr_y(count),'<','color',[0.4660 0.6740 0.1880], 'MarkerSize',20,'LineWidth',2)
%     plot(spr_x(count),-spr_y(count),'>k', 'MarkerSize',20,'LineWidth',2)
%     plot(freq_x(count),-freq_y(count),'v','color',[0.4940 0.1840 0.5560], 'MarkerSize',20,'LineWidth',2)
%     legend({'Mics position','Source position','Calculated position (DAS)','Calculated position (MVDR)','Calculated position (SPR)','Calculated position (PBM)'},'Location','east')
%     title('Calculated position for the',[num2str(count),'recording'])
%     xlabel('x coordinate [m]')
%     ylabel('y coordinate [m]')
    err_DAS(count) = sqrt((0-x_p(ind_max_dasx))^2 + (0-y_p(ind_max_dasy))^2);
    err_MVDR(count) = sqrt((0-x_p(ind_max_mvdrx))^2 + (0-y_p(ind_max_mvdry))^2);
    err_SPR(count) = sqrt((0-x_p(ind_max_sprx))^2 + (0-y_p(ind_max_spry))^2);
    err_FREQ(count) = sqrt((0-x_p(ind_max_freqx))^2 + (0-y_p(ind_max_freqx))^2);
    errores(count,:) = [err_DAS(count), err_MVDR(count), err_SPR(count), err_FREQ(count)]';
end

% Statistical analysis of synchronization 
t_synchro = zeros(n,4);
for syn = 1:n
    for synn = 1:4
        if synn ~= 3
            t_synchro(syn,synn) = syncro_times(syn,3)-syncro_times(syn,synn);
        else
            t_synchro(syn,synn) = 0;
        end
    end
end
%% Find maximum value from a threshold
figure(900)
set(gcf,'Position',[750 0 700 500])
grid on
plot(errores(:,1),'^','color',[0 0.4470 0.7410], 'MarkerSize',20,'LineWidth',2)
hold on 
plot(errores(:,2)','<','color',[0.4660 0.6740 0.1880], 'MarkerSize',20,'LineWidth',2)
plot(errores(:,3)','>k', 'MarkerSize',20,'LineWidth',2)
plot(errores(:,4)','vr', 'MarkerSize',20,'LineWidth',2)
legend('DAS','MVDR','SPR','PBM')
title('Localization error vs synchronization error in node 1')
xlabel('Error in synchronization in node 1 [s]')
ylabel('Error in localization of source [m]')

error_mean = [mean(errores(:,1)) mean(errores(:,2)) mean(errores(:,3)) mean(errores(:,4))];
error_mad = [mad(errores(:,1)) mad(errores(:,2)) mad(errores(:,3)) mad(errores(:,4))];
error_std= [std(errores(:,1)) std(errores(:,2)) std(errores(:,3)) std(errores(:,4))];