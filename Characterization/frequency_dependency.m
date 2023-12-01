%% Program to generate an acoustic energy map of a WASN
close all
clear
clc
nframes = 1024;         % window size
c = 343;                % speed of sound 
fs = 48000;             % sampling frequency
ts = 1/fs;              % sampling time 
%% Proof points & mics position
x_mic1 = 0.32;          % input('x coordinate of mic1: ');
y_mic1 = 0.32;          % input('y coordinate of mic1: ');
n_mics = 4;             % number of microphones
% proof points definition
x_proof_limit = 0.4;    % grid x-limit for proof points
y_proof_limit = 0.4;    % grid y-limit for proof points
grid_size = 200;        % grid size
x = -x_proof_limit:2*x_proof_limit/grid_size:x_proof_limit;
y = -y_proof_limit:2*x_proof_limit/grid_size:y_proof_limit;

%% Signal simulation 
n_sources = input('Number of sources: ');
X_S = zeros(n_sources,2);

for i = 1:n_sources
    x_ss = input(['X coordinate of the ',num2str(i),' source: ']);
    y_ss = input(['Y coordinate of the ',num2str(i),' source: ']);
    X_S(i,1) = x_ss;
    X_S(i,2) = y_ss;
    %omega = input(['Frequency of the ',num2str(i),' source: ']);
    %lambda = c/omega;
    %n_lambdas = 10;
    
end

samples = 1024;
defase = 0;
sec = samples*ts;             % number of seconds in the signal
t = 0:ts:sec-ts;                            % time vector (1 second)
N = length(t);

n_proof_x = 45;         % number of proof poinst in x axis
n_proof_y = 45;         % number of proof poinst in y axis
n_proof = n_proof_x*n_proof_y;  %total proof points in a squared grid
[X_mics,X_p,x_p,y_p] = proof(x_mic1,y_mic1,x,y,n_mics,n_proof_x,n_proof_y); % proof points and mics location 


count = 0;
ind_beg = 100;
step = 100;
ind_end = 4000;

time_das_vect = zeros((ind_end-ind_beg)/step+1,n_proof);
time_mvdr_vect = zeros((ind_end-ind_beg)/step+1,n_proof);
time_pbm_vect = zeros((ind_end-ind_beg)/step+1,n_proof);
time_spr_vect = zeros((ind_end-ind_beg)/step+1,n_proof);
SS = zeros(n_sources,N);
for omega = ind_beg:step:ind_end
    count = count+1;
    
    s = cos(2*pi*omega*t);
    SS(1,:) = s;
    [signals,times] = signal_sim(n_mics,SS,X_S,X_mics,N,n_sources,c,t);
    fig_input = figure('Name',num2str(100+1),'visible','off');
    %plot(times(1,:),real(ifft(signals(1,:))),'+',times(2,:),real(ifft(signals(2,:))),'*',times(3,:),real(ifft(signals(3,:))),'o',times(4,:),real(ifft(signals(4,:))),'-')
    plot(t,real(ifft(signals(1,:))),'+',t,real(ifft(signals(2,:))),'*',t,real(ifft(signals(3,:))),'o',t,real(ifft(signals(4,:))),'-')
    xlim([0 t(1,N)])
    legend({'Node 1','Node 2','Node 3','Node 4'},'Location','southwest')
    legend('boxoff')
    xlabel('Time [s]')
    ylabel('Amplitude')
    title('Input signals: frequency ',[num2str(omega),' Hz.'])
    saveas(fig_input,['input',num2str(omega),'.png'])   
%end
    %% Acoustic map generation
    E_MVDR = zeros(1,n_proof);
    E_FREQ = zeros(1,n_proof);
    E_SPR = zeros(1,n_proof);
    E_DAS = zeros(1,n_proof);
    Pow = zeros(1,n_proof);
    amp_out = 10;                           %post-amplification for beamformer output (original: 1)
    phase_diff_threshold = 3*pi/180;      %mask threshold in degrees (original:0.01)
    w = [0 1:N/2+1 (-N/2+1):-1]/N*fs;
    for i = 1:n_proof 
        [T_p,w_c] = steering_vector(X_p,X_mics,c,n_mics,i,N,w);
        o_f_DAS = zeros(1,N);
        o_f_MVDR = zeros(1,N);
        o_f_FREQ = zeros(1,N);
        o_f_FREQ(1) = signals(1,1);
        Pow_n = 0;
        [E_MVDR(1,i),E_FREQ(1,i),E_DAS(1,i),o_DAS,o_FREQ,o_MVDR,time_das,time_mvdr,time_pbm] = mapping(signals,N,n_mics,i,w_c,phase_diff_threshold,t);
        time_das_vect(count,i) = time_das;
        time_mvdr_vect(count,i) = time_mvdr;
        time_pbm_vect(count,i) = time_pbm;
        tic
        for h = 1:n_mics
            for j = h+1:n_mics
                q_rm = sqrt((X_p(i,1)-X_mics(h,1))^2+(X_p(i,2)-X_mics(h,2))^2);
                q_rn = sqrt((X_p(i,1)-X_mics(j,1))^2+(X_p(i,2)-X_mics(j,2))^2);
                tau = (q_rm-q_rn)/c;
    
                R_tmp = signals(h,:).*conj(signals(j,:));
                R_tmpabs=abs(R_tmp);
                R_tmpabs(R_tmpabs == 0) = eps;
                R_nm = (R_tmp./R_tmpabs).*exp(1i*w(1:end-1)*tau);
                Pow_n = Pow_n + sum(R_nm);
            end
        end
        Pow(i) = Pow_n;
        [val_srp(count),val_srp_ind(count)]  = max(Pow);
        E_SPR(i) = norm(Pow(i));
        time_spr_vect(count,i) = toc;
    end
    
    %% Plotting energy maps
    fig_maps = figure('Name',num2str(200+count),'visible','off');
    %figure(200+count)
    subplot(2,2,1)
    scatter3(X_p(:,1),X_p(:,2),E_MVDR,'*')
    subtitle('MVDR beamformer','FontSize',12)
    subplot(2,2,2)
    scatter3(X_p(:,1),X_p(:,2),E_DAS,'*')
    subtitle('DAS beamformer','FontSize',12)
    sgtitle('Acoustic energy mapping')
    subplot(2,2,3)
    scatter3(X_p(:,1),X_p(:,2),E_SPR,'*')
    subtitle('SPR beamformer','FontSize',12)
    subplot(2,2,4)
    scatter3(X_p(:,1),X_p(:,2),E_FREQ,'*')
    subtitle('PBM beamformer','FontSize',12)
    sgtitle(['Acoustic energy mapping for a sample size of ',[num2str(n_proof),' samples']])
    saveas(fig_maps,['maps',num2str(count),'.png'])

    E_DAS_slides = zeros(n_proof_x,n_proof_y);
    E_MVDR_slides = zeros(n_proof_x,n_proof_y);
    E_SPR_slides = zeros(n_proof_x,n_proof_y);
    E_FREQ_slides = zeros(n_proof_x,n_proof_y);
    for i = 1:n_proof_x
        E_DAS_slides(i,:) = E_DAS((i-1)*n_proof_x+1:i*n_proof_x);
        E_MVDR_slides(i,:) = E_MVDR((i-1)*n_proof_x+1:i*n_proof_x);
        E_SPR_slides(i,:) = E_SPR((i-1)*n_proof_x+1:i*n_proof_x);
        E_FREQ_slides(i,:) = E_FREQ((i-1)*n_proof_x+1:i*n_proof_x);
    end
    x_label = string(round(-x_proof_limit:2*x_proof_limit/(n_proof_x-1):x_proof_limit,2));
    y_label = string(round(-x_proof_limit:2*x_proof_limit/(n_proof_x-1):x_proof_limit,2));
    
    fig_das = figure('Name',num2str(300+count),'visible','off');
    %figure(300+count)
    set(gcf,'Position',[200 200 700 600])
    h = heatmap(x_label,y_label,E_DAS_slides,'Colormap',parula);
    h.XLabel = 'X coordinate'; h.YLabel = 'Y coordinate';     
    h.Title = ['Heatmap DAS for a frequency of ',[num2str(omega),' Hz']];
    saveas(fig_das,['das',num2str(count),'.png'])
    fig_spr = figure('Name',num2str(400+count),'visible','off');
    %figure(400+count)
    set(gcf,'Position',[200 200 700 600])
    k = heatmap(x_label,y_label,E_SPR_slides,'Colormap',parula);
    k.XLabel = 'X coordinate'; k.YLabel = 'Y coordinate'; 
    k.Title = ['Heatmap SRP for a frequency of ',[num2str(omega),' Hz']];
    saveas(fig_spr,['spr',num2str(count),'.png'])
    fig_pbm = figure('Name',num2str(500+count),'visible','off');
    %figure(500+count)
    set(gcf,'Position',[200 200 700 600])
    k = heatmap(x_label,y_label,E_FREQ_slides,'Colormap',parula);
    k.XLabel = 'X coordinate'; k.YLabel = 'Y coordinate'; 
    k.Title = ['Heatmap PBM for a frequency of ',[num2str(omega),' Hz']];
    saveas(fig_pbm,['pbm',num2str(count),'.png'])
    fig_mvdr = figure('Name',num2str(600+count),'visible','off');
    %figure(600+count)
    set(gcf,'Position',[200 200 700 600])
    k = heatmap(x_label,y_label,E_MVDR_slides,'Colormap',parula);
    k.XLabel = 'X coordinate'; k.YLabel = 'Y coordinate'; 
    k.Title = ['Heatmap MVDR for a frequency of ',[num2str(omega),' Hz']];
    saveas(fig_mvdr,['mvdr',num2str(count),'.png'])
end