%% Program to generate an acoustic energy map of a WASN
close all
clear
clc
nframes = 1024;
c = 343;
fs = 48000;
ts = 1/fs;
%% Proof points, mics position & steering vector
x_mic1 = 0.32;        %input('x coordinate of mic1: ');
y_mic1 = 0.325;       %input('y coordinate of mic1: ');
n_mics = 4;
x_proof_limit = 0.5;
y_proof_limit = 0.5;
grid_size = 200;
x = -x_proof_limit:2*x_proof_limit/grid_size:x_proof_limit;
y = -y_proof_limit:2*x_proof_limit/grid_size:y_proof_limit;
n_proof_x = 35;
n_proof_y = 35;
n_proof = n_proof_x*n_proof_y;
[X_mics,X_p] = proof(x_mic1,y_mic1,x,y,n_mics,n_proof_x,n_proof_y);

%% Signal simulation 
n_sources = input('Number of sources: ');
ome = 100;                              % angular frequency
lambda = c/ome;
n_lambdas = 3;
sec = n_lambdas*lambda/c;                                  % number of seconds in the signal
t = 0:ts:sec-ts;                            % time vector (1 second)
N = length(t);
SS = zeros(n_sources,N);
X_S = zeros(n_sources,2);
for i = 1:n_sources
    x_ss = input(['X coordinate of the ',num2str(i),' source: ']);
    y_ss = input(['X coordinate of the ',num2str(i),' source: ']);
    X_S(i,1) = x_ss;
    X_S(i,2) = y_ss;
    omega = input(['Frequency of the ',num2str(i),' source: ']);
    s = cos(2*pi*omega*t); 
    SS(i,:) = s;
end
[signals] = signal_sim(n_mics,SS,X_S,fs,X_mics,N,n_sources,c);
figure(1)
plot(X_p(:,1),X_p(:,2),'k*')
hold on 
plot(X_mics(:,1),X_mics(:,2),'r*')
plot(X_S(:,1),X_S(:,2),'bo')
legend('Proof points','Mics position','Source position')

figure(2)
plot(t,real(ifft(signals(1,:))),t,real(ifft(signals(2,:))),t,real(ifft(signals(3,:))),t,real(ifft(signals(4,:))))
legend('mic1','mic2','mic3','mic4')
%% Acoustic map generation
E_MVDR = zeros(1,n_proof);
E_FREQ = zeros(1,n_proof);
amp_out = 10;                           %post-amplification for beamformer output (original: 1)
phase_diff_threshold = 5*pi/180;      %mask threshold in degrees (original:0.01)

for i = 1:n_proof
    w = [0 1:N/2+1 (-N/2+1):-1]/N*fs;
    [T_p,w_c] = steering_vector(X_p,X_mics,c,n_mics,i,N,w);
    o_f_MVDR = zeros(1,N);
    o_f_FREQ = zeros(1,N);
    o_f_FREQ(1) = signals(1,1);
    [E_MVDR(1,i),E_FREQ(1,i)] = mapping(signals,N,n_mics,i,w_c,phase_diff_threshold);
end

figure(1000)
subplot(1,2,1)
scatter3(X_p(:,1),X_p(:,2),E_MVDR,'*')
subtitle('MVDR beamformer','FontSize',12)
subplot(1,2,2)
scatter3(X_p(:,1),X_p(:,2),E_FREQ,'*')
subtitle('Frequency mask (beamformer)','FontSize',12)
sgtitle('Acoustic energy mapping')