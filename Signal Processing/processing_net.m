%% Adquiring signals
function [signals,times,sync_index] = processing_net(times_full,datas_full,i,sync_index)
    % Folder 
    [t] = deal(times_full(:,1));
    [data1,data2,data3,data4] = deal(datas_full(:,1),datas_full(:,2),datas_full(:,3),datas_full(:,4));
    %% processing the signals
    index1 = find(data1~=0,1,'first');
    index2 = find(data2~=0,1,'first');
    index3 = find(data3~=0,1,'first');
    index4 = find(data4~=0,1,'first');
    data1(data1==0)=mean(data1);
    data2(data2==0)=mean(data2);
    data3(data3==0)=mean(data3);
    data4(data4==0)=mean(data4);

    data1_c = data1 - mean(data1);  
    data2_c = data2 - mean(data2);
    data3_c = data3 - mean(data3);
    data4_c = data4 - mean(data4);

    % plot full input signals
%     figure(1000+i)
%     plot(t,data1_c)
%     hold on
%     plot(t,data2_c)
%     plot(t,data3_c)
%     plot(t,data4_c)
%     ylabel('Amplitud')
%     xlabel('Time [s]')
%     title('Input signals for the',[num2str(i),' recording'])
%     legend('data1','data2','data3','data4')
    % Encontrar un máximo general y normalizar con respecto a ese!!!!
    max_data = max([max(abs(data1_c)) max(abs(data2_c)) max(abs(data3_c)) max(abs(data4_c))]);
    data1_c = data1_c/max_data;  
    data2_c = data2_c/max_data;
    data3_c = data3_c/max_data;
    data4_c = data4_c/max_data;
    %error = max([index1 index2 index3 index4])-min([index1 index2 index3 index4]);
    %data2_c(error+1:end) = data2_c(1:end-error);
    %% 
%     figure(3000+i)
%     plot(t,data1_c,'k')
%     hold on
%     plot(t,data2_c,'b')
%     plot(t,data3_c,'g')
%     plot(t,data4_c,'r')
%     ylabel('Amplitud')
%     xlabel('Time [s]')
%     title('Adjusted input signals for the',[num2str(i),' recording'])
%     legend('data1','data2','data3','data4')
    
    
    %% RMS for the windows
    min_ind = length(t);
    rms_s = zeros(floor(min_ind/1024*2),4);
    for rms_ind = 1:min_ind/1024
        rms_s(rms_ind,1) = rms(data1_c((rms_ind-1)*1024+1:rms_ind*1024));
        rms_s(rms_ind,2) = rms(data2_c((rms_ind-1)*1024+1:rms_ind*1024));
        rms_s(rms_ind,3) = rms(data3_c((rms_ind-1)*1024+1:rms_ind*1024));
        rms_s(rms_ind,4) = rms(data4_c((rms_ind-1)*1024+1:rms_ind*1024));
    end
    max_rms = [max(rms_s(:,1)) max(rms_s(:,2)) max(rms_s(:,3)) max(rms_s(:,4))];
    for k = 1:4
        indexes_data(k) = find(rms_s(:,k)>max_rms(k)/2,1,'first');
    end
    minn = min(indexes_data(indexes_data>10));
    %% Crop data & time
    data1_c = data1_c((minn-1)*1024:(minn+3)*1024-1);
    data2_c = data2_c((minn-1)*1024:(minn+3)*1024-1);
    data3_c = data3_c((minn-1)*1024:(minn+3)*1024-1);
    data4_c = data4_c((minn-1)*1024:(minn+3)*1024-1);
    t = t((minn-1)*1024:(minn+3)*1024-1);
%     t1 = t1((minn-1)*1024:(minn+3)*1024-1);
%     t2 = t2((minn-1)*1024:(minn+3)*1024-1);
%     t3 = t3((minn-1)*1024:(minn+3)*1024-1);
%     t4 = t4((minn-1)*1024:(minn+3)*1024-1);
    try
        sync_index(i,1) = find(data2_c>max(data2_c)/5,1,'first');
        [pks,locs] = findpeaks(data2_c(sync_index(i,1):end));
        sync_index(i,1) = sync_index(i,1)+locs(1)-1;
    catch
        sync_index(i,1) = length(data2_c);
    end
    try
        sync_index(i,2) = find(data3_c>max(data3_c)/5,1,'first');
        [pks,locs] = findpeaks(data3_c(sync_index(i,2):end));
        sync_index(i,2) = sync_index(i,2)+locs(1)-1;
    catch
        sync_index(i,2) = length(data3_c);
    end
    try
        sync_index(i,3) = find(data1_c>max(data1_c)/5,1,'first');
        [pks,locs] = findpeaks(data1_c(sync_index(i,3):end));
        sync_index(i,3) = sync_index(i,3)+locs(1)-1;
    catch
        sync_index(i,3) = length(data1_c);
    end
    try
        sync_index(i,4) = find(data4_c>max(data4_c)/5,1,'first');
        [pks,locs] = findpeaks(data4_c(sync_index(i,4):end));
        sync_index(i,4) = sync_index(i,4)+locs(1)-1;
    catch
        sync_index(i,4) = length(data4_c);
    end
    
    %% Post procesamiento
    %[data1_c,data2_c,data3_c,data4_c,t1,t2,t3,t4] = post_processing(data1_c,data2_c,data3_c,data4_c,t1,t2,t3,t4);
    %% Results
%     figure(4000+i)
%     plot(t,data1_c,'k')
%     hold on
%     plot(t,data2_c,'b')
%     plot(t,data3_c,'g')
%     plot(t,data4_c,'r')
%     ylabel('Amplitud')
%     xlabel('Time [s]')
%     title('Adjusted input signals for the',[num2str(i),' recording'])
%     legend('data1','data2','data3','data4')
    signals = [data1_c  data4_c data2_c data3_c]';
    times = [t t t t]'; 

    for m = 1:4
        signals(m,:) = fft(signals(m,:).*(hann(length(data1_c))'));
    end
    