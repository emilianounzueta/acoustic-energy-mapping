function [time,full_data] = signal_adquisicion_net(folder,recording_folder)
    %recording = recording_folder;
    %folder = folders(1,:);
    %% Opening audio, timestamp & synchro files
    nframes = 1024; %this is assumed
    fs = 48000;
    %for mic1
    %recording = 'kernel8/netjack'; folder = 'prueba01';
    
    path_sec1 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/sec.txt';
    path_sec1 = strrep(path_sec1,'folder',folder);
    path_sec1 = strrep(path_sec1,'recording',recording_folder);
    timefid = fopen(path_sec1);
    timedata1 = fread(timefid, Inf, "uint64");
    fclose(timefid);
    
    % for mic1
    path_wav1 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/rasp1.wav';
    path_wav1 = strrep(path_wav1,'folder',folder);
    path_wav1 = strrep(path_wav1,'recording',recording_folder);
    [data1, ~] = audioread(path_wav1);
    % for mic2
    path_wav2 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/rasp2.wav';
    path_wav2 = strrep(path_wav2,'folder',folder);
    path_wav2 = strrep(path_wav2,'recording',recording_folder);
    [data2, ~] = audioread(path_wav2);
    % for mic3
    path_wav3 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/rasp3.wav';
    path_wav3 = strrep(path_wav3,'folder',folder);
    path_wav3 = strrep(path_wav3,'recording',recording_folder);
    [data3, ~] = audioread(path_wav3);
    % for mic4 
    path_wav4 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/rasp4.wav';
    path_wav4 = strrep(path_wav4,'folder',folder);
    path_wav4 = strrep(path_wav4,'recording',recording_folder);
    [data4, ~] = audioread(path_wav4);

    %% Processing the timestamp vectors
    timedata = timedata1(1) + timedata1(2);
    timedata = timedata(1:end)/1000000; % so that timedata is in seconds
    

    %% Creating time vectors
    periodo = 1/fs;
    for i = 1:length(data1)
        timedata(i) =  timedata(1) + periodo*i;
    end

    %% Results
    time(1:length(timedata),1) = timedata;
    full_data = zeros(max([length(data1) length(data2) length(data3) length(data4)]),4);
    full_data(1:length(data1),1) = data1;
    full_data(1:length(data2),2) = data2;
    full_data(1:length(data3),3) = data3;
    full_data(1:length(data4),4) = data4;
end 