function [full_times,full_data,begin_end] = signal_adquisicion(folder,recording_folder)
    %recording = recording_folder;
    %folder = folders(1,:);
    %% Opening audio, timestamp & synchro files
    nframes = 1024; %this is assumed
    %for mic1
    path_sec1 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/sec1.txt';
    path_sec1 = strrep(path_sec1,'folder',folder);
    path_sec1 = strrep(path_sec1,'recording',recording_folder);
    timefid = fopen(path_sec1);
    timedata1 = fread(timefid, Inf, "uint64");
    fclose(timefid);
    path_usec1 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/usec1.txt';
    path_usec1 = strrep(path_usec1,'folder',folder);
    path_usec1 = strrep(path_usec1,'recording',recording_folder);
    timefid = fopen(path_usec1);
    usec1 = fread(timefid, Inf, "uint64");
    fclose(timefid);
    %for mic2
    path_sec2 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/sec2.txt';
    path_sec2 = strrep(path_sec2,'folder',folder);
    path_sec2 = strrep(path_sec2,'recording',recording_folder);
    timefid = fopen(path_sec2);
    timedata2 = fread(timefid, Inf, "uint64");
    fclose(timefid);
    path_usec2 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/usec2.txt';
    path_usec2 = strrep(path_usec2,'folder',folder);
    path_usec2 = strrep(path_usec2,'recording',recording_folder);
    timefid = fopen(path_usec2);
    usec2 = fread(timefid, Inf, "uint64");
    fclose(timefid);
    %for mic3
    path_sec3 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/sec3.txt';
    path_sec3 = strrep(path_sec3,'folder',folder);
    path_sec3 = strrep(path_sec3,'recording',recording_folder);
    timefid = fopen(path_sec3);
    timedata3 = fread(timefid, Inf, "uint64");
    fclose(timefid);
    path_usec3 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/usec3.txt';
    path_usec3 = strrep(path_usec3,'folder',folder);
    path_usec3 = strrep(path_usec3,'recording',recording_folder);
    timefid = fopen(path_usec3);
    usec3 = fread(timefid, Inf, "uint64");
    fclose(timefid);
    %for mic4
    path_sec4 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/sec4.txt';
    path_sec4 = strrep(path_sec4,'folder',folder);
    path_sec4 = strrep(path_sec4,'recording',recording_folder);
    timefid = fopen(path_sec4);
    timedata4 = fread(timefid, Inf, "uint64");
    fclose(timefid);
    path_usec4 = '/Users/emiliano/Desktop/Modelo Mapeo Acústico/recording/folder/usec4.txt';
    path_usec4 = strrep(path_usec4,'folder',folder);
    path_usec4 = strrep(path_usec4,'recording',recording_folder);
    timefid = fopen(path_usec4);
    usec4 = fread(timefid, Inf, "uint64");
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
    %% for mic1
    synchro_vector1 = timedata1((timedata1>1.6e15 & timedata1<1.7e15) | timedata1 < 1000);
    for i = 1:length(synchro_vector1)
       if synchro_vector1(i) < 10000
           synchro_vector1(i) = synchro_vector1(i-1) + (synchro_vector1(i-1)-synchro_vector1(i-2));
       end 
    end
    synchro_vector1 = synchro_vector1 + usec1;
    index_synchro1 = find((timedata1>1.6e15 & timedata1<1.7e15) | timedata1 < 1000);
    synchro_avr1 = zeros(1,length(index_synchro1));
    synchro_avr1(1) = index_synchro1(2)-index_synchro1(1);
    for i = 2:length(index_synchro1)
        synchro_avr1(i) = index_synchro1(i)-index_synchro1(i-1);
        lim_index = floor(mean2(synchro_avr1(i-1)));
        if (synchro_avr1(i) > lim_index*1.1) || (synchro_avr1(i) < lim_index*0.9)
            index_synchro1(i) = index_synchro1(i-1) + (index_synchro1(i-1)-index_synchro1(i-2));
        end
    end
    timedata1 = timedata1(timedata1>1.7e15 | timedata1<1.6e15);
    time_avrg = timedata1(2:end) - timedata1(1:end-1);
    for i = 2:length(timedata1)
        lim_index_time(i) = floor(mean2(time_avrg(1:i-1)));
        if (time_avrg(i-1) > lim_index_time(i)*1.1) || (time_avrg(i-1) < lim_index_time(i)*0.9)
            timedata1(i) = timedata1(i-1) + (timedata1(i-1)-timedata1(i-2));
        end
    end

    %% for mic2
    synchro_vector2 = timedata2((timedata2>1.6e15 & timedata2<1.7e15) | timedata2 < 1000);
    for i = 1:length(synchro_vector2)
        if synchro_vector2(i) < 1000000
           synchro_vector2(i) = synchro_vector2(i-1) + (synchro_vector2(i-1)-synchro_vector2(i-2));
        end 
    end
    synchro_vector2 = synchro_vector2 + usec2;
    index_synchro2 = find((timedata2>1.6e15 & timedata2<1.7e15) | timedata2 < 1000);
    synchro_avr2 = zeros(1,length(index_synchro2));
    synchro_avr2(1) = index_synchro2(2)-index_synchro2(1);
    for i = 2:length(index_synchro2)
        synchro_avr2(i) = index_synchro2(i)-index_synchro2(i-1);
        lim_index = floor(mean2(synchro_avr2(i-1)));
        if (synchro_avr2(i) > lim_index*1.1) || (synchro_avr2(i) < lim_index*0.9)
            index_synchro2(i) = index_synchro2(i-1) + (index_synchro2(i-1)-index_synchro2(i-2));
        end
    end
    timedata2 = timedata2(timedata2>1.7e15 | timedata2<1.6e15);
    time_avrg = timedata2(2:end) - timedata2(1:end-1);
    for i = 2:length(timedata2)
        lim_index_time(i) = floor(mean2(time_avrg(1:i-1)));
        if (time_avrg(i-1) > lim_index_time(i)*1.1) || (time_avrg(i-1) < lim_index_time(i)*0.9)
            timedata2(i) = timedata2(i-1) + (timedata2(i-1)-timedata2(i-2));
        end
    end
    
    %% for mic3
    synchro_vector3 = timedata3((timedata3>1.6e15 & timedata3<1.7e15) | timedata3 < 1000);
    for i = 1:length(synchro_vector3)
       if synchro_vector3(i) < 1000000
           synchro_vector3(i) = synchro_vector3(i-1) + (synchro_vector3(i-1)-synchro_vector3(i-2));
       end 
    end
    synchro_vector3 = synchro_vector3 + usec3;
    index_synchro3 = find((timedata3>1.6e15 & timedata3<1.7e15) | timedata3 < 1000);
    synchro_avr3 = zeros(1,length(index_synchro3));
    synchro_avr3(1) = index_synchro3(2)-index_synchro3(1);
    for i = 2:length(index_synchro3)
        synchro_avr3(i) = index_synchro3(i)-index_synchro3(i-1);
        lim_index = floor(mean2(synchro_avr3(i-1)));
        if (synchro_avr3(i) > lim_index*1.1) || (synchro_avr3(i) < lim_index*0.9)
            index_synchro3(i) = index_synchro3(i-1) + (index_synchro3(i-1)-index_synchro3(i-2));
        end
    end
    timedata3 = timedata3(timedata3>1.7e15 | timedata3<1.6e15);
    time_avrg = timedata3(2:end) - timedata3(1:end-1);
    for i = 2:length(timedata3)
        lim_index_time(i) = floor(mean2(time_avrg(1:i-1)));
        if (time_avrg(i-1) > lim_index_time(i)*1.1) || (time_avrg(i-1) < lim_index_time(i)*0.9)
            timedata3(i) = timedata3(i-1) + (timedata3(i-1)-timedata3(i-2));
        end
    end

    %% for mic4
    synchro_vector4 = timedata4((timedata4>1.6e15 & timedata4<1.7e15) | timedata4 < 1000);
    for i = 1:length(synchro_vector4)
       if synchro_vector4(i) < 1000000
           synchro_vector4(i) = synchro_vector4(i-1) +(synchro_vector4(i-1) - synchro_vector4(i-2));
       end 
    end
    synchro_vector4 = synchro_vector4 + usec4;
    index_synchro4 = find((timedata4>1.6e15 & timedata4<1.7e15) | timedata4 < 1000);
    synchro_avr4 = zeros(1,length(index_synchro4));
    synchro_avr4(1) = index_synchro4(2)-index_synchro4(1);
    for i = 2:length(index_synchro4)
        synchro_avr4(i) = index_synchro4(i)-index_synchro4(i-1);
        lim_index = floor(mean2(synchro_avr4(i-1)));
        if (synchro_avr4(i) > lim_index*1.1) || (synchro_avr4(i) < lim_index*0.9)
            index_synchro4(i) = index_synchro4(i-1) + (index_synchro4(i-1)-index_synchro4(i-2));
        end
    end
    timedata4 = timedata4(timedata4>1.7e15 | timedata4<1.6e15);
    time_avrg = timedata4(2:end) - timedata4(1:end-1);
    for i = 2:length(timedata4)
        lim_index_time(i) = floor(mean2(time_avrg(1:i-1)));
        if (time_avrg(i-1) > lim_index_time(i)*1.1) || (time_avrg(i-1) < lim_index_time(i)*0.9)
            timedata4(i) = timedata4(i-1) + (timedata4(i-1)-timedata4(i-2));
        end
    end

    %% Adjusting timestamps
    %for mic1
    for i = 1:length(index_synchro1) 
      this_index_ini = index_synchro1(i);

      this_index_fin = 0;
      if i == length(index_synchro1)
        this_index_fin = length(timedata1);
        timedata1(this_index_ini:this_index_fin) = timedata1(this_index_ini:this_index_fin) - timedata1(this_index_ini) + synchro_vector1(i);
      else
        this_index_fin = index_synchro1(i+1)-1;

        this_diff1 = (timedata1(this_index_fin+1) - timedata1(this_index_ini) + synchro_vector1(i)) - synchro_vector1(i+1);

        this_l = this_index_fin-this_index_ini;

        synch_corrected1 = ((0:(this_l+1))'/(this_l+1))*this_diff1;

        timedata1(this_index_ini:this_index_fin) = timedata1(this_index_ini:this_index_fin) - timedata1(this_index_ini) + synchro_vector1(i) - synch_corrected1(1:end-1);
      end
    end
    %for mic2
    for i = 1:length(index_synchro2) 
      this_index_ini = index_synchro2(i);
      this_index_fin = 0;
      if i == length(index_synchro2)
        this_index_fin = length(timedata2);
        timedata2(this_index_ini:this_index_fin) = timedata2(this_index_ini:this_index_fin) - timedata2(this_index_ini) + synchro_vector2(i);
      else
        this_index_fin = index_synchro2(i+1)-1;

        this_diff2 = (timedata2(this_index_fin+1) - timedata2(this_index_ini) + synchro_vector2(i)) - synchro_vector2(i+1);

        this_l = this_index_fin-this_index_ini;

        synch_corrected2 = ((0:(this_l+1))'/(this_l+1))*this_diff2;

        timedata2(this_index_ini:this_index_fin) = timedata2(this_index_ini:this_index_fin) - timedata2(this_index_ini) + synchro_vector2(i) - synch_corrected2(1:end-1);
      end
    end
    %for mic3
    for i = 1:length(index_synchro3) 
      this_index_ini = index_synchro3(i);

      this_index_fin = 0;
      if i == length(index_synchro3)
        this_index_fin = length(timedata3);
        timedata3(this_index_ini:this_index_fin) = timedata3(this_index_ini:this_index_fin) - timedata3(this_index_ini) + synchro_vector3(i);
      else
        this_index_fin = index_synchro3(i+1)-1;

        this_diff3 = (timedata3(this_index_fin+1) - timedata3(this_index_ini) + synchro_vector3(i)) - synchro_vector3(i+1);

        this_l = this_index_fin-this_index_ini;

        synch_corrected3 = ((0:(this_l+1))'/(this_l+1))*this_diff3;

        timedata3(this_index_ini:this_index_fin) = timedata3(this_index_ini:this_index_fin) - timedata3(this_index_ini) + synchro_vector3(i) - synch_corrected3(1:end-1);
      end
    end
    %for mic4
    for i = 1:length(index_synchro4) 
      this_index_ini = index_synchro4(i);
      this_index_fin = 0;
      if i == length(index_synchro4)
        this_index_fin = length(timedata4);
        timedata4(this_index_ini:this_index_fin) = timedata4(this_index_ini:this_index_fin) - timedata4(this_index_ini) + synchro_vector4(i);
      else
        this_index_fin = index_synchro4(i+1)-1;

        this_diff4 = (timedata4(this_index_fin+1) - timedata4(this_index_ini) + synchro_vector4(i)) - synchro_vector4(i+1);

        this_l = this_index_fin-this_index_ini;

        synch_corrected4 = ((0:(this_l+1))'/(this_l+1))*this_diff4;

        timedata4(this_index_ini:this_index_fin) = timedata4(this_index_ini:this_index_fin) - timedata4(this_index_ini) + synchro_vector4(i) - synch_corrected4(1:end-1);
      end
    end

    %new timestamp vectors
    timedata1 = timedata1(1:end)/1000000; % so that timedata is in seconds
    timedata2 = timedata2(1:end)/1000000; % so that timedata is in seconds
    timedata3 = timedata3(1:end)/1000000; % so that timedata is in seconds
    timedata4 = timedata4(1:end)/1000000; % so that timedata is in seconds
    

    %% Creating time vectors
    %for mic1
    if index_synchro1(1) ~= 1
        new_dataIndex1 = (index_synchro1(1)-1)*1024+1;
        timedata1 = timedata1(index_synchro1(1):end);
        data1 = data1(new_dataIndex1:end);
    end
    t_length1 = length(timedata1) * nframes;
    t1 = zeros(t_length1,1);
    
    
    for i = 1:length(timedata1)
      t_start = timedata1(i);
      if(i ~= length(timedata1))
        t_end = timedata1(i+1);
      else
        t_end = t_start + (timedata1(i)-timedata1(i-1));
      end

      t_delta = (t_end-t_start)/nframes;
      t_values = t_start:t_delta:(t_end-t_delta);

      t_range = (((i-1)*nframes)+1):i*nframes;

      t1(t_range) = t_values;
    end
    
    %for mic2
    if index_synchro2(1) ~= 1
        new_dataIndex2 = (index_synchro2(1)-1)*1024+1;
        timedata2 = timedata2(index_synchro2(1):end);
        data2 = data2(new_dataIndex2:end);
    end
    t_length2 = length(timedata2) * nframes;
    t2 = zeros(t_length2,1);
    
    
    for i = 1:length(timedata2)
      t_start = timedata2(i);
      if(i ~= length(timedata2))
        t_end = timedata2(i+1);
      else
        t_end = t_start + (timedata2(i)-timedata2(i-1));
      end

      t_delta = (t_end-t_start)/nframes;
      t_values = t_start:t_delta:(t_end-t_delta);

      t_range = (((i-1)*nframes)+1):i*nframes;

      t2(t_range) = t_values;
    end

    %for mic3
    if index_synchro3(1) ~= 1
        new_dataIndex3 = (index_synchro3(1)-1)*1024+1;
        timedata3 = timedata3(index_synchro3(1):end);
        data3 = data3(new_dataIndex3:end);
    end
    t_length3 = length(timedata3) * nframes;
    t3 = zeros(t_length3,1);
    
    for i = 1:length(timedata3)
      t_start = timedata3(i);
      if(i ~= length(timedata3))
        t_end = timedata3(i+1);
      else
        t_end = t_start + (timedata3(i)-timedata3(i-1));
      end

      t_delta = (t_end-t_start)/nframes;
      t_values = t_start:t_delta:(t_end-t_delta);

      t_range = (((i-1)*nframes)+1):i*nframes;

      t3(t_range) = t_values;
    end

    %for mic4
    if index_synchro4(1) ~= 1
        new_dataIndex4 = (index_synchro4(1)-1)*1024+1;
        timedata4 = timedata4(index_synchro4(1):end);
        data4 = data4(new_dataIndex4:end);
    end
    t_length4 = length(timedata4) * nframes;
    t4 = zeros(t_length4,1);
    for i = 1:length(timedata4)
      t_start = timedata4(i);
      if(i ~= length(timedata4))
        t_end = timedata4(i+1);
      else
        t_end = t_start + (timedata4(i)-timedata4(i-1));
      end
      t_delta = (t_end-t_start)/nframes;
      t_values = t_start:t_delta:(t_end-t_delta);

      t_range = (((i-1)*nframes)+1):i*nframes;

      t4(t_range) = t_values;
    end


    %% Results
    full_times = zeros(max([length(t1) length(t2) length(t3) length(t4)]),4);
    full_times(1:length(t1),1) = t1;
    full_times(1:length(t2),2) = t2;
    full_times(1:length(t3),3) = t3;
    full_times(1:length(t4),4) = t4;
    full_data = zeros(max([length(data1) length(data2) length(data3) length(data4)]),4);
    full_data(1:length(data1),1) = data1;
    full_data(1:length(data2),2) = data2;
    full_data(1:length(data3),3) = data3;
    full_data(1:length(data4),4) = data4;
    begin_end = zeros(4,2);
    begin_end(1,1) = t1(1);
    begin_end(2,1) = t2(1);
    begin_end(3,1) = t3(1);
    begin_end(4,1) = t4(1);
    begin_end(1,2) = t1(end);
    begin_end(2,2) = t2(end);
    begin_end(3,2) = t3(end);
    begin_end(4,2) = t4(end);
end 