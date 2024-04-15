% Automatic detection program for coseismic water level oscillations
% Set Parameters
clc;clear all
% Parameters for STA/LTA
sta_seconds = 5; % STA time window (Second) water:5 earthquake:5
lta_seconds =6000; % LTA time window (Second) water:6000 earthquake:6000
water_on =  13; % Event triggers "ON" with STA/LTA ratio exceeds water:13 earthquake:30
water_off = 8; % Event triggers "OFF" when STA/LTA ratio drops below  water:8 earthquake:15
etq_on_e = 30;
etq_off_e = 15;
minimum_event_duration_seconds = 2; % Trigger must be on at least 2 secs
water_detection_params = [1 sta_seconds lta_seconds water_on water_off ...
    minimum_event_duration_seconds];
etq_detection_params = [40 sta_seconds lta_seconds etq_on etq_off ...
    minimum_event_duration_seconds];

% Parameters for Spectral detection
duration = 8000;
window = 1200;
noverlap = window/2;

f_tsecs=[];
f_min=0.02;
f_max=0.075;
fs=1;
z_th=3;

% Input File 

path_e = fullfile('input path');
namelist_e = dir([path, 'input_filename']);

eventname = namelist_e.name;
e = load([path, eventname]);

path_w = fullfile('input path');
namelist_w = dir([path, 'input_filename']);

eventname = namelist_w.name;
w = load([path, eventname]);

% Output File

out_file = fullfile('output path');

% Earthquake Detection
% run the STA/LTA detector. lta_mode = 'frozen' means the LTA stops
% updating when trigger is "ON".
%--------------------------------------------------------------------------
[etq_sta,etq_LTA,etq_sta_to_lta, etq_ta_num] = watersta_lta(e, 'edp', etq_detection_params, ...
    'lta_mode', 'frozen');
etq_ta_secs=zeros(size(etq_ta_num));
for j=1:size(etq_ta_num,1)
    etq_ta_secs(j,1)=e.t(etq_ta_num(j,1));
    etq_ta_secs(j,2)=e.t(etq_ta_num(j,2));
end

save([out_file,'etq_trigger_outfile.mat'],'etq_sta','etq_LTA','etq_sta_to_lta','etq_ta_secs','etq_ta_num')

% Water Detection
[w_sta,w_LTA,w_sta_to_lta, w_ta_num] = watersta_lta(w, 'edp', water_detection_params, ...
    'lta_mode', 'frozen');
w_ta_secs=zeros(size(w_ta_num));
for j=1:size(w_ta_num,1)
    w_ta_secs(j,1)=w.t(w_ta_num(j,1));
    w_ta_secs(j,2)=w.t(w_ta_num(j,2));
end

% Deletion of misidentification based on seismic detection
record=[];
for k=1:size(w_ta_secs,1)
    for m=1:size(etq_ta_secs,1)
        time_w=datevec(w_ta_secs(k,1));  
        time_e=datevec(etq_ta_secs(m,1));
        day=time_w(3);
        if time_e(3)==day
            difft=etime(time_w,time_e);
            if abs(difft)<=2400 % It is considered to identify the hydrological response when time differences are less than 2400s 
               record=[record;w_ta_secs(k,:)];
               break;
            end
        end
    end
end
w_ta_secs=record;

save([out_file,'w_trigger_outfile.mat'],'w_sta','w_LTA','w_sta_to_lta','w_ta_secs');


% Spectral Detection
eresult=load('input_triggerfile');
wresult=load('input_trigerfile');

%---Find the detected seismic events with water level responses
eta_secs=eresult.ta_secs;
eB=[];
for i=1:size(eresult.ta_secs,1)
    A=datevec(eresult.ta_secs(i,1));
    eB=[eB;A];
end
wC=[];
for i=1:size(wresult.ta_secs,1)
    A=datevec(wresult.ta_secs(i,1));
    wC=[wC;A];
end
for i=1:size(eta_secs,1)
    day=eB(i,3);
    for j=1:size(wresult.ta_secs,1)
        if wC(j,3)==day
            e=etime(wC(j,:),eB(i,:));
            if abs(e)<=duration
                eresult.ta_secs(i,:)=100;
            end
        end
    end
end
eresult.ta_secs(eresult.ta_secs(:,1)==100,:)=[];
% ---cut the data window with no water level response but with a detected seismic
% event
for q=1:size(eresult.ta_secs,1)
    if eresult.ta_secs(q,1)<w.t(1)
        continue;
    end
    tdiff=abs(w.t-eresult.ta_secs(q,1));
    count_s=find(tdiff==min(tdiff));
    wtime_s=w.t(count_s);
    tdiff=abs(w.t-eresult.ta_secs(q,2));
    count_end=find(tdiff==min(tdiff));
    
    if (count_end+duration)>length(w.t)
        t_cut=w.t(count_s-duration-1000:end);
        data_cut=w.Pre_filt(count_s-duration-1000:end);
    else
        t_cut=w.t(count_s-duration:count_end+duration);
        data_cut=w.Pre_filt(count_s-duration:count_end+duration);
    end
    f_tsecs=spectral_detect(data_cut,t_cut,wtime_s,f_min,f_max,window,noverlap,fs,z_th,duration);
end
      
save([out_file,'w_trigger_outfile.mat'],'f_tsecs','-append');