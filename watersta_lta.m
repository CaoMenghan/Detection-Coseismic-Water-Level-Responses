function [sta, LTA, sta_to_lta, ta_secs] = watersta_lta(w,varargin)
    l_v = length(w.Pre_filt); % Length of time series
    w.Pre_filt=reshape(w.Pre_filt,length(w.Pre_filt),[]);
 
    % Check varargin 
    nv = numel(varargin);

    for p = 1:2:nv-1
        v1 = varargin{p};
        v2 = varargin{p+1};
        switch lower(v1)
            case 'edp'
                if isnumeric(v2) && numel(v2) == 6
                    Fs = v2(1);
                    l_sta = round(v2(2)*Fs);    % STA window length
                    l_lta = round(v2(3)*Fs);    % LTA window length
                    th_on = v2(4);       % Trigger on theshold
                    th_off = v2(5);      % Trigger off threshold
                    minimum_duration_days = v2(6);  % Minimum event duration v2(5)/86400
                else
                    error('STA_LTA: Wrong format for input ''edp''')
                end
            case 'lta_mode'
                switch lower(v2)
                    case {'freeze','frozen'}
                        lta_mode = 'frozen';
                    case {'continue','continuous'}
                        lta_mode = 'continuous';
                    case {'grow','growing'}
                        lta_mode = 'grow';
                    otherwise
                      error('STA_LTA: Wrong format for input ''lta_mode''')
                end
            otherwise
                warning('property name not recognized')
                
        end
    end
    % Initialize waveform data
    y = (w.Pre_filt).*(w.Pre_filt);
    Number=1:1:l_v;
    eventnum = 0;
    EVENT_ON = false;
    RECORD=false;
    
    eventstart=0;
    eventend=0;
    index1=0;
    wrong_end=0;
    
    trig_array = [];
    sta = ones(size(y));
    lta = sta;
    sta_to_lta = sta;

    
    %-----Process for data '0'---
    % index=(w.Pre_filt==0);
    % C=diff(double(index));
    % index2=(find(C==-1));
    % index3=(find(C==1));
    % if length(index3)==length(index2)+1
    %     index2=[index2;index3(end)];
    % end
    % dif=index2-index3;
    % index5=index2(dif<=(l_lta+l_sta)*0.2);
    % for i=1:length(index4)
    %     index(index4(i)+1:index5(i))=0;
    % end
    % D=diff(double(index));
    % index6=(D==-1);
    % index6=[index6;index6(end)];


    
    
    
      
    %-------------------initialize for first l_lta samples-----------------
    sta(1:l_sta) = cumsum(y(1:l_sta))/l_sta;
    lta(1:l_lta) = cumsum(y(1:l_lta))/l_lta;
    for count=l_lta+1:l_lta+l_sta
        lta(count)=(1/l_lta)*y(count)+(1-1/l_lta)*lta(count-1);
    end
    lta(count)=lta(count-l_sta);
    
    for count = l_sta+1:l_lta+l_sta
        sta(count)=(1/l_sta)*y(count)+(1-1/l_sta)*sta(count-1);% recuisive STA/LTA
    end
    sta_to_lta(1:l_lta+l_sta) = sta(1:l_lta+l_sta)./lta(1:l_lta+l_sta);
    
    for count=l_lta+l_sta+1:length(y) 
        lta(count)=(1/l_lta)*y(count-l_sta)+(1-1/l_lta)*lta(count-1);
         sta(count)=(1/l_sta)*y(count)+(1-1/l_sta)*sta(count-1);
    end
    LTA=lta;
 
    for count=l_lta+l_sta+1:length(y)

        if EVENT_ON && strcmp(lta_mode,'frozen') % When an event occurs and it is frozen, make the 'lta' value equal to the frozen value
            LTA(count) = lta_freeze_level; % freeze LTA is event is triggering
        else
            if eventnum>0 && (count-index1)<1.5*(l_sta+l_lta)&&strcmp(lta_mode,'frozen')  
              LTA(count)=lta_freeze_level;
            end
        end    
        %---------------------Process for data '0'----------------------
       %      if index(count)==1
       %          LTA(count)=LTA(count-1);
       %      end
       %      if index6(count)==1
       %          wrong_end=count;
       %      end
       %      if count-wrong_end<(l_sta+l_lta)*0.8 % The effect of data '0' is considered negligible when the length of '0' is longer than 20% of window
       %          LTA(count)=LTA(count-1);
       %      end
       %--------------------------------------------------------------------------
       sta_to_lta(count) = sta(count)/LTA(count);     
        if ~EVENT_ON && ~RECORD&&sta_to_lta(count) >=th_off
           lta_freeze_level = LTA(count);
           RECORD = false;
        end
        if ~EVENT_ON&&RECORD&&sta_to_lta(count)<=th_off
            RECORD = false;
        end

        if ~EVENT_ON && sta_to_lta(count) >= th_on 
            EVENT_ON = true;
            eventstart = Number(count); 
        end
        if EVENT_ON && ((sta_to_lta(count) <= th_off) || count == length(y))
            EVENT_ON = false;
            RECORD = false;
            eventend = Number(count); 
            if strcmp(lta_mode,'frozen') % unfreeze the lta
                LTA(count)=lta_freeze_level;   
            end
            index1=Number(count);
            if (eventend-eventstart)>=minimum_duration_days*Fs
                eventnum = eventnum + 1; 
                trig_array(eventnum, 1) = eventstart;
                trig_array(eventnum, 2) = eventend;
                eventstart = 0;
                eventend = 0;
            end 
        end
    end
   
    
    % Combining events
    if eventnum>0
    eventnum_end=eventnum;
    newtrig_array=[];
    newtrig_array(1,1)=trig_array(1,1);
    newtrig_array(1,2)=trig_array(1,2);
    k=1;
    for j=1:1:eventnum-1
    if (trig_array(j+1, 1)-trig_array(j,2))<=1800*Fs
        newtrig_array(k,2)=trig_array(j+1,2);
        eventnum_end=eventnum_end-1;
    else
    k=k+1;
    newtrig_array(k,1)=trig_array(j+1,1);
    newtrig_array(k,2)=trig_array(j+1,2);
    end
    end
      
    ta_secs = newtrig_array-Number(1);
    else    

    ta_secs = (trig_array-Number(1)); 
    end
end   
    