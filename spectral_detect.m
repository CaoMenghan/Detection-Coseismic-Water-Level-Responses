function [f_tsecs] = spectral_detect(data_cut, t_cut, wtime_s, f_min, f_max, window, noverlap, fs, z_th, duration)
%Detection in the frequency domain, including two parts.
%Step 1: calculate the power spectral density for each window;
%Step 2: determine whether there is a response
f_tsecs=[];
n=length(data_cut);
f_len=window/2+1;% if nfft(2^(ceil(log2(length(window))))) is odd,len=(nfft+1)/2 
f=linspace(f_min,f_max,f_len);
[s,~,~,p]=spectrogram(data_cut,window,noverlap,f,fs);
t=linspace(t_cut(1),t_cut(n-window+1),size(s,2));

[raw,N]=size(s);
sml=zeros(1,N);
for j=1:N
    for i=1:raw
        sml(1,j)=sml(1,j)+p(i,j);
    end
end
miu=mean(sml);
sigma=std(sml);

j=1;
z=[];
for k=1:length(sml)
    z(j)=(sml(k)-miu)/sigma;
    j=j+1;
end
j=1;result=[];
for k=1:length(z)
    tdiff=etime(datevec(t(k)),datevec(wtime_s));
    if z(k)>=z_th && abs(tdiff)<=duration
        result(j)=k;
        j=j+1;
    end
end
if size(result,2)>1
    for i=1:size(result,2)
            new=[t(result(i)),t(result(i))];
            f_tsecs=[f_tsecs;new];
    end
else
    if logical(result)
        f_tsecs=[f_tsecs;[t(result),t(result)]];
    end
end
end