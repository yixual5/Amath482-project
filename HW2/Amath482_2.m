clc;close all;clear all
%part 1
load handel
v=y';
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('signal of Interest, v(n)')

n=length(v); %how many points do we have
t = (1:length(v))/Fs; % follow what instruction did to define t
L=length(v)/Fs; %around 9 seconds
k = (2*pi/L)*[0:(n-1)/2 -(n-1)/2:-1]; % rescale to L domain
ks = fftshift(k); %easy for ploting.


%% Spectrograms for varying window sizes(5,1,0.2 with dt = 0.1)
figure(2)

a_vec = [10 1 0.1];
for i = 1:length(a_vec)
    a = a_vec(i);
    tslide=0:0.1:L;
    Vgt_spec = zeros(length(tslide),n);%gaussian window
    Vmgt_spec = zeros(length(tslide),n);%mexican hat wavelet
    Vsgt_spec = zeros(length(tslide),n);%step function
    for j=1:length(tslide)
        g=exp(-a*(t-tslide(j)).^2);
        m = 2/(pi^1/4*sqrt(3*a))*(1-((t-tslide(j))/a).^2).*exp(-((t-tslide(j))/a).^2);
        s = zeros(1,length(v));
        temp = find(tslide(j)-a<t & t< tslide(j)+a);
        s(1,temp(1):temp(end)) = ones(1,length(temp));
        Vg=g.*v;
        Vmg=m.*v;
        Vsg=s.*v;
        Vgt=fft(Vg);
        Vmgt=fft(Vmg);
        Vsgt=fft(Vsg);
        Vgt_spec(j,:) = fftshift(abs(Vgt)); % match with ks
        Vmgt_spec(j,:) = fftshift(abs(Vmgt));
        Vsgt_spec(j,:) = fftshift(abs(Vsgt));
    end
    
    figure(i)
    pcolor(tslide,ks,Vgt_spec.'), 
    shading interp 
    title(['Gabor filtering(Gaussian window) with a = ',num2str(a)],'Fontsize',16)
    xlabel('Time(s)')
    ylabel('Frequency(¦Ø)')
    colormap(hot) 
    
    figure(i+3)
    pcolor(tslide,ks,Vmgt_spec.'), 
    shading interp 
    title(['Gabor filtering(Mexican hat window) with a = ',num2str(a)],'Fontsize',16)
    xlabel('Time(s)')
    ylabel('Frequency(¦Ø)')
    colormap(hot) 
    
    figure(i+6)
    pcolor(tslide,ks,Vsgt_spec.'), 
    shading interp 
    title(['Gabor filtering(step function) with a = ',num2str(a)],'Fontsize',16)
    xlabel('Time(s)')
    ylabel('Frequency(¦Ø)')
    colormap(hot) 
    
end

%% Spectrograms for varying dt[0.01 1] with a =1
a=1;
dt_vec = [0.01 1];
for jj = 1:length(dt_vec)
    dt = dt_vec(jj);
    tslide=0:dt:L;
    Vgt_spec = zeros(length(tslide),n);
    for j=1:length(tslide)
        g=exp(-a*(t-tslide(j)).^2); 
        Vg=g.*v; 
        Vgt=fft(Vg); 
        Vgt_spec(j,:) = fftshift(abs(Vgt)); 
    end
    
    figure()
    pcolor(tslide,ks,Vgt_spec.'), 
    shading interp 
    title(['dt = ',num2str(dt)],'Fontsize',16)
    xlabel('Time(s)')
    ylabel('Frequency(¦Ø)')
    colormap(hot) 
end

%% part 2
%% piano
[y,Fs] = audioread('music1.wav');
tr_piano=length(y)/Fs; % record time in seconds
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)');

n=length(y); %how many points do we have
t = (1:length(y))/Fs; % follow what instruction did to define t
L=tr_piano; %around 16 seconds
k = (2*pi/L)*[0:(n-1)/2 -n/2:-1]; % rescale to L domain
ks = fftshift(k); %easy for ploting
v = y';

a=20;
dt=0.1;
tslide=0:dt:L;
music_score =  zeros(1,length(tslide));
Vgt_spec = zeros(length(tslide),n);%gaussian window
    for j=1:length(tslide)
        g=exp(-a*(t-tslide(j)).^2);
        Vg=g.*v;
        Vgt=fft(Vg);
        Vgt_spec(j,:) = fftshift(abs(Vgt))/max(abs(Vgt)); % match with ks
        [mxv,idx] = max(abs(Vgt));
        music_score(j) = k(idx)/(2*pi); % no shifting of k and convert it to Hz
    end
plot(tslide,music_score)
title('Frequency for the piano without overtones')
xlabel('Time(sec)')
ylabel('Frequency(Hz')
 figure()
    pcolor(tslide,ks/(2*pi),Vgt_spec.'), 
    shading interp 
    title('spectrograms of piano','Fontsize',16)
    xlabel('Time(s)')
    ylabel('Frequency(¦Ø)')
    set(gca,'Ylim',[200 500],'Fontsize',16) 
    colormap(hot) 
    
%% recorder
figure(2)
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs; % record time in seconds
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
n=length(y); %how many points do we have
t = (1:length(y))/Fs; % follow what instruction did to define t
L=tr_rec; %around 14 seconds
k = (2*pi/L)*[0:(n-1)/2 -n/2:-1]; % rescale to L domain
ks = fftshift(k); %easy for ploting
v = y';

a=20;
dt=0.1;
tslide=0:dt:L;
music_score =  zeros(1,length(tslide));
Vgt_spec = zeros(length(tslide),n);%gaussian window
    for j=1:length(tslide)
        g=exp(-a*(t-tslide(j)).^2);
        Vg=g.*v;
        Vgt=fft(Vg);
        Vgt_spec(j,:) = fftshift(abs(Vgt))/max(abs(Vgt)); % match with ks
        [mxv,idx] = max(abs(Vgt));
        music_score(j) = k(idx)/(2*pi); % no shifting of k and convert it to Hz
    end
plot(tslide,music_score)
title('Frequency for the recorder without overtones')
xlabel('Time(sec)')
ylabel('Frequency(Hz')
 figure()
    pcolor(tslide,ks/(2*pi),Vgt_spec.'), 
    shading interp 
    title('spectrograms of recorder','Fontsize',16)
    xlabel('Time(s)')
    ylabel('Frequency(¦Ø)')
    set(gca,'Ylim',[750 1100],'Fontsize',16) 
    
    colormap(hot) 