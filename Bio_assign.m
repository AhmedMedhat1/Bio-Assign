close all
clc
%%  

%Plotting the input raw signal
yt=amp'; %amp is the vector containg the input data .. yt is the same datain a row instead of column

fsample=256;
N=2000;
t=(1:N);

plot(t,yt);
title('Input Signal (Time domain)')
xlabel('Time (sec)')
ylabel('Amplitude')

%%

% Design of the filter
forder=50;
fcutoff=45;
wn=fcutoff/(fsample/2); % Normalized cutoff frequency % Max=1 =Nyquist Rate = Half Sample Rate OR pi in radian/Sample
[b,a]=butter(forder,wn,'low');
[H,w]=freqz(b,a,512); %what's this function
w=w*fsample/(2*pi); %transform from w as radian/Sample to f as mulitple of fs
figure
plot(w,20*log10(abs(H))); %dB scale
grid on
title('Filter frequency response')
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')


%%

% Applying the filter on the given signal
clc
ytfiltered=filter(b,a,yt);%first method
figure
subplot(2,1,1); plot(t,yt)
title('Unfiltered signal')
xlabel('Time (sec)')

subplot(2,1,2); plot(t,ytfiltered);
title('Filtered signal')
xlabel('Time (sec)')


%%

%Differentiation of the filtered signal
xt=ytfiltered;
for k=3:N-2
y(k-2)=(1/8)*fsample*(-xt(k-1)-2*xt(k-1)+2*xt(k+1)+xt(k+2));
end
y(N)=0;
y(N-1)=0;
figure
subplot(2,1,1); 
plot(t,y);

title('Output signal (Time domain)')
xlabel('Time (sec)')
yff=abs(fft(y));
subplot(2,1,2); 
plot(f,yff(1:N/2));
title('Frequency domain')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

%%

%Sqauring of the the signal after differentiation
ysq=power(y,2);
plot(t,ysq)
title('After Squaring (Time domain)')
xlabel('Time (sec)')



%%
%Smoothing squared signal using moving average window
ws=25; %windows
ys=0;
ysi=ysq;
for k=N:-1:ws
    for j=1:ws-1
    
        ysi(k)=ysi(k)+ysi(k-j);
    end
    ys(k)=(ysi(k))/ws;
end
    plot(t,ys);
    title('At Window Size 25');

%%
%Threshold
peak=max(ys);
threshold =0.6*peak;
plot([1 2000], [threshold threshold]);
hold on
plot (t,ys);
hold off
% l=1;
% pks=zeros([1 N]);
% for i=1:N
%     if ys(i)>threshold
%         if ys(i)> max(pks)
%             pks(i)=ys(i);
%             
%         end
%     else plot(t(i),max(pks),'*');
%          pks=zeros([1 N]);
%     
% 
%             




    






















