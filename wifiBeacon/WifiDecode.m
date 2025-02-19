% % Open the file containing the received samples
% f2 = fopen('rx.dat', 'rb');
% 
% % read data from the file
% tmp = fread(f2, 'float32');
% 
% % close the file
% fclose(f2);
% 
% % since the USRP stores the data in an interleaved fashion
% % with real followed by imaginary samples 
% % create a vector of half the length of the received values to store the
% % data. Make every other sample the real part and the remaining samples the
% % imaginary part
% rx = zeros(length(tmp)/2,1);
% rx = tmp(1:2:end)+j*tmp(2:2:end);
% 
rx = rx_new;

%%

%real:
LTS=[0,0,0,0,0,1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 0, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1,0,0,0,0,0,0];

%LTS for CFO and channel estimation
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

STS=[0, 0, 1+j, 0,0,0,-1-j, 0,0,0,1+j,0,0,0,-1-j, 0,0,0,1+j,0,0,0,0,0,0,0,-1-j,0,0,0,1+j,0,0,0,1+j,0,0,0,1+j,0,0,0,1+j,0,0];

%%
% Cross correlates to find start of lts
[Ryx, lags] = xcorr(rx, [lts_t(33:end),lts_t,lts_t]);
[mm, ii] = max(abs(Ryx));
idx = lags(ii) + 1 + 32;

% ltx_rx is the two lts symbols. We are excluding the first half lts
% symbol because it is not necessary for our calculations
lts_rx = rx(idx:idx+length(lts_t)*2 - 1);  

plot(lags, abs(Ryx))

%Schmidl-Cox Algorithm
angle_sum  = 0;
for n = 1: 64
   angle_sum = angle_sum + angle(lts_rx(64+n)/lts_rx(n)); 
end

f_delta_hat = angle_sum/(64*64);

% Samples 1-128 are the lts, samples 129-208 are the signal field
rx_phase_corrected = zeros(10000,1);

% Correct phase offset
for k = 1:10000
    rx_phase_corrected(k) = rx(idx+k-1)/exp(1j*f_delta_hat*(k-1));
end

%%
% Channel Estimation
lts_f = [realmax 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 realmax realmax realmax realmax realmax realmax realmax realmax realmax realmax realmax 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

headers_rx = zeros(64,2);
HEADERS_rx = zeros(64,2);  % FFT of headers_rx
H_ests = zeros(64,2);

for n = 1:2
    headers_rx(:,n) = rx_phase_corrected((n-1)*64+1:(n-1)*64+64);
    HEADERS_rx(:,n) = fft(headers_rx(:,n));
    H_ests(:,n) = HEADERS_rx(:,n)./(lts_f');
end

H_est = mean(H_ests,2);

%%

% Isolate the 64 sample signal field from the 80 sample circular
% convolution
signal_field_rx = rx_phase_corrected(145:208);

% Convert to frequency domain
SIGNAL_FIELD_rx = fft(signal_field_rx);

% Correct for channel
SIGNAL_FIELD_CORRECTED = SIGNAL_FIELD_rx./H_est;

% Sets unreasonably large values to 0, these subcarriers don't hold data
for k=1:64
   if(abs(SIGNAL_FIELD_CORRECTED(k)) > 20)
       SIGNAL_FIELD_CORRECTED(k) = 0;
   end
end

% Snaps subcarriers to 1 and -1
SIGNAL_FIELD_SIGNED = sign(real(SIGNAL_FIELD_CORRECTED));

plot(real(X_est),imag(X_est),'.');
title('Signal Field');
plot(real(SIGNAL_FIELD_CORRECTED));
figure;
plot(real(SIGNAL_FIELD_SIGNED));
title('Signal Field Adjusted');

msg_block_num = 15;

ys_rx = zeros(64,msg_block_num);
Ys_rx = zeros(64,msg_block_num); % FFT of ys_rx
X_est = [];
f_drift = zeros(msg_block_num,1);

for n = 1:msg_block_num
   ys_rx(:,n) = rx_phase_corrected((n-1)*80+240+17:240+80*n); 
   Ys_rx(:,n) = fft(ys_rx(:,n));
   ys_H(:,n)  = [Ys_rx(:,n)./(H_est)];
   
   % Calculates drift based off of pilot subcarriers at 7,21,44, and 58
   f_drift(n) = 0.25*(angle(ys_H(7,n)/1)+angle(ys_H(21,n)/-1)+angle(ys_H(44,n)/1)+angle(ys_H(58,n)/1));
   X_est = [X_est; ys_H(:,n)./exp(j*f_drift(n))];
   
end

for k=1:length(X_est)
   if abs(X_est(k))>20
       X_est(k) = NaN;
   end
end

figure;
plot(real(X_est),imag(X_est),'.');
xlabel('I');
ylabel('Q');
title('Constellation Plot');

figure;
plot(abs(H_est));
title('Channel Estimation');