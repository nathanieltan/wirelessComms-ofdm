% Open the file containing the received samples
f2 = fopen('rx.dat', 'rb');

% read data from the file
tmp = fread(f2, 'float32');

% close the file
fclose(f2);

% since the USRP stores the data in an interleaved fashion
% with real followed by imaginary samples 
% create a vector of half the length of the received values to store the
% data. Make every other sample the real part and the remaining samples the
% imaginary part
rx = zeros(length(tmp)/2,1);
rx = tmp(1:2:end)+j*tmp(2:2:end);

HEADER_tx = header;

% Cross correlates to find start of lts
[Ryx, lags] = xcorr(rx, [timingHeader;timingHeader;timingHeader]);
[mm, ii] = max(abs(Ryx));
idx = lags(ii) + 1 - 64;
lts_rx = rx(idx:idx+length(timingHeader)*3 - 1);

%Schmidl-Cox Algorithm
angle_sum  = 0;
for n = 1: 64
   angle_sum = angle_sum + angle(lts_rx(128+n)/lts_rx(64+n)); 
end

f_delta_hat = angle_sum/(64*64);

rx_phase_corrected = zeros(length(rx),1);
% Correct phase offset
for k = 1:length(rx(idx:end));
    rx_phase_corrected(idx+k-1) = rx(idx+k-1)/exp(1j*f_delta_hat*(k-1));
end


% Channel Estimation
headers_start = idx+length(timingHeader)*3;

cyclic_headers_rx = rx_phase_corrected(headers_start:headers_start + 7999);
headers_rx = zeros(64,100);
HEADERS_rx = zeros(64,100);  % FFT of headers_rx
H_ests = zeros(64,100);

for n = 1:10
    headers_rx(:,n) = cyclic_headers_rx(80*(n-1)+17:80*n); 
    HEADERS_rx(:,n) = fft(headers_rx(:,n));
    H_ests(:,n) = HEADERS_rx(:,n)./HEADER_tx;
end

H_est = mean(H_ests,2);

% Data Estimation
msg_block_num = length(message)/64;

y_start = headers_start + 800;

cyclic_y_rx = rx_phase_corrected(y_start:y_start+msg_block_num*80-1);
ys_rx = zeros(64,msg_block_num);
Ys_rx = zeros(64,msg_block_num); % FFT of ys_rx
X_est = [];
f_drift = zeros(msg_block_num,1);

for n = 1:msg_block_num
   ys_rx(:,n) = cyclic_y_rx((n-1)*80+17:80*n); 
   Ys_rx(:,n) = fft(ys_rx(:,n));
   f_drift(n) = 0.25*(angle(Ys_rx(7,n)/message(64*(n-1)+7))+angle(Ys_rx(21,n)/message(64*(n-1)+21))+angle(Ys_rx(44,n)/message(64*(n-1)+44))+angle(Ys_rx(58,n)/message(64*(n-1)+58)));
   Ys_rx(:,n) = Ys_rx(:,n)./exp(j*f_drift(n));
   X_est  = [X_est;Ys_rx(:,n)./(H_est)];
end


X_adjust = sign(real(X_est));

errors = 0;

for k = 1:msg_block_num
    for n = 1:64
        if(real(X_adjust(n+64*(k-1))) ~= message(n+64*(k-1)))
            errors = errors+1;
        end
    end
end

percentError = errors/(length(X_adjust))
plot(message);
hold on;
plot(X_adjust, 'o');

legend('X', 'Xhat');