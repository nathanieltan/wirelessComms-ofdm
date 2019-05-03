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
%%
LTS=[1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 0, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1];
STS=[0, 0, 1+j, 0,0,0,-1-j, 0,0,0,1+j,0,0,0,-1-j, 0,0,0,1+j,0,0,0,0,0,0,0,-1-j,0,0,0,1+j,0,0,0,1+j,0,0,0,1+j,0,0,0,1+j,0,0];


% Cross correlates to find start of lts
[Ryx, lags] = xcorr(beacon, [LTS]);
[mm, ii] = max(abs(Ryx));
idx = lags(ii) + 1 - 64;
lts_rx = beacon(idx:idx+length(LTS*3) - 1);

%Schmidl-Cox Algorithm
angle_sum  = 0;
for n = 1: 64
   angle_sum = angle_sum + angle(lts_rx(128+n)/lts_rx(64+n)); 
end

f_delta_hat = angle_sum/(64*64);

rx_phase_corrected = zeros(length(rx),1);
