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
lts_rx = beacon(idx:idx+length(LTS*2.5) - 1);

plot(real(lts_rx))