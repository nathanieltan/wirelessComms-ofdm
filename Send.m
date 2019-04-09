%%Kimberly Winter                       4/7/19
%Send

%Generate message to send through channel using generateMessage
bufferSize=50;
headerFin=[];

for i=1:100
    headerFin=[headerFin; header];
end

preamble=timingHeader;

preamble=[preamble; preamble; preamble]./2;
mess2send=generateMessage(preamble, headerFin,message,bufferSize);

write_usrp_data_file(mess2send);