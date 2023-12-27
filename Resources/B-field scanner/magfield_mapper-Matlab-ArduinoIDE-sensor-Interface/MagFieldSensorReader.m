clear all
close all
clc
s = serial('COM7', 'BaudRate', 9600);
if isvalid(s)
    fopen(s);
    try
        collectionDuration = 100; % How Many Loops will the code run
        rcd=collectionDuration+1;
        Data = [];
        startTime = now; 
        while (now - startTime) * 24*3600 <= rcd
            raw_data = fscanf(s);   
            dataPoint = str2double(raw_data);
            Data = [Data; dataPoint];
            pause(0.05);
        end
    catch exception
        disp(['An error occurred: ', exception.message]);
    end
    fclose(s);
    delete(s);
else
    disp('Failed to open the serial port. Check the COM port and Arduino connection.');
end
Data = Data';
Data = Data(2:length(Data));
limit = length(Data);
Xrange = 1:3:limit;
Yrange = 2:3:limit;
Zrange = 3:3:limit;
FinalData(1,:) = Data(Xrange);
FinalData(2,:) = Data(Yrange);
FinalData(3,:) = Data(Zrange);

%%
figure
FinalData = FinalData./10000;
hold on
plot(1:100,FinalData(1,:))
plot(1:100,FinalData(2,:))
plot(1:100,FinalData(3,:))
legend('X-Axis','Y-Axis','Z-Axis','Location','southwest')
title('Oscilating a magnet towards and away from the sensor')
xlabel('Time (Not necessarily seconds)')
ylabel('Teslas?')
grid on
hold off