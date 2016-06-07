trail_st{1,:} = ['2016-05-30 10:54:35';'2016-06-01 10:39:40'];
trail_st{2,:} = ['2016-05-24 12:47:48';'2016-05-25 13:32:09';'2016-05-26 13:12:08';'2016-05-27 13:40:36';'2016-05-28 14:25:36';'2016-05-29 13:33:34';...
    '2016-05-30 13:48:23';'2016-05-31 13:13:46';'2016-06-01 13:38:26'];
trail_st{3,:} = ['2016-05-23 22:10:13';'2016-05-24 22:00:06';'2016-05-25 22:33:13'];

RT = [];
Alpha = [];
Beta = [];
Theta = [];
HeartRate = [];

AF3 = 1;
F7 = 2;
F3 = 3;
FC5 = 4;
T7 = 5;
P7 = 6;
O1 = 7;
O2 = 8;
P8 = 9;
T8 = 10;
FC6 = 11;
F4 = 12;
F8 = 13;
AF4 = 14;

fftlength = 128*2;
sampling_frequency = 127.94;
f=[sampling_frequency/fftlength:sampling_frequency/fftlength:sampling_frequency]; % frequency index for the spectral array

thetaIndex = find(f>=4 & f<8);
alphaIndex = find(f>=8 & f<13);
lowBetaIndex = find(f>=13 & f<16);
highBetaIndex = find(f>=16 & f<25);
gammaIndex = find(f>=25 & f<40);
totIndex = find(f>=6 & f<=40);
            

for number = 1:3
    if(number==1)
        Phase = 'Morning'; % Phase = 'Morning','Noon','Night'
    else if(number==2)
            Phase = 'Noon';
        else if number == 3
                Phase = 'Night';
            end
        end
    end
        
    eeg_list = ls([Phase,'/eeg2016*.mat']);

    % pvt_file
    directory = ['D:\Research\Fatigue Detection\Experiment\PVT\Fatigue detection\Qianyi\',Phase];
    pvt_list = ls(directory);

    %deltaIndex = find(f>=1 & f<4);

    %outdata = [];

    for trial = 1:size(eeg_list,1)
        load([Phase,'/',eeg_list(trial,:)]);
        eeg_start_time = eeg.start_time;

        pvt_filename = [directory,'\',pvt_list(trial+2,:),'\trial.mat'];
        load(pvt_filename);
        pvt_start_time = datenum(trail_st{number}(trial,:));

        
        
        datestring = datestr(datenum(eeg_start_time));

        if(datestring(4:6)=='May')
            hr_filename = ['May',num2str(str2num(datestring(1:2))),'.mat'];
        else
            if (datestring(4:6)=='Jun')
                hr_filename = ['June',num2str(str2num(datestring(1:2))),'.mat'];
            end
        end

        load(['D:\Research\Fatigue Detection\Experiment\Fitbit\Qianyi\HeartRate','\',hr_filename]);
        time_hr = datenum(HR(:,1:6));
        time_hr = (time_hr-eeg_start_time)*24*3600;
        
        disp([datestr(eeg.start_time),',',pvt_list(trial+2,:),',',hr_filename]);

        time_pvt = (data.et+39.5)/(60*60*24)+pvt_start_time; % 39.5s is the time difference between two computers
        %rt = data.rt;
        data.rt(data.sp==0)=[];
        time_pvt(data.sp==0)=[];
        time_pvt = (time_pvt-eeg_start_time)*24*3600;
        
        eeg.alpha(:,1) = eeg.alpha(:,1)-fftlength/2;
        eeg.theta(:,1) = eeg.theta(:,1)-fftlength/2;
        eeg.lowBeta(:,1) = eeg.lowBeta(:,1)-fftlength/2;
        eeg.highBeta(:,1) = eeg.highBeta(:,1)-fftlength/2;

        time_eeg = [];

        %time_eeg_freq = eeg_start_time+(eeg.alpha(:,1)-fftlength/2)/(128*60*60*24);
        time_eeg(1) = 0;
        for i=2:eeg.pnts
            if eeg.data(1,i)-eeg.data(1,i-1)==1 || eeg.data(1,i)-eeg.data(1,i-1)==-128
                time_eeg(i) = time_eeg(i-1)+1/sampling_frequency;
            else
                if(eeg.data(1,i)-eeg.data(1,i-1)>0)
                    time_eeg(i) = time_eeg(i-1)+(eeg.data(1,i)-eeg.data(1,i-1))/sampling_frequency;
                else
                    if(eeg.data(1,i)-eeg.data(1,i-1)==0)
                        time_eeg(i) = time_eeg(i-1);
                    else
                        time_eeg(i) = time_eeg(i-1)+(eeg.data(1,i)-eeg.data(1,i-1)+129)/sampling_frequency;
                    end
                end
            end
        end
        
        time_eeg_freq = time_eeg(eeg.alpha(:,1));

%         d = diff(time_eeg);
%         time_eeg(find(d==0))=[];
%         eeg.filt(find(d==0),:)=[];
%         eeg.data(:,find(d==0))=[];


        step = 128*2;
        gyro_x = eeg.data(18,:)';
        gyro_y = eeg.data(19,:)';
        STD = movingstd(gyro_x,step)+movingstd(gyro_y,step);
        Moving = find(STD>5);
    
        
        %flag = ismember(eeg.alpha(:,1)-fftlength/2,Moving);

        % index = zeros(size(time_pvt));
        % index_gyro = zeros(size(time_pvt));
        % for i=1:length(time_pvt)
        %     [~,index(i)] = min(abs(time_eeg_freq-time_pvt(i)));
        %     [~,index_gyro(i)] = min(abs(time_eeg-time_pvt(i)));
        % end

        % figure
        % plot(time_pvt,rt);
        % title('PVT reaction time');

        minute = 60;
        second = 1;
        fast10 = [];
        slow10 = [];
        meanRT = [];
        lapses = [];
        fs = [];

        alpha = [];
        theta = [];
        beta = [];

        rt = [];
        hr = [];
        
        hr_interval = find(time_hr>=time_pvt(1) & time_hr<time_pvt(end));
%         figure;
%         subplot(2,1,1)
%         plot(time_hr(hr_interval),HR(hr_interval,7));
%         txt1 = [Phase,',',num2str(trial)];
%         text(1,HR(hr_interval(1),7),txt1,'fontsize',14);
%         axis([time_hr(hr_interval(1)) time_hr(hr_interval(end)) 60 100]);
%         subplot(2,1,2);
%         plot(time_pvt,data.rt);
%         axis([time_pvt(1) time_pvt(end) 200 1000]);

        i=1;
        for time = time_pvt(1):5*minute:time_pvt(end)
            pvt_interval = find(time_pvt>=time-0.1*minute & time_pvt<time+4.9*minute);
            eeg_interval = find(time_eeg_freq>=time-0.1*minute & time_eeg_freq<time+4.9*minute);
            hr_interval = find(time_hr>=time-0.1*minute & time_hr<time+4.9*minute);
            
            Moving_index = find(ismember(eeg.alpha(eeg_interval,1),Moving));
            eeg_interval(Moving_index) = [];
            
            if(isempty(eeg_interval)|| isempty(hr_interval))
                continue;
            end
            

    %         if(data.rt(index)<100)
    %             continue;
    %         end

%             fftlength = 2^nextpow2(length(eeg_interval)); % make the window for sample length fftlength, 2 seconds in this case
%             hanning = [1:fftlength]';
%             hanning_in = 2* pi() * (hanning - (fftlength+1)/2)/(fftlength+1); %rescaled x-axis to match sample length
%             hanning = (sin(hanning_in)./hanning_in).^2; % sinc^2
%             hanning = repmat(hanning, 1, size(eeg.raw,2)); % match to number of channels

            
            
%             fftdata = padarray(eeg.filt(eeg_interval,:),[fftlength-length(eeg_interval) 0],0,'post');
%             %fftdata(1,:)=[];
% 
%             spectrum = fft(fftdata .* hanning); % apply window to HP filtered data
%             spectrum = sqrt(spectrum .* conj(spectrum));

            alpha(i,:) = mean(eeg.alpha(eeg_interval,:));
            beta(i,:) = mean(eeg.lowBeta(eeg_interval,:)+eeg.highBeta(eeg_interval,:));
            theta(i,:) = mean(eeg.theta(eeg_interval,:));


            hr(i,:) = [mean(HR(hr_interval,7)),var(HR(hr_interval,7)),max(HR(hr_interval,7)),min(HR(hr_interval,7))];

            a = sort(data.rt(pvt_interval));
            b = a(a>100);
            L = length(b);

            rt(i,1) = trial;
            rt(i,2) = mean(b);
            rt(i,3) = mean(b(round(0.9*L):end));
            rt(i,4) = mean(b(1:max(round(0.1*L),1)));
            rt(i,5) = length(eeg_interval);

            i = i+1;
        end

        RT = [RT;rt];
        Alpha = [Alpha;alpha];
        Beta = [Beta;beta];
        Theta = [Theta;theta];
        HeartRate = [HeartRate;hr];

%         figure;
%         suptitle(datestr(pvt_start_time));
%         subplot(3,2,1)
%         plot(rt);
%         %axis([0 length(fast10) 0 1])
%         title('Reaction Time');
%         subplot(3,2,2)
%         plot(sum(Theta,2)./(sum(Theta,2)+sum(Alpha,2)+sum(Beta,2)));
%         title('Theta');
%         subplot(3,2,3)
%         plot(sum(Alpha(:,[O1,O2]),2)./(sum(Theta(:,[O1,O2]),2)+sum(Alpha(:,[O1,O2]),2)+sum(Beta(:,[O1,O2]),2)));
%         title('Alpha');
%         subplot(3,2,4)
%         plot(sum(Beta,2)./(sum(Theta,2)+sum(Alpha,2)+sum(Beta,2)));
%         %axis([0 length(fast10) 0 1])
%         title('Beta');
%         subplot(3,2,5)
%         plot((sum(Theta,2)+sum(Alpha,2))./sum(Beta,2));
%         title('(Theta+Alpha/Beta)');
%         subplot(3,2,6)
%         plot(sum(Alpha,2));
%         title('Alpha');
        %axis([0 length(fast10) 0 1])
        comb1 = [AF3,AF4];
        comb2 = [O1,O2];
        comb3 = [AF3,AF4,T7,T8];
        comb4 = [AF3,AF4,T7,T8,F7,F8,F3,F4,FC5,FC6];
        comb5 = [AF3,AF4,F3,F4];

        bad_channel = [];
        for channel_no = 25:38
            if(length(find(eeg.data(channel_no,:)>=3))/length(eeg.data(channel_no,:))<0.5)
                bad_channel = [bad_channel;channel_no-24];
            end
        end
        
        bad_channel;

        comb1 = setdiff(comb1,bad_channel);
        comb2 = setdiff(comb2,bad_channel);
        comb3 = setdiff(comb3,bad_channel);
        comb4 = setdiff(comb4,bad_channel);
        comb5 = setdiff(comb5,bad_channel);

        if(number==2)
            trial = trial+2;
        else if(number==3)
                trial = trial+11;
            end
        end
        
        Result(trial,1) = corr(rt(:,2),(sum(alpha(:,comb2),2)+sum(theta(:,comb2),2))./sum(beta(:,comb2),2));
        Result(trial,2) = corr(rt(:,2),hr(:,1));
        Result(trial,3) = corr(rt(:,2),hr(:,2));
        Result(trial,4) = corr((sum(alpha(:,comb2),2)+sum(theta(:,comb2),2))./sum(beta(:,comb2),2),hr(:,1));
        Result(trial,5) = corr((sum(alpha(:,comb2),2)+sum(theta(:,comb2),2))./sum(beta(:,comb2),2),hr(:,2));
        

%         Result(trial,1) = corr(rt(:,2),(sum(alpha(:,comb1),2)+sum(theta(:,comb1),2))./sum(beta(:,comb1),2));
%         Result(trial,2) = corr(rt(:,2),(sum(alpha(:,comb2),2)+sum(theta(:,comb2),2))./sum(beta(:,comb2),2));
%         Result(trial,3) = corr(rt(:,2),(sum(alpha(:,comb3),2)+sum(theta(:,comb3),2))./sum(beta(:,comb3),2));
%         Result(trial,4) = corr(rt(:,2),(sum(alpha(:,comb4),2)+sum(theta(:,comb4),2))./sum(beta(:,comb4),2));
%         Result(trial,5) = corr(rt(:,2),(sum(alpha(:,comb5),2)+sum(theta(:,comb5),2))./sum(beta(:,comb5),2));

%         Result(trial,6) = corr(rt(:,2),sum(theta(:,comb1),2)./(sum(alpha(:,comb1),2)+sum(theta(:,comb1),2)+sum(beta(:,comb1),2)));
%         Result(trial,7) = corr(rt(:,2),sum(theta(:,comb2),2)./(sum(alpha(:,comb2),2)+sum(theta(:,comb2),2)+sum(beta(:,comb2),2)));
%         Result(trial,8) = corr(rt(:,2),sum(theta(:,comb3),2)./(sum(alpha(:,comb3),2)+sum(theta(:,comb3),2)+sum(beta(:,comb3),2)));
%         Result(trial,9) = corr(rt(:,2),sum(theta(:,comb4),2)./(sum(alpha(:,comb4),2)+sum(theta(:,comb4),2)+sum(beta(:,comb4),2)));
%         Result(trial,10) = corr(rt(:,2),sum(theta(:,comb5),2)./(sum(alpha(:,comb5),2)+sum(theta(:,comb5),2)+sum(beta(:,comb5),2)));
% 
%         Result(trial,11) = corr(rt(:,3),(sum(alpha(:,comb1),2)+sum(theta(:,comb1),2))./sum(beta(:,comb1),2));
%         Result(trial,12) = corr(rt(:,3),(sum(alpha(:,comb2),2)+sum(theta(:,comb2),2))./sum(beta(:,comb2),2));
%         Result(trial,13) = corr(rt(:,3),(sum(alpha(:,comb3),2)+sum(theta(:,comb3),2))./sum(beta(:,comb3),2));
%         Result(trial,14) = corr(rt(:,3),(sum(alpha(:,comb4),2)+sum(theta(:,comb4),2))./sum(beta(:,comb4),2));
%         Result(trial,15) = corr(rt(:,3),(sum(alpha(:,comb5),2)+sum(theta(:,comb5),2))./sum(beta(:,comb5),2));
% 
%         Result(trial,16) = corr(rt(:,3),sum(theta(:,comb1),2)./(sum(alpha(:,comb1),2)+sum(theta(:,comb1),2)+sum(beta(:,comb1),2)));
%         Result(trial,17) = corr(rt(:,3),sum(theta(:,comb2),2)./(sum(alpha(:,comb2),2)+sum(theta(:,comb2),2)+sum(beta(:,comb2),2)));
%         Result(trial,18) = corr(rt(:,3),sum(theta(:,comb3),2)./(sum(alpha(:,comb3),2)+sum(theta(:,comb3),2)+sum(beta(:,comb3),2)));
%         Result(trial,19) = corr(rt(:,3),sum(theta(:,comb4),2)./(sum(alpha(:,comb4),2)+sum(theta(:,comb4),2)+sum(beta(:,comb4),2)));
%         Result(trial,20) = corr(rt(:,3),sum(theta(:,comb5),2)./(sum(alpha(:,comb5),2)+sum(theta(:,comb5),2)+sum(beta(:,comb5),2)));
% 
%         Result(trial,21) = corr(rt(:,4),(sum(alpha(:,comb1),2)+sum(theta(:,comb1),2))./sum(beta(:,comb1),2));
%         Result(trial,22) = corr(rt(:,4),(sum(alpha(:,comb2),2)+sum(theta(:,comb2),2))./sum(beta(:,comb2),2));
%         Result(trial,23) = corr(rt(:,4),(sum(alpha(:,comb3),2)+sum(theta(:,comb3),2))./sum(beta(:,comb3),2));
%         Result(trial,24) = corr(rt(:,4),(sum(alpha(:,comb4),2)+sum(theta(:,comb4),2))./sum(beta(:,comb4),2));
%         Result(trial,25) = corr(rt(:,4),(sum(alpha(:,comb5),2)+sum(theta(:,comb5),2))./sum(beta(:,comb5),2));
% 
%         Result(trial,26) = corr(rt(:,4),sum(theta(:,comb1),2)./(sum(alpha(:,comb1),2)+sum(theta(:,comb1),2)+sum(beta(:,comb1),2)));
%         Result(trial,27) = corr(rt(:,4),sum(theta(:,comb2),2)./(sum(alpha(:,comb2),2)+sum(theta(:,comb2),2)+sum(beta(:,comb2),2)));
%         Result(trial,28) = corr(rt(:,4),sum(theta(:,comb3),2)./(sum(alpha(:,comb3),2)+sum(theta(:,comb3),2)+sum(beta(:,comb3),2)));
%         Result(trial,29) = corr(rt(:,4),sum(theta(:,comb4),2)./(sum(alpha(:,comb4),2)+sum(theta(:,comb4),2)+sum(beta(:,comb4),2)));
%         Result(trial,30) = corr(rt(:,4),sum(theta(:,comb5),2)./(sum(alpha(:,comb5),2)+sum(theta(:,comb5),2)+sum(beta(:,comb5),2)));
%         
%         Result(trial,31) = corr(rt(:,2),sum(alpha(:,comb2),2)./(sum(alpha(:,comb2),2)+sum(theta(:,comb2),2)+sum(beta(:,comb2),2)));
%         Result(trial,32) = corr(rt(:,3),sum(alpha(:,comb2),2)./(sum(alpha(:,comb2),2)+sum(theta(:,comb2),2)+sum(beta(:,comb2),2)));
%         Result(trial,33) = corr(rt(:,4),sum(alpha(:,comb2),2)./(sum(alpha(:,comb2),2)+sum(theta(:,comb2),2)+sum(beta(:,comb2),2)));
%         
%         Result(trial,34) = stats.rt_mean;
%         Result(trial,35) = stats.rt_std;
%         Result(trial,36) = stats.rt_fastest_10;
%         Result(trial,37) = stats.rt_slowest_10;
        
        
        
        x = (sum(alpha(:,comb2),2)+sum(theta(:,comb2),2))./sum(beta(:,comb2),2);
        p(trial,:) = polyfit(x,rt(:,2),1);
        
%         figure;
%         plot(x,rt(:,2),'ro');
%         hold on;
%         plot(min(x):max(x),polyval(p(trial,:),min(x):max(x)));

        if(number==2)
            trial = trial-2;
        else if(number==3)
                trial = trial-11;
            end
        end

    end
end

trial = 15;

Result(trial,1) = corr(RT(:,2),(sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2))./sum(Beta(:,comb2),2));
Result(trial,2) = corr(RT(:,2),HeartRate(:,1));
Result(trial,3) = corr(RT(:,2),HeartRate(:,2));
Result(trial,4) = corr((sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2))./sum(Beta(:,comb2),2),HeartRate(:,1));
Result(trial,5) = corr((sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2))./sum(Beta(:,comb2),2),HeartRate(:,2));
        
% Result(trial,1) = corr(RT(:,2),(sum(Alpha(:,comb1),2)+sum(Theta(:,comb1),2))./sum(Beta(:,comb1),2));
% Result(trial,2) = corr(RT(:,2),(sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2))./sum(Beta(:,comb2),2));
% Result(trial,3) = corr(RT(:,2),(sum(Alpha(:,comb3),2)+sum(Theta(:,comb3),2))./sum(Beta(:,comb3),2));
% Result(trial,4) = corr(RT(:,2),(sum(Alpha(:,comb4),2)+sum(Theta(:,comb4),2))./sum(Beta(:,comb4),2));
% Result(trial,5) = corr(RT(:,2),(sum(Alpha(:,comb5),2)+sum(Theta(:,comb5),2))./sum(Beta(:,comb5),2));
% 
% Result(trial,6) = corr(RT(:,2),sum(Theta(:,comb1),2)./(sum(Alpha(:,comb1),2)+sum(Theta(:,comb1),2)+sum(Beta(:,comb1),2)));
% Result(trial,7) = corr(RT(:,2),sum(Theta(:,comb2),2)./(sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2)+sum(Beta(:,comb2),2)));
% Result(trial,8) = corr(RT(:,2),sum(Theta(:,comb3),2)./(sum(Alpha(:,comb3),2)+sum(Theta(:,comb3),2)+sum(Beta(:,comb3),2)));
% Result(trial,9) = corr(RT(:,2),sum(Theta(:,comb4),2)./(sum(Alpha(:,comb4),2)+sum(Theta(:,comb4),2)+sum(Beta(:,comb4),2)));
% Result(trial,10) = corr(RT(:,2),sum(Theta(:,comb5),2)./(sum(Alpha(:,comb5),2)+sum(Theta(:,comb5),2)+sum(Beta(:,comb5),2)));
% 
% Result(trial,11) = corr(RT(:,3),(sum(Alpha(:,comb1),2)+sum(Theta(:,comb1),2))./sum(Beta(:,comb1),2));
% Result(trial,12) = corr(RT(:,3),(sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2))./sum(Beta(:,comb2),2));
% Result(trial,13) = corr(RT(:,3),(sum(Alpha(:,comb3),2)+sum(Theta(:,comb3),2))./sum(Beta(:,comb3),2));
% Result(trial,14) = corr(RT(:,3),(sum(Alpha(:,comb4),2)+sum(Theta(:,comb4),2))./sum(Beta(:,comb4),2));
% Result(trial,15) = corr(RT(:,3),(sum(Alpha(:,comb5),2)+sum(Theta(:,comb5),2))./sum(Beta(:,comb5),2));
% 
% Result(trial,16) = corr(RT(:,3),sum(Theta(:,comb1),2)./(sum(Alpha(:,comb1),2)+sum(Theta(:,comb1),2)+sum(Beta(:,comb1),2)));
% Result(trial,17) = corr(RT(:,3),sum(Theta(:,comb2),2)./(sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2)+sum(Beta(:,comb2),2)));
% Result(trial,18) = corr(RT(:,3),sum(Theta(:,comb3),2)./(sum(Alpha(:,comb3),2)+sum(Theta(:,comb3),2)+sum(Beta(:,comb3),2)));
% Result(trial,19) = corr(RT(:,3),sum(Theta(:,comb4),2)./(sum(Alpha(:,comb4),2)+sum(Theta(:,comb4),2)+sum(Beta(:,comb4),2)));
% Result(trial,20) = corr(RT(:,3),sum(Theta(:,comb5),2)./(sum(Alpha(:,comb5),2)+sum(Theta(:,comb5),2)+sum(Beta(:,comb5),2)));
% 
% Result(trial,21) = corr(RT(:,4),(sum(Alpha(:,comb1),2)+sum(Theta(:,comb1),2))./sum(Beta(:,comb1),2));
% Result(trial,22) = corr(RT(:,4),(sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2))./sum(Beta(:,comb2),2));
% Result(trial,23) = corr(RT(:,4),(sum(Alpha(:,comb3),2)+sum(Theta(:,comb3),2))./sum(Beta(:,comb3),2));
% Result(trial,24) = corr(RT(:,4),(sum(Alpha(:,comb4),2)+sum(Theta(:,comb4),2))./sum(Beta(:,comb4),2));
% Result(trial,25) = corr(RT(:,4),(sum(Alpha(:,comb5),2)+sum(Theta(:,comb5),2))./sum(Beta(:,comb5),2));
% 
% Result(trial,26) = corr(RT(:,4),sum(Theta(:,comb1),2)./(sum(Alpha(:,comb1),2)+sum(Theta(:,comb1),2)+sum(Beta(:,comb1),2)));
% Result(trial,27) = corr(RT(:,4),sum(Theta(:,comb2),2)./(sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2)+sum(Beta(:,comb2),2)));
% Result(trial,28) = corr(RT(:,4),sum(Theta(:,comb3),2)./(sum(Alpha(:,comb3),2)+sum(Theta(:,comb3),2)+sum(Beta(:,comb3),2)));
% Result(trial,29) = corr(RT(:,4),sum(Theta(:,comb4),2)./(sum(Alpha(:,comb4),2)+sum(Theta(:,comb4),2)+sum(Beta(:,comb4),2)));
% Result(trial,30) = corr(RT(:,4),sum(Theta(:,comb5),2)./(sum(Alpha(:,comb5),2)+sum(Theta(:,comb5),2)+sum(Beta(:,comb5),2)));
% 
% Result(trial,31) = corr(RT(:,2),sum(Alpha(:,comb2),2)./(sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2)+sum(Beta(:,comb2),2)));
% Result(trial,32) = corr(RT(:,3),sum(Alpha(:,comb2),2)./(sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2)+sum(Beta(:,comb2),2)));
% Result(trial,33) = corr(RT(:,4),sum(Alpha(:,comb2),2)./(sum(Alpha(:,comb2),2)+sum(Theta(:,comb2),2)+sum(Beta(:,comb2),2)));

%%
Sleep = [7+55/60;8+7/60;3+47/60;8+8/60;8+53/60;6+57/60;4+44/60;9+30/60;7+55/60;8+10/60;8+7/60;6+20/60;3+47/60;8+8/60;];

% for i=1:30
%     figure;
%     bar(Sleep,Result(3:11,i));
% end

find(diff(RT(:,1))>0);
start_index = [1;find(diff(RT(:,1))~=0)+1];
end_index = [find(diff(RT(:,1))~=0);size(RT,1)];
%end_index(11) = end_index(11)-1;
for index = 3:10
    %figure;
    %figure;
    x = (sum(Alpha(start_index(index):end_index(index),comb2),2)+sum(Theta(start_index(index):end_index(index),comb2),2))./sum(Beta(start_index(index):end_index(index),comb2),2);
    plot(x,RT(start_index(index):end_index(index),2),'o');
    hold on;
    plot(min(x):0.01:max(x),polyval(p(index,:),min(x):0.01:max(x)));
    %txt1 = ['Sleep: ',num2str(Sleep(index))];
    txt1 = ['HR: ',num2str(mean(HeartRate(start_index(index):end_index(index),1)))];
    if(index==7||index==4||index==10)
        text(max(x),polyval(p(index,:),max(x)),txt1,'fontsize',10);
    else
        text(min(x),polyval(p(index,:),min(x)),txt1,'fontsize',10);
    end
    %legend(['y=',num2str(p(index,1)),'*x+',num2str(p(index,2))]);
end

%%

find(diff(RT(:,1))>0);
start_index = [1;find(diff(RT(:,1))~=0)+1];
end_index = [find(diff(RT(:,1))~=0);size(RT,1)];
%end_index(11) = end_index(11)-1;
for index = 2:2
    %figure;
    x = HeartRate(start_index(index):end_index(index),2);
    plot(x,RT(start_index(index):end_index(index),2),'o');
    hold on;
    %plot(min(x):0.01:max(x),polyval(p(index,:),min(x):0.01:max(x)),'c');
    %legend(['y=',num2str(p(index,1)),'*x+',num2str(p(index,2))]);
end

%%
