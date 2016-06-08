% Phase = 'Morning','Noon','Afternoon','Early Night','Night'
for number = 1:5
    if(number==1)
        Phase = 'Morning'; 
    else if(number==2)
            Phase = 'Noon';
        else if number == 3
                Phase = 'Afternoon';
            else if number == 4
                    Phase = 'Early Night';
                else
                    Phase = 'Night';
                end
            end
        end
    end

    % eeg_file
    eeg_list = ls([Phase,'/*.edf']);

    for trial = 1:size(eeg_list,1)
        eeg_filename = eeg_list(trial,:);
        disp(eeg_filename);
        eeg = pop_biosig([Phase,'\',eeg_filename]); %pop_biosig() - import data files into EEGLAB using BIOSIG toolbox
        
        % get time vector
        datestring = eeg_filename(length(eeg_filename)-20:length(eeg_filename)-4);
        
%         if(length(eeg_filename)==31)
%             datestring = eeg_filename(length(eeg_filename)-21:length(eeg_filename)-5);
%         else
%             datestring = eeg_filename(length(eeg_filename)-20:length(eeg_filename)-4);
%         end
        datestring = ['20',datestring(7:8),'-',datestring(4:5),'-',datestring(1:2),' ',datestring(10:11),':',datestring(13:14),':',datestring(16:17)];
        eeg_start_time = datenum(datestring,'yyyy-mm-dd HH:MM:SS');

        % remove DC offset
        % IIR_TC = 256;
        % EEG_data = OUTEEG.data(3:16,:)';
        % [rows columns]= size(EEG_data);
        % AC_EEG_data = zeros(rows-1, columns);
        % back = EEG_data(1,:);
        % 
        % for r = 2 : rows
        %     back= (back*(IIR_TC-1)+EEG_data(r,:))/IIR_TC;
        %     AC_EEG_data(r-1,:) = EEG_data(r,:)- back;
        % end
        % 
        % time(1)=[];

        eeg.raw = eeg.data(3:16,:)';
        
        med = median(eeg.raw,2); % remove median of each sample
        eeg.raw = eeg.raw - repmat(med, 1, 14);
        for j=2:size(eeg.raw,1) % limit slew rate
            del = eeg.raw(j,:) - eeg.raw(j-1,:);
            del = min(del, ones(1,14)*15);
            del = max(del, -ones(1,14)*15);
            eeg.raw(j,:) = eeg.raw(j-1,:) + del;
        end

        % High pass filter
        a = 0.0078125; % HPF filter coeffs
        b = 0.9921875;
        preVal = zeros(1,14);
        eeg.filt = zeros(size(eeg.raw));
        for j=2:size(eeg.raw,1)
            preVal = a * eeg.raw(j,:) + b * preVal;
            eeg.filt(j,:) = eeg.raw(j,:) - preVal;
        end
        % end HPF
        
        % fft
        fftlength = 128*2; % make the window for sample length fftlength, 2 seconds in this case
        hanning = [1:fftlength]';
        hanning_in = 2* pi() * (hanning - (fftlength+1)/2)/(fftlength+1); %rescaled x-axis to match sample length
        hanning = (sin(hanning_in)./hanning_in).^2; % sinc^2
        hanning = repmat(hanning, 1, size(eeg.raw,2)); % match to number of channels

        f=[128/fftlength:128/fftlength:128]; % frequency index for the spectral array
        %deltaIndex = find(f>=1 & f<4);
        thetaIndex = find(f>=4 & f<8);
        alphaIndex = find(f>=8 & f<13);
        lowBetaIndex = find(f>=13 & f<16);
        highBetaIndex = find(f>=16 & f<25);
        gammaIndex = find(f>=25 & f<40);
        totIndex = find(f>=6 & f<=40);
        
        %eeg.delta = [];
        eeg.theta = [];
        eeg.alpha = [];
        eeg.lowBeta = [];
        eeg.highBeta = [];
        eeg.gamma = []; 
        eeg.tot = []; 
        eeg.totmed = []; 

        for k = fftlength:32:size(eeg.filt,1) % step through every quarter second starting at first possible sample
            spectrum = fft(eeg.filt(k-fftlength+1:k,:) .* hanning); % apply window to HP filtered data
            spectrum = sqrt(spectrum .* conj(spectrum)); % get magnitude
            %eeg.delta = [eeg.delta; k sum(spectrum(deltaIndex,:))];
            eeg.theta = [eeg.theta; k sum(spectrum(thetaIndex,:))]; % append total spectral power in band, including sample index k
            eeg.alpha = [eeg.alpha; k sum(spectrum(alphaIndex,:))];
            eeg.lowBeta = [eeg.lowBeta; k sum(spectrum(lowBetaIndex,:))];
            eeg.highBeta = [eeg.highBeta; k sum(spectrum(highBetaIndex,:))];
            eeg.gamma = [eeg.gamma; k sum(spectrum(gammaIndex,:))];
            eeg.tot = [eeg.tot; k sum(spectrum(totIndex,:))];
        end

        eeg.dp_delta = [];
        eeg.dp_theta = [];
        eeg.dp_alpha = [];
        eeg.dp_lowBeta = [];
        eeg.dp_higBeta = [];
        %eeg.dp_gamma = []; 
        gyro_x = eeg.data(18,:)';
        gyro_y = eeg.data(19,:)';
        
        for channel = 1:14
            wpt = wpdec(eeg.filt(:,channel),4,'db4');

            eeg.dp_delta(:,channel+1) = wpcoef(wpt,[4 0]); %0-4 Hz
            eeg.dp_theta(:,channel+1) = wpcoef(wpt,[4 1]); %4-8 Hz
            eeg.dp_alpha(:,channel+1) = wpcoef(wpt,[4 2]); %8-12Hz
            eeg.dp_lowBeta(:,channel+1) = wpcoef(wpt,[4 3]);%12-16Hz
            
            highBeta = wpcoef(wpt,[2 1]);
            eeg.dp_highBeta(:,channel+1) = highBeta([1:4:length(highBeta)]);%16-32Hz
            %eeg.dp_gamma = [eeg.gamma; k sum(spectrum(gammaIndex,:))];
        end
        
        eeg.dp_delta(:,1) = 8:16:length(eeg.dp_delta)*16-8;
        eeg.dp_theta(:,1) = 8:16:length(eeg.dp_delta)*16-8;
        eeg.dp_alpha(:,1) = 8:16:length(eeg.dp_delta)*16-8;
        eeg.dp_lowBeta(:,1) = 8:16:length(eeg.dp_delta)*16-8;
        eeg.dp_highBeta(:,1) = 8:16:length(eeg.dp_highBeta)*16-8;

        save([Phase,'/eeg', [datestring(1:4),datestring(6:7),datestring(9:10),datestring(12:13),datestring(15:16),datestring(18:19)],'.mat'],'eeg');
    end
end