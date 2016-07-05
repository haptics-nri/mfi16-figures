function [swind] = proj321_OA(aX,aY,aZ)
%Takes in SubjectNum and TrialNum; returns a321

% output = f321(x,y,z)
% Synthesizes the columns of a matrix into 1 column vector while approximately preserving the power spectral density. Samplingfrequency f_s, given in Hz.

nf = 64;
no = nf/2;

x = buffer(aX(:), no); % 10-sample nonoverlapping frames of data
y = buffer(aY(:), no);
z = buffer(aZ(:), no);

% Define the sqrt of the hamming window
w=sqrt(hamming(nf+1)); w(end)=[];      % for now always use sqrt hamming window
w=w/sqrt(sum(w(1:nf/2:nf).^2));           % normalize to give overall gain of 1

% Loop over each frame of source data, to mimic the sequential
% arrival of each single frame of data:
x_left = [];            % overflow samples that do not fill the number of samples in a window
y_left = [];
z_left = [];
optX = [];               % portion of the previous sample that we want to use for overlap
optY = [];
optZ = [];
accel_last = zeros(nf,1);   % portion of the previous sample that we need to use to overlap and add.
accel_total = [];           % full time history of output

    for cnt = 1:size(x,2) % Loop over each source frame (column)
        acqX = x(:,cnt);
        acqY = y(:,cnt);
        acqZ = z(:,cnt);
        
        [X,x_left,optX] = buffer([x_left;acqX],nf,no,optX);
        [Y,y_left,optY] = buffer([y_left;acqY],nf,no,optY);
        [Z,z_left,optZ] = buffer([z_left;acqZ],nf,no,optZ);
        
        x_w = X.*w;
        y_w = Y.*w;
        z_w = Z.*w;


        sMat = [x_w,y_w,z_w];
        f_s = 10000;

        [N_w numChan] = size(sMat); % Obtain windowsize and number of channels.

        Sens_thresh = floor(f_s/2); % Sensitivity threshhold set to 400 Hz.

        N_w_odd = mod(N_w,2);
        if N_w_odd
            idx_keep = min([floor(Sens_thresh*N_w/f_s) (N_w+1)/2]);
            SmatHalftemp = fft(sMat,[],1)/N_w;
            SmatHalf = zeros((N_w+1)/2,numChan);
            SmatHalf(1:idx_keep,:) = SmatHalftemp(1:idx_keep,:);
        else
            idx_keep = min([ceil(Sens_thresh*N_w/f_s) N_w/2+1]);
            SmatHalftemp = fft(sMat,[],1)/N_w;
            SmatHalf = zeros(N_w/2+1,numChan);
            SmatHalf(1:idx_keep,:) = SmatHalftemp(1:idx_keep,:);
        end
        %STFT window but save only half since the transform is conjugately 
        %symmetric for real data (and we WANT real output data).
        %Square the absolute value of each component to get the energy of each
        %component.
        SSmatHalf = SmatHalf.*conj(SmatHalf);

        %Add together to obtain total energy of signal in window.
        S2windHalf = sum(SSmatHalf,2);

        %Take squareroot to obtain the equivalent total fourier transform of 
        %the window without phase.
        absSwindHalf = realsqrt(S2windHalf);

        %CHOICE OF PHASE.
        %--------------------------------------------
        % 
        SumforPhase = sum(SmatHalf,2);
        % phaseS =atan(imag(SumforPhase)./real(SumforPhase));
        % Nhalf = length(phaseS);
        % for k = 1:Nhalf
        %     if sin(phaseS(k))*imag(SumforPhase(k))+cos(phaseS(k))*real(SumforPhase(k)) > 0
        %     else
        %         phaseS(k) = phaseS(k)+pi;
        %     end
        % end
        phaseS = angle(SumforPhase);

        % CHOICE 2 
        % %--------------------------------------------
        % NumComp = 2^(numChan-1);
        % NumSamp = length(SmatHalf(:,1));
        % PhaseVec = zeros(NumSamp,1);
        % 
        % for k = 1:NumSamp
        %     CompVec=([SmatHalf(k,1)+SmatHalf(k,2)+SmatHalf(k,3);
        %               SmatHalf(k,1)+SmatHalf(k,2)-SmatHalf(k,3);
        %               SmatHalf(k,1)-SmatHalf(k,2)+SmatHalf(k,3);
        %              -SmatHalf(k,1)+SmatHalf(k,2)+SmatHalf(k,3)]);
        %     [MaxVec Idx] =max(abs(CompVec));
        %     phaseS(k)=angle(CompVec(Idx));
        % end
        % 

        %--------------------------------------------    


%         absSwindHalf(1)=0; % Set DC-component to zero for mechanical reasons.

        SwindHalf = absSwindHalf .* exp(1i * phaseS);

        % Make transform conjugately symmetric.
        if N_w_odd
                SwindFin = [SwindHalf; conj(flipud(SwindHalf(2:end)))];
            else
                SwindFin = [SwindHalf; conj(flipud(SwindHalf(2:end-1)))];
        end

        accel = ifft(SwindFin,'symmetric') * (N_w);
        
        accel_w = accel.*w;
        accel_out = accel_last+accel_w;
        accel_out = accel_out(1:no);
        
        accel_last = [accel_w(end-no+1:end);zeros(no,1)];
        accel_total = [accel_total;accel_out(1:no)];
    end
    swind = accel_total;
end

