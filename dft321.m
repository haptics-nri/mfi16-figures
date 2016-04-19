% DFT321 on Nx3 arrays
% also returns freq. domain data
function [out, freq] = dft321(in)

    % unity-gain first-order high-pass filter
    %a = [.02-1 0];
    %b = [1 .02-1];

    A = fft(bsxfun(@minus, in, mean(in)));
    %A = fft(filtfilt(b, a, in));
    A_s = sqrt(sum(A .* conj(A), 2));
    theta = angle(sum(A,2));
    
    freq = A_s .* exp(1j * theta);
    out = abs(ifft(freq));
    freq = abs(freq);

end
