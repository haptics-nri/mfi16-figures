% DFT321 on Nx3 arrays
% also returns freq. domain data
function [out, freq] = dft321(in)

    A = fft(bsxfun(@minus, in, mean(in)));
    A_s = sqrt(sum(A .* conj(A), 2));
    theta = angle(sum(A,2));
    
    freq = A_s .* exp(1j * theta);
    out = abs(ifft(freq));
    freq = abs(freq);

end
