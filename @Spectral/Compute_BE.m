function S = Compute_BE(S)
        
        S.wave = S.mat\S.wave';
        S.wave = S.wave';
        % fft to momentum space
        S.wave = fft(S.wave);

        % Full step in momentum space
        S.wave = S.K .* S.wave; 

        % ifft back
        S.wave = ifft(S.wave);
        
        % final half-step in real space
        S.wave = S.mat\S.wave';
        
        S.wave = S.wave';
        
        
end 