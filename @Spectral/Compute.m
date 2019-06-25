function S = Compute(S)

        S.wave = S.wave .* S.R;

        % fft to momentum space
        S.wave = fftn(S.wave);

        % Full step in momentum space
        S.wave = S.wave .* S.K; 

        % ifft back
        S.wave = ifftn(S.wave);

        % final half-step in real space
        S.wave =  S.wave .* S.R;
        
        S.density = real(conj(S.wave).*S.wave);
        
end 