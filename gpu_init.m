
function [psi, properties, Po]=gpu_init(psi, properties, method, Po)
    psi.wave = gpuArray(psi.wave);
    Po.pot = gpuArray(Po.pot);
    if method == 0  
        psi.matrix = gpuArray(psi.matrix);
    elseif method == 1
        psi.R = gpuArray(psi.R);
        psi.K = gpuArray(psi.K);
        properties.ham = gpuArray(properties.ham);
        properties.mass = gpuArray(properties.mass);
        properties.mu = gpuArray(properties.mu);
    end
    
end