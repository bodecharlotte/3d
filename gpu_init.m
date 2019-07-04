
function [psi, properties, Po]=gpu_init(psi, properties, Po)
    psi.wave = gpuArray(psi.wave);
    Po.pot = gpuArray(Po.pot);
    psi.R = gpuArray(psi.R);
    psi.K = gpuArray(psi.K);
    properties.ham = gpuArray(properties.ham);
    properties.mass = gpuArray(properties.mass);
    properties.mu = gpuArray(properties.mu);  
end