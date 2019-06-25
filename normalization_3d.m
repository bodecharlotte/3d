clear all
close all
    
    % c positive/negative defocusing/focusing
    c=[0,3.137, 12.54, 31.371, 62.742, 156.855, 313.71, 627.42, 1254.8];
    dt = [0.005, 0.001, 0.0005,0.00025];
    
    xt_struct=struct('xmax', 8, 'xmin', -8,...
        'ymax', 6, 'ymin', -6, ...
        'zmax', 4, 'zmin', -4, ...
        'res', 128, 'dt',dt(2), 'tf',20);
    myParam=struct('timesteps', xt_struct.tf/xt_struct.dt,'dx', (xt_struct.xmax-xt_struct.xmin)/xt_struct.res,...
        'dy', (xt_struct.ymax-xt_struct.ymin)/xt_struct.res,...
        'dz', (xt_struct.zmax-xt_struct.zmin)/xt_struct.res,...
        'x',xt_struct.xmin + xt_struct.xmax/xt_struct.res : (xt_struct.xmax-xt_struct.xmin)/xt_struct.res : xt_struct.xmax,...
        'y', xt_struct.ymin + xt_struct.ymax/xt_struct.res : (xt_struct.ymax-xt_struct.ymin)/xt_struct.res : xt_struct.ymax,...
        'z', xt_struct.zmin + xt_struct.zmax/xt_struct.res : (xt_struct.zmax-xt_struct.zmin)/xt_struct.res : xt_struct.zmax,...
        'dk', pi/(xt_struct.xmax-xt_struct.xmin),'dw', pi/(xt_struct.ymax-xt_struct.ymin),'dv', pi/(xt_struct.zmax-xt_struct.zmin),...
        'k', cat(2, 0:xt_struct.res/2 - 1 , -xt_struct.res/2 : -1) * pi/(xt_struct.xmax-xt_struct.xmin),...
        'w', cat(2, 0:xt_struct.res/2 - 1 , -xt_struct.res/2 : -1) * pi/(xt_struct.ymax-xt_struct.ymin),...
        'v', cat(2, 0:xt_struct.res/2 - 1 , -xt_struct.res/2 : -1) * pi/(xt_struct.zmax-xt_struct.zmin),...
        'c', 200,'n_atoms', 2500 ,'imtime', true, 'nonlinear', true);
    
    %Initialize and Compute

    tic
    [opr, properties] = init(xt_struct, myParam, 0); 
    [opr, properties] = gpu_init(opr, properties);
    [opr, density, properties]=split_op(properties, xt_struct,myParam, opr);
    gpuTime= toc;
    
    %Plot
    figure(2)
    subplot(2,1,1)
    meshc(properties.X(:,:), properties.Y(:,:), density(:,:))
    xlabel('x')
    ylabel('y')
    title('Density, Z = 0')
    %title(['time t= ', num2str(t)]) 
    subplot(2,1,2)
    plot(properties.ham(1:end)) 
    title('Energy Evolution imaginary timestepping')    
    
    %Energy and chemical potential
    [E_rel] = err(myParam, properties);
    ix = find(properties.relEn < 1e-12, 1, 'first');
    if isempty(ix)==false
        z=[ix,properties.ham(ix), properties.mu(ix)]
    else
        z=[properties.ham(end), properties.mu(end)]
    end
    
    % Save to files
    fprintf('Ground state after %i timesteps\n', ix);
    dlmwrite('Box_3d_200.csv', density , 'delimiter',',','-append')
    dlmwrite('Box_3d_200.csv', xt_struct , 'delimiter',',','-append')
    dlmwrite('En_3d_200.csv', [z,xt_struct.tf, xt_struct.dt], 'delimiter',',','-append')
    dlmwrite('En_3d_200.csv', [properties.rel_diff(end), properties.relEn(end),xt_struct.tf, xt_struct.dt], 'delimiter',',','-append')
    disp('GPU time:'), disp(gpuTime);



            
function [opr, properties]=gpu_init(opr, properties)
    opr.w = gpuArray(opr.w);
    opr.V = gpuArray(opr.V);
    opr.K = gpuArray(opr.K);
    opr.R = gpuArray(opr.R);
    properties.ham = gpuArray(properties.ham);
    properties.mass = gpuArray(properties.mass);
    properties.mu = gpuArray(properties.mu);
end

function opr =  Box_pot(xt_struct, myParam)
opr. V = zeros(xt_struct.res, xt_struct.res, xt_struct.res);
    for i = 1: length(myParam.x)
        if abs(myParam.x(i)) >=2
            opr.V(i,:,:)= 1000;
        end 
    end 
    for i = 1: length(myParam.y)
        if abs(myParam.y(i)) >=2
            opr.V(:,i,:)= 1000;
        end 
    end 
    for i = 1: length(myParam.z)
        if abs(myParam.z(i)) >=2
            opr.V(:,:,i)= 1000;
        end 
    end 
end

%initialize
function [opr, properties]=init(xt_struct,myParam, voffset)
    
    w0 = 4;
    delta = 1;
    r0 = 1;
    gamma_x = 1;
    gamma_y = 1;
    gamma_z = 2;
    [X,Y,Z] = meshgrid(myParam.x, myParam.y, myParam.z); 
    [k,w,v] = meshgrid(myParam.k, myParam.w, myParam.v); 
    properties.X=X;
    properties.Y=Y;
    properties.Z=Z;    
    properties.k=k;
    properties.w=w;
    properties.v=v;
    opr =  Box_pot(xt_struct, myParam);%0.5*((gamma_x^2 * X.^2+gamma_y.^2*Y.^2 + gamma_z^2 * Z.^2));% +w0*exp(-delta*((X-r0).^2+Y.^2));
    opr.w =  (gamma_x*gamma_y*gamma_z)^(1/4)/(pi)^(3/4)*exp(-(gamma_x*X.^2+gamma_y.^2*Y.^2 + gamma_z * Z.^2)/2);
    density = real(conj(opr.w).*opr.w);
    if myParam.imtime == true
        M = myParam.dx * myParam.dy * myParam.dz * sum(density(:));
        opr.w = sqrt(1/M) * opr.w ; %1/(pi)^(1/4) * exp(-(myParam.x).^2/2);%1/(4*pi)^(1/4) * exp(-(myParam.x).^2/8);
    end
    
    %save all initial properties
   
    [M,H,mu] = conserved(myParam, opr);
    
    properties.mass(1) = M;
    properties.ham(1) = H;
    properties.mu(1) = mu;

    if (myParam.imtime)
        opr.K = exp(-0.5*(k.^2+w.^2+v.^2)*xt_struct.dt);
        if (myParam.nonlinear)
            opr.R = exp(-0.5*(opr.V + myParam.c*abs(opr.w).^2)*xt_struct.dt);
        else
            opr.R = exp(-0.5*(opr.V)*xt_struct.dt);
        end
    else
        opr.K = exp(-1i*0.5*(k.^2+w.^2+v.^2)*xt_struct.dt);
        if (myParam.nonlinear)
            opr.R = exp(-1i*0.5*(opr.V + myParam.c*abs(opr.w).^2)*xt_struct.dt);
        else 
            opr.R = exp(-1i*0.5*(opr.V)*xt_struct.dt);
        end
    end
end

function [opr, density, properties]=split_op(properties,xt_struct,myParam, opr)
    ln=2;
    count = 0;
    w_old = opr.w;
    t = 0;
    E_old = properties.ham(1);
    %while properties.diff>1e-8
    for i = 1:myParam.timesteps%:-1:1
        t=t+i*xt_struct.dt;
       
        % Half-step in real space
        opr.w = opr.w .* opr.R;

        % fft to momentum space
        opr.w = fftn(opr.w);

        % Full step in momentum space
        opr.w = opr.w .* opr.K;

        % ifft back
        opr.w = ifftn(opr.w);

        % final half-step in real space
        opr.w = opr.w .* opr.R;
        
        % density for plotting and potential
        density = real(conj(opr.w).*opr.w);

        % renormalizing for imaginary time
        if (myParam.imtime == true)
            M = myParam.dx * myParam.dy * myParam.dz * sum(density(:));
            opr.w = sqrt(1/M) * opr.w;
  
        end
        %properties.wave(ln,:,:,:)=opr.w;
        [M,H, mu] = conserved(myParam, opr);
      
        properties.mass(ln) = M;
        properties.ham(ln) = H;
        properties.mu(ln) = mu;
        
         %compute l_inf norm
        if count==0
            properties.diff = abs(w_old-opr.w);
            properties.rel_diff(ln) = max(max(max(properties.diff)));
            rele = max(max(max(abs(E_old-H))));
            properties.relEn(ln) =rele;
        else 
            if mod(i, 2^count)==0 %compute l_inf norm of wave with same difference than before
                properties.diff = abs(w_old-opr.w);
                properties.rel_diff(ln) = max(max(max(properties.diff)));
                rele = max(max(max(abs(E_old-H))));
                properties.relEn(ln) =rele;    
            end 
            if(rele < 1e-14)
                z=[properties.ham(end), properties.mu(end)]
                break
            end
        end
        %reduce timestep
        %half timestep after 2000 steps
        if (mod(i,2000)==0)
            xt_struct.dt = xt_struct.dt/2;
            opr.R = exp(-0.5*(opr.V + myParam.c*abs(opr.w).^2)*xt_struct.dt);
            opr.K = exp(-0.5*(myParam.k.^2+myParam.w.^2+ myParam.v.^2)*xt_struct.dt);
            count = count +1;
        end
        
        % Plot while computing 
        if mod(i,100000)==0
            figure(1)
            subplot(2,1,1)
            meshc(properties.X(:,:), properties.Y(:,:), density(:,:))
            xlabel('x')
            ylabel('y')
            hold off
            title(['time t= ', num2str(t)])
            drawnow
            subplot(2,1,2)
            plot(properties.ham(1:ln))
            hold off 
            title('Energy Evolution imaginary timestepping')
            drawnow     
        end
        
        w_old = opr.w;
        E_old = H;
        ln=ln+1;     
    end  
    
end


function [M, H, mu] = conserved( myParam, opr)
    %
    % use the Trapezoidal rule to compute the mass M and hamiltonian H
    %             M = \int |u|^2
    %             H = \int |u_x|^2 - K/2 |u|^4
    %
    
    arg = real(conj(opr.w).*opr.w);
    M = myParam.dx * myParam.dy * myParam.dz * sum(arg(:));
    
    
    arg1 = real(conj(opr.w).*ifftn((myParam.k.^2+ myParam.w.^2 + myParam.v.^2) .* fftn(opr.w)));
    arg2 = real(conj(opr.w).*opr.V .* opr.w);
    arg3 = real(conj(opr.w) .* (myParam.c * abs(opr.w).^2) .* opr.w);
    
    if myParam.nonlinear ==true
        H = myParam.dx*myParam.dy*myParam.dz*sum(sum(sum(0.5*arg1+arg2+0.5*arg3))); 
    else
        H = myParam.dx*myParam.dy*myParam.dz*sum(sum(sum(0.5*arg1+arg2)));
    end
    mu = myParam.dx*myParam.dy*myParam.dz*sum(sum(sum(0.5*arg1+arg2+arg3)));
end
function plot_results(xt_struct, myParam, opr, properties, density)
    %Plot surface
    fig1= figure(1);
    subplot(2,1,1)
    meshc(properties.X(:,:), properties.Z(:,:),properties.op_pulse_initial(:,:))
    xlabel('x')
    ylabel('z')
    zlabel('density')
    title('Wavefunction Density of initial wave')
    hold on 
    subplot(2,1,2)
    meshc(properties.X(:,:), properties.Z(:,:), density(:,:))
    xlabel('x')
    ylabel('z')
    zlabel('density')
    title(['Wavefunction Density after ', num2str(myParam.timesteps), ' timesteps'])
    pngFileName= sprintf('mesh32%d_%i.png', myParam.c, xt_struct.dt);
    saveas(fig1, pngFileName)
    %plot energy evolution
    fig2 = figure();
    plot(0:xt_struct.dt:xt_struct.tf, properties.ham)
    xlabel('timesteps')
    ylabel('energy')
    title(['Energyplot after ', num2str(myParam.timesteps), ' timesteps'])
    pngFileName2= sprintf('energy3d%d_%i.png', myParam.c, xt_struct.dt);
    saveas(fig2, pngFileName2)
    properties.delta_E=max(properties.ham)-min(properties.ham);
    if myParam.imtime == false
        fprintf('energy conservation error: %i\n', properties.delta_E)
    end
    
end