function S = mat_spec(S, P, Po)
    S.mat = diag((eye(1,P.res)+P.timestep*Po.pot + P.nonlinear*P.timestep*abs(S.wold).^2), 0);
    S.mat = sparse(S.mat);
end