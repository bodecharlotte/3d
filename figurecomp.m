figure(5001)
subplot(3,1,1)
plot(reshape(P.X(end/2+1,:,end/2+1), [1,P.res]), reshape(a(end/2+1,:,end/2+1),[1,P.res]))
hold on
plot(reshape(P.X(end/2+1,:,end/2+1), [1,P.res]), reshape(psi.density(end/2+1,:,end/2+1),[1,P.res]))
xlabel('x')
ylabel('\psi(\cdot, 0,0)')
legend('hard', 'smoothed')
subplot(3,1,2)

plot(reshape(P.Y(:,end/2+1,end/2+1), [1,P.res]), reshape(a(:,end/2+1,end/2+1),[1,P.res]))
hold on
plot(reshape(P.Y(:,end/2+1,end/2+1), [1,P.res]), reshape(psi.density(:,end/2+1,end/2+1),[1,P.res]))
xlabel('y')
ylabel('\psi(0,\cdot, 0)')
legend('hard', 'smoothed')

subplot(3,1,3)
plot(reshape(P.Z(end/2+1,end/2+1,:), [1,P.res]), reshape(a(end/2+1,end/2+1,:),[1,P.res]))
hold on
plot(reshape(P.Z(end/2+1,end/2+1,:), [1,P.res]), reshape(psi.density(end/2+1,end/2+1,:),[1,P.res]))
xlabel('z')
ylabel('\psi(0,0,\cdot)')
legend('hard', 'smoothed')