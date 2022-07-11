%Developed and tested by Aditya Natu%

function [r] = Pathline(U,V,x1,x2,y1,y2,t1,t2,n)

linspacep = @(a1,a2,N) (N==1)*((a1+a2)/2) +(N~=1)*linspace(a1,a2,N);

if t1==t2
        
x = linspace(x1,x2,101);
X = repmat(x,n,1);
dx = (x2-x1)/100;

Y = zeros(n,101);
Y(:,1) = linspacep(y1,y2,n)';

    t=t1;
for i = 2:1:101
    
    k1 = dx*V(X(:,i-1),Y(:,i-1),t*ones(n,1))./U(X(:,i-1),Y(:,i-1),t*ones(n,1));
    k2 = dx*V(X(:,i-1)+dx/2,Y(:,i-1)+k1/2,t*ones(n,1))./U(X(:,i-1)+dx/2,Y(:,i-1)+k1/2,t*ones(n,1));
    k3 = dx*V(X(:,i-1)+dx/2,Y(:,i-1)+k2/2,t*ones(n,1))./U(X(:,i-1)+dx/2,Y(:,i-1)+k2/2,t*ones(n,1));
    k4 = dx*V(X(:,i-1)+dx,Y(:,i-1)+k2,t*ones(n,1))./U(X(:,i-1)+dx,Y(:,i-1)+k2,t*ones(n,1));
    
    Y(:,i) = Y(:,i-1) +k1/6 + k2/3 + k3/3 + k4/6 ;    
end

[X0,Y0]=meshgrid(linspace(x1,x2,51),linspace(min(min(Y)),max(max(Y)),51));
quiver(X0,Y0,U(X0,Y0,t*ones(size(X0))),V(X0,Y0,t*ones(size(X0))))
hold on
plot(x,Y')
hold off
xlabel('x')
ylabel('y')
title('Pathlines')
drawnow
    
    else
l = 10*(t2-t1);

X = zeros(n,l+1);
Y = zeros(n,l+1);

t = linspace(t1,t2,l+1);
dt = 0.1;

X(:,1) = x1*ones(n,1);
Y(:,1) = (linspacep(y1,y2,n))';

for i = 2:1:l+1
    
    k1 = dt*U(X(:,i-1),Y(:,i-1),t(i-1)*ones(n,1));
    l1 = dt*V(X(:,i-1),Y(:,i-1),t(i-1)*ones(n,1));
    
    k2 = dt*U(X(:,i-1)+k1/2,Y(:,i-1)+l1/2,(t(i-1)+dt/2)*ones(n,1));
    l2 = dt*V(X(:,i-1)+k1/2,Y(:,i-1)+l1/2,(t(i-1)+dt/2)*ones(n,1));
    
    k3 = dt*U(X(:,i-1)+k2/2,Y(:,i-1)+l2/2,(t(i-1)+dt/2)*ones(n,1));
    l3 = dt*V(X(:,i-1)+k2/2,Y(:,i-1)+l2/2,(t(i-1)+dt/2)*ones(n,1));
    
    k4 = dt*U(X(:,i-1)+k3,Y(:,i-1)+l3,(t(i-1)+dt)*ones(n,1));
    l4 = dt*V(X(:,i-1)+k3,Y(:,i-1)+l3,(t(i-1)+dt)*ones(n,1));
    
    X(:,i) = X(:, i-1) + k1/6 + k2/3 + k3/3 + k4/6;
    Y(:,i) = Y(:, i-1) + l1/6 + l2/3 + l3/3 + l4/6;
    
    xplot = X(1:n,1:i);
    yplot = Y(1:n,1:i);
    
    [X0,Y0]=meshgrid(linspace(min(min(xplot)),max(max(xplot)),51),linspace(min(min(yplot)),max(max(yplot)),51));
    quiver(X0,Y0,U(X0,Y0,(t(i))*ones(size(X0))),V(X0,Y0,(t(i))*ones(size(X0))))
    hold on
    plot(xplot',yplot')
    hold off
    xlabel('x')
    ylabel('y')
    title('Pathlines at t = '+string(t(i)))
    drawnow
    
end

r='Plotting complete';

end
end
