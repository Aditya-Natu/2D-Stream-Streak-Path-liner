%Developed and tested by Aditya Natu%

function [r] = Streamline(U,V,x1,x2,y1,y2,t1,t2,n)

linspacep = @(a1,a2,N) (N==1)*((a1+a2)/2) +(N~=1)*linspace(a1,a2,N);

x = linspace(x1,x2,101);
X = repmat(x,n,1);
dx = (x2-x1)/100;

Y = zeros(n,101);
Y(:,1) = linspacep(y1,y2,n)';

if t1==t2    
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
title('Streamlines')
drawnow

    
else
    
for t = t1: 0.1 :t2
    
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
title('Streamlines at t = '+string(t))

drawnow

end
end
r='Plotting complete';
end
