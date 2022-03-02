theta = 0:0.01:2*pi;
figure(1);

for c=[0.15, 0.35, 0.5, 0.75, 1.0, 1.2, 1.5]
    polarplot(theta,ftbsa(c,theta));
    hold on;
end
title('FTBS Scheme - Amplitude Factor');
legend('cfl=0.15', 'cfl=0.35', 'cfl=0.5', 'cfl=0.75', 'cfl=1.0', 'cfl=1.2', 'cfl=1.5')
hold off;
figure(2);

for c=[0.15, 0.35, 0.5, 0.75, 1.0, 1.2, 1.5]
    polarplot(theta,explfa(c,theta));
    hold on;
end
title('Explicit Lax-Friedrichs Scheme - Amplitude Factor')
legend('cfl=0.15', 'cfl=0.35', 'cfl=0.5', 'cfl=0.75', 'cfl=1.0', 'cfl=1.2', 'cfl=1.5')
hold off;
figure(3);

for c=[0.15, 0.35, 0.5, 0.75, 1.0, 1.2, 1.5]
    polarplot(theta,implfa(c,theta));
    hold on;
end
title('Implicit Lax-Friedrichs Scheme - Amplitude Factor')
legend('cfl=0.15', 'cfl=0.35', 'cfl=0.5', 'cfl=0.75', 'cfl=1.0', 'cfl=1.2', 'cfl=1.5')
hold off;
figure(4);

for c=[0.15, 0.35, 0.5, 0.75, 1.0, 1.2, 1.5]
    polarplot(theta,explfp(c,theta));
    hold on;
end
title('Explicit Lax-Friedrichs Scheme - Phase Factor')
legend('cfl=0.15', 'cfl=0.35', 'cfl=0.5', 'cfl=0.75', 'cfl=1.0', 'cfl=1.2', 'cfl=1.5')
hold off;
figure(5);

for c=[0.15, 0.35, 0.5, 0.75, 1.0, 1.2, 1.5]
    polarplot(theta,implfp(c,theta));
    hold on;
end
title('Implicit Lax-Friedrichs Scheme - Phase Factor')
legend('cfl=0.15', 'cfl=0.35', 'cfl=0.5', 'cfl=0.75', 'cfl=1.0', 'cfl=1.2', 'cfl=1.5')
hold off;


function val=ftbsa(c,theta)
         val = atan(-c*sin(theta)./(1-c+c*cos(theta)))./(-theta*c);
end

function val = explfa(c, theta)
        val = 1 + sin(theta).^2 * (c^2 - 1);
end

function val = implfa(c, theta)
        val = cos(theta).^2./(1 + c*sin(theta).^2);
end

function val = explfp(c, theta)
        val = atan(-c*tan(theta))./(-theta*c);
end

function val = implfp(c, theta)
        val = atan(-c*sin(theta))./(-theta*c)
end


    