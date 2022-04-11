

figure(1);
plot(FTBSSchemet4cfl0(1:101), FTBSSchemet4cfl0(102:202));
hold on;

plot(Exacta2(1:101), Exacta2(102:202));
title('FTBS - CFL = 0.5 v/s Exact Soln at t=4');
xlabel('x')
ylabel('u')
legend('FTBS', 'Exact');
hold off;

figure(2);
plot(FTCSSchemet4cfl0(1:101), FTCSSchemet4cfl0(102:202));
hold on;

plot(Exacta2(1:101), Exacta2(102:202));
title('FTCS - CFL = 0.5 v/s Exact Soln at t = 4');
xlabel('x')
ylabel('u')
legend('FTCS', 'Exact');
hold off;

figure(3);
plot(Leapfrogt4cfl0(1:101), Leapfrogt4cfl0(102:202));
hold on;

plot(Exacta2(1:101), Exacta2(102:202));
title('LeapFrog - CFL = 0.5 v/s Exact Soln at t=4');
xlabel('x')
ylabel('u')
legend('LeapFrog', 'Exact');
hold off;



figure(4);
plot(FTBSSchemet4cfl1(1:101), FTBSSchemet4cfl1(102:202));
hold on;

plot(Exacta2(1:101), Exacta2(102:202));
title('FTBS - CFL = 1 v/s Exact Soln at t=4');
xlabel('x')
ylabel('u')
legend('FTBS', 'Exact');
hold off;

figure(5);
plot(FTCSSchemet4cfl1(1:101), FTCSSchemet4cfl1(102:202));
hold on;

plot(Exacta2(1:101), Exacta2(102:202));
title('FTCS - CFL = 1 v/s Exact Soln at t = 4');
xlabel('x')
ylabel('u')
legend('FTCS', 'Exact');
hold off;

figure(6);
plot(Leapfrogt4cfl1(1:101), Leapfrogt4cfl1(102:202));
hold on;

plot(Exacta2(1:101), Exacta2(102:202));
title('LeapFrog - CFL = 1 v/s Exact Soln at t=4');
xlabel('x')
ylabel('u')
legend('LeapFrog', 'Exact');
hold off;



figure(7);
plot(FTBSSchemet4cfl3(1:101), FTBSSchemet4cfl3(102:202));
hold on;

plot(Exacta2(1:101), Exacta2(102:202));
title('FTBS - CFL = 3 v/s Exact Soln at t=4');
xlabel('x')
ylabel('u')
legend('FTBS', 'Exact');
hold off;

figure(8);
plot(FTCSSchemet4cfl3(1:101), FTCSSchemet4cfl3(102:202));
hold on;

plot(Exacta2(1:101), Exacta2(102:202));
title('FTCS - CFL = 3 v/s Exact Soln at t = 4');
xlabel('x')
ylabel('u')
legend('FTCS', 'Exact');
hold off;

figure(9);
plot(Leapfrogt4cfl3(1:101), Leapfrogt4cfl3(102:202));
hold on;

plot(Exacta2(1:101), Exacta2(102:202));
title('LeapFrog - CFL = 3 v/s Exact Soln at t=4');
xlabel('x')
ylabel('u')
legend('LeapFrog', 'Exact');
hold off;



