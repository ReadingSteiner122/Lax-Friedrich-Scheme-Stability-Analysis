x = 0:2001;
figure(1)
plot(x, s_error_FTBS);
hold on;
plot(x, errorFTBS_05);
plot(x, errorFTBS810);
plot(x, errorFTBS811);
legend('FTBS - CFL=1, P=500', 'FTBS - CFL=0.5, P=500', 'FTBS - CFL=0.5, P=80');
title('FTBS - Absolute Error');
xlabel('Number of Iterations');
ylabel('Absolute Error');
hold off;

figure(2);
plot(x, s_error_FTCS);
hold on;
plot(x, errorFTCS_05);
plot(x, errorFTCS810);
plot(x, errorFTCS811);
legend('FTCS - CFL=1, P=500', 'FTCS - CFL=0.5, P=500', 'FTCS - CFL=0.5, P=80');
title('FTCS - Absolute Error');
xlabel('Number of Iterations');
ylabel('Absolute Error');
hold off;

figure(3);
plot(x, s_error_Leapfrog);
hold on;
plot(x, errorLeapfrog_05);
plot(x, errorLeapfrog810);
plot(x, errorLeapfrog811);
legend('Leapfrog - CFL=1, P=500', 'Leapfrog - CFL=0.5, P=500', 'Leapfrog - CFL=0.5, P=80');
title('Leapfrog - Absolute Error');
xlabel('Number of Iterations');
ylabel('Absolute Error');
hold off;