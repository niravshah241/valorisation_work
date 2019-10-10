load('c11_range')
load('condition_number')
load('pressure_error')
load('velocity_error')

figure()
semilogy(c11_range,condition_number,'-*');
title('c11 vs condition number (semilog scale)');
xlabel('c11');
ylabel('condition number');
axis tight;

figure()
semilogy(c11_range,pressure_error,'-*');
title('c11 vs Pressure L^2 error (semilog scale)');
xlabel('c11');
ylabel('Pressure L^2 error');
axis tight;

figure()
semilogy(c11_range,velocity_error,'-*');
title('c11 vs Velocity L^2 error (semilog scale)');
xlabel('c11');
ylabel('Velocity L^2 error');
axis tight;