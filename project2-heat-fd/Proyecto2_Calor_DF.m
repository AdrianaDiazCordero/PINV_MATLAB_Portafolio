%PROYECTO 2

% EJERCICIO 1
D = 1;   %coeficiente de difusion
L = pi;
M = 40;  % tambien sera 20, 40
N = 400;  % tambien sera 25, 50 y 100
T = 0.1; % tambien en 1
%N = 1000*T+1;
%discretizacion espacial
delta_x = L/(M + 1);
x = linspace(0, L, M+2);
%discretizacion del tiempo
delta_t = T/(N+1);
t = linspace(0,T,N+1);
d = D*delta_t/(delta_x)^2;
Aprin = sparse(1:M, 1:M, (1-2*d), M, M);
Asup = sparse(1:M-1, 2:M, d, M, M); 
Ainf = sparse(2:M, 1:M-1, d, M, M);
A = Aprin + Asup + Ainf;
%disp(full(A))
f= @(x) 10.*sin(2*x);
%f0 =2/pi*x.*(x<=pi/2)+2/pi*(pi-x).*(x>pi/2);
u = zeros(M,N+2);
u(:,1) = f(x(2:M+1));
%u(:,1) = f0(2:M+1);
for j = 2:N+2 
    u(:,j) = A * u(:,j-1);
end
disc_x =x(2:end-1)';
tempfinal = u(:,end);
datos = [disc_x,tempfinal];

%save('udif_fin.mat',u(:,end))
%writematrix(u(:,end),'u_dif_finT1.txt')
%writematrix(datos,'datosdif_fin_p2T1.txt','Delimiter','tab')
%save("datos_dif_finitas.mat","disc_x","tempfinal")

% solucion exacta
uexact =@(x,t) 10.*exp(-4*t).*sin(2*x);
%calculo del error relativo
sol_exacta_final = uexact(x(2:M+1),T);
sol_num_final = u(:,end)';
errorRelativo = max(abs(sol_exacta_final-sol_num_final))/max(abs(sol_exacta_final));
% Display the relative error
fprintf('El error relativo es: %.10f\n', errorRelativo);
fprintf('Delta x es: %.4f\n', delta_x);
fprintf('Delta t es: %.4f\n', delta_t);
fprintf('Para un valor de d= %.4f\n', d);


% EJERCICIO 2.1
% Problema inverso 
K = 10;
T_0 = T - K*delta_t;
u_inv = u;
Ainv = inv(A);
for j =N+2:-1:(N+2 - K)
    u_inv(:,j-1) = Ainv * u_inv(:,j);
end

x_internas = x(2:M+1);          % Tamaño M (filas de u_inv)
t_plot = linspace(0, T, N+2);   % Tamaño N+2 (columnas de u_inv)
[T_mesh, X_mesh] = meshgrid(t_plot, x_internas);

figure;
subplot(1, 3, 1);
surf(X_mesh, T_mesh, u_inv);
% Display the inverse solution
xlabel('x');
ylabel('t');
zlabel('u_{inv}');
title('Solucion del Problema Inverso');
colorbar;
shading interp; % Para suavizar los colores
colorbar;
view(3);
clim([-10 10]);

u_exacta_plot = uexact(X_mesh, T_mesh);
subplot(1, 3, 2);
surf(X_mesh, T_mesh, u_exacta_plot);
% Display the exact solution
xlabel('x');
ylabel('t');
zlabel('u_{exact}');
title('Solucion Exacta');
colorbar;
shading interp; % Para suavizar los colores
view(3);
clim([-10 10]);

subplot(1, 3, 3);
plot(x_internas, u_inv(:, end-K), 'r.-', 'LineWidth', 1.5); 
hold on;
plot(x_internas, sol_exacta_final, 'b.--');
%legend('Recuperada (Malla 35)', 'Original (Malla 40 interpolada)');
title(['Comparación (K=' num2str(K) ')']);
grid on;

% EJERCICIO 2.2
% Problema inverso con ruido
epsilon = 0.0001;
% Add noise to the numerical solution
ruido = epsilon *(2*randn(size(u_inv)) - 1);
u_con_ruido = u(:,end) + ruido;
u_inv_ruido = u_con_ruido;
for j =N+2:-1:(N+2 - K)
    u_inv_ruido(:,j-1) = Ainv * u_inv_ruido(:,j);
end

figure;
subplot(1, 3, 1);
surf(X_mesh, T_mesh, u_inv_ruido);
% Display the inverse solution
xlabel('x');
ylabel('t');
zlabel('u_{inv} ruido');
title('Solucion del Problema Inverso con ruido');
colorbar;
shading interp; % Para suavizar los colores
colorbar;
view(3);
clim([-10 10]);

subplot(1, 3, 2);
surf(X_mesh, T_mesh, u_exacta_plot);
% Display the exact solution
xlabel('x');
ylabel('t');
zlabel('u_{exact}');
title('Solucion Exacta');
colorbar;
shading interp; % Para suavizar los colores
view(3);
clim([-10 10]);

subplot(1, 3, 3);
plot(x_internas, u_inv_ruido(:, end-K), 'r.-', 'LineWidth', 1.5); 
hold on;
plot(x_internas, sol_exacta_final, 'b.--');
%legend('Recuperada (Malla 35)', 'Original (Malla 40 interpolada)');
title(['Comparación (K=' num2str(K) ')']);
grid on;

% EJERCICIO 2.3
% Problema inverso con el valor final exacto

u_final_exacto = uexact(x(2:M+1), T);
u_inv_exacta = zeros(M, N+2);
u_inv_exacta(:, N+2) = u_final_exacto;
for j =N+2:-1:(N+2 - K)
    u_inv_exacta(:,j-1) = Ainv * u_inv_exacta(:,j);
end

figure;
subplot(1, 3, 1);
surf(X_mesh, T_mesh, u_inv_exacta);
% Display the inverse solution
xlabel('x');
ylabel('t');
zlabel('u_{inv} exacta');
title('Solucion del Problema Inverso con final exacto');
colorbar;
shading interp; % Para suavizar los colores
colorbar;
view(3);
clim([-10 10]);

subplot(1, 3, 2);
surf(X_mesh, T_mesh, u_exacta_plot);
% Display the exact solution
xlabel('x');
ylabel('t');
zlabel('u_{exact}');
title('Solucion Exacta');
colorbar;
shading interp; % Para suavizar los colores
view(3);
clim([-10 10]);

subplot(1, 3, 3);
plot(x_internas, u_inv_exacta(:, end-K), 'r.-', 'LineWidth', 1.5); 
hold on;
plot(x_internas, sol_exacta_final, 'b.--');
%legend('Recuperada (Malla 35)', 'Original (Malla 40 interpolada)');
title(['Comparación (K=' num2str(K) ')']);
grid on;

% EJERCICIO 3
A_tilde = sparse(1:M, 1:M, (1+2*d), M, M) +sparse(1:M-1, 2:M, -d, M, M) + sparse(2:M, 1:M-1, -d, M, M);

u_inv2 = u;
for j =N+2:-1:(N+2 - K)
    u_inv2(:,j-1) = A_tilde * u_inv2(:,j);
end

figure;
subplot(1, 3, 1);
surf(X_mesh, T_mesh, u_inv2);
% Display the inverse solution
xlabel('x');
ylabel('t');
zlabel('u_{inv}');
title('Solucion del Problema Inverso Metodo 2');
colorbar;
shading interp; % Para suavizar los colores
colorbar;
view(3);
clim([-10 10]);

subplot(1, 3, 2);
surf(X_mesh, T_mesh, u_exacta_plot);
% Display the exact solution
xlabel('x');
ylabel('t');
zlabel('u_{exact}');
title('Solucion Exacta');
colorbar;
shading interp; % Para suavizar los colores
view(3);
clim([-10 10]);

subplot(1, 3, 3);
plot(x_internas, u_inv2(:, end-K), 'r.-', 'LineWidth', 1.5); 
hold on;
plot(x_internas, sol_exacta_final, 'b.--');
%legend('Recuperada (Malla 35)', 'Original (Malla 40 interpolada)');
title(['Comparación (K=' num2str(K) ')']);
grid on;

%K = 15;
u_inv_ruido2 = u_con_ruido;
for j =N+2:-1:(N+2 - K)
    u_inv_ruido2(:,j-1) = A_tilde * u_inv_ruido2(:,j);
end

figure;
subplot(1, 3, 1);
surf(X_mesh, T_mesh, u_inv_ruido2);
% Display the inverse solution
xlabel('x');
ylabel('t');
zlabel('u_{inv} ruido');
title('Solucion del Problema Inverso con ruido Metodo 2');
colorbar;
shading interp; % Para suavizar los colores
colorbar;
view(3);
clim([-10 10]);

subplot(1, 3, 2);
surf(X_mesh, T_mesh, u_exacta_plot);
% Display the exact solution
xlabel('x');
ylabel('t');
zlabel('u_{exact}');
title('Solucion Exacta');
colorbar;
shading interp; % Para suavizar los colores
view(3);
clim([-10 10]);

subplot(1, 3, 3);
plot(x_internas, u_inv_ruido2(:, end-K), 'r.-', 'LineWidth', 1.5); 
hold on;
plot(x_internas, sol_exacta_final, 'b.--');
%legend('Recuperada (Malla 35)', 'Original (Malla 40 interpolada)');
title(['Comparación (K=' num2str(K) ')']);
grid on;

u_inv_exacta2 = zeros(M, N+2);
u_inv_exacta2(:, N+2) = u_final_exacto;
for j =N+2:-1:(N+2 - K)
    u_inv_exacta2(:,j-1) = A_tilde * u_inv_exacta2(:,j);
end

figure;
subplot(1, 3, 1);
surf(X_mesh, T_mesh, u_inv_exacta2);
% Display the inverse solution
xlabel('x');
ylabel('t');
zlabel('u_{inv} exacta');
title('Solucion del Problema Inverso con final exacto Metodo 2');
colorbar;
shading interp; % Para suavizar los colores
colorbar;
view(3);
clim([-10 10]);

subplot(1, 3, 2);
surf(X_mesh, T_mesh, u_exacta_plot);
% Display the exact solution
xlabel('x');
ylabel('t');
zlabel('u_{exact}');
title('Solucion Exacta');
colorbar;
shading interp; % Para suavizar los colores
view(3);
clim([-10 10]);

subplot(1, 3, 3);
plot(x_internas, u_inv_exacta2(:, end-K), 'r.-', 'LineWidth', 1.5); 
hold on;
plot(x_internas, sol_exacta_final, 'b.--');
%legend('Recuperada (Malla 35)', 'Original (Malla 40 interpolada)');
title(['Comparación (K=' num2str(K) ')']);
grid on;