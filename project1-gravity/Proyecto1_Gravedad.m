F_z = importdata("CampoGravitacional1.mat");
t = F_z(1,:);
gravF = F_z(2,:)';
n = numel(t);

%para graficar
x_plot = linspace (0,1,100);
%rng(1);

%perturbando los datos de gravF
eps = 0.001;
ruido = (2*rand(size(gravF)) - 1) * eps;
gravF_ruido = gravF + ruido;

f1 = figure('Name', 'Comparacion de datos sin y con ruido');
ks = 1:12;
norma_res_r = zeros(size(ks));
norma_res_sinr = zeros(size(ks));
norma_sol = zeros(size(ks));
norma_sol_sinr = zeros(size(ks));

%for k=[2,4,5,7,9,10]
for i=1:length(ks)
    k = ks(i);

    A = zeros(n, k+1);

    for j=1:(k+1)
        for row=1:n
            f = @(x) x.^(j-1)./((x-t(row)).^2+1).^(3/2);
            A(row,j) = integral(f,0,1);
        end
    end
    %calculos de coeficientes sin ruido
    coefs = (A.'*A)\(A.'*gravF);
    p_inv = flip(coefs);
    Lamda_pred = polyval(p_inv,x_plot) ;
    condicionamiento = cond(A.'*A);
    fprintf('Para k = %d, el condicionamiento de la matriz es: %e\n', k, condicionamiento);
    norma_res_sinr(i) = norm(A*coefs - gravF);
    norma_sol_sinr(i) = norm(coefs);
    %calculo de coeficientes con ruido
    coefs_ruido = (A.'*A)\(A.'*gravF_ruido);
    p_invruido = flip(coefs_ruido);
    Lamda_predruido = polyval(p_invruido,x_plot) ;

    norma_res_r(i) = norm(A*coefs_ruido - gravF_ruido);
    norma_sol(i) = norm(coefs_ruido);

    % Gráfica 1: Datos sin perturbar 
    subplot(2, 1, 1); % 2 filas, 1 columna, posición 1
    plot(x_plot, Lamda_pred, 'DisplayName', sprintf('masa para k=%d',k))
    hold on
    
    % Gráfica 2: Datos perturbados
    subplot(2, 1, 2); % 2 filas, 1 columna, posición 2
    plot(x_plot, Lamda_predruido, 'DisplayName', sprintf('masa para k=%d',k))
    hold on
    
    % Calculo de otros datos
    residuo = norm(A*coefs - gravF);
    fprintf('Para k = %d, el residuo es: %e\n', k, residuo);
   % residuo_ruido = norm(A*coefs_ruido - gravF_ruido);
    fprintf('Para k = %d, el residuo con ruido es: %e\n', k, norma_res_r(i));
    
    %calculo del condicionamiento de A
    condici = cond(A);
    fprintf('Para k = %d, el número de condición es: %e\n', k, condici);
    fprintf('Los coeficientes encontrados son: \n');
    disp(coefs)
    fprintf('Los coeficientes encontrados con ruido son: \n');
    disp(coefs_ruido)
end

%plot(t,gravF)
subplot(2, 1, 1);
title('Reconstrucción: Datos Originales');
ylabel('\lambda(x)');
legend('Location', 'best');
grid on;

subplot(2, 1, 2);
title(['Reconstrucción: Datos con Ruido (\epsilon = ', num2str(eps), ')']);
ylabel('\lambda(x)');
legend('Location', 'best');
grid on;
hold off;

figure
plot(t,gravF, 'DisplayName', sprintf('valores de Fz'));
hold on 
plot(t, A*coefs, 'DisplayName', sprintf('valores aproximados'));
plot(t, A*coefs_ruido, 'DisplayName', sprintf('valores aproximados para el ruido'));
legend('Location', 'best');
hold off

disp('-------------------')
[min_err_sinr, idx_sinr] = min(norma_res_sinr);
fprintf('El k óptimo para la estabilidad sin ruido es: %d\n', ks(idx_sinr));

[min_err, idx] = min(norma_res_r);
fprintf('El k óptimo basado en estabilidad frente al ruido es: %d\n', ks(idx));

% Graficar la Curva L
figure('Name', 'Curva L');
loglog(norma_res_r, norma_sol, '-s', 'LineWidth', 2);
%loglog(norma_res_sinr, norma_sol_sinr, '-o', 'LineWidth', 2);
text(norma_res_r, norma_sol, arrayfun(@(x) num2str(x), ks, 'UniformOutput', false), 'VerticalAlignment','bottom');
xlabel('Norma del Residuo');
ylabel('Norma de la Solución');
title('Método de la Curva L (Escala Logarítmica)');
grid on;

% Graficar la Curva L datos sin ruido
figure('Name', 'Curva L (datos sin ruido)');
loglog(norma_res_sinr, norma_sol_sinr, '-o', 'LineWidth', 2);
text(norma_res_sinr, norma_sol_sinr, arrayfun(@(x) num2str(x), ks, 'UniformOutput', false), 'VerticalAlignment','bottom');
xlabel('Norma del Residuo');
ylabel('Norma de la Solución');
title('Método de la Curva L (datos sin ruido)');
grid on;


%-------------------------------------------------
% VALIDACION CRUZADA
%------------------------------------------------
fprintf('Iniciando validación cruzada...\n');
porcent = 0.8; %80 porciento para entrenamiento
data_train = round(porcent*numel(gravF));
posiciones = randperm(numel(gravF),data_train);
F_train = gravF(posiciones);
F_test = setdiff(gravF,F_train);
t_train = t(posiciones);
t_test = setdiff(t,t_train);

n_train = numel(t_train);

%para graficar
x_plot = linspace (0,1,100);

%perturbando los datos de gravF
ruido_train = (2*rand(size(F_train)) - 1) * eps;
F_train_ruido = F_train + ruido_train;

f1 = figure('Name', 'Validacion cruzada');
ks = 1:12;
res_train = zeros(size(ks));
res_train_noisy = zeros(size(ks));
norma_sol_train = zeros(size(ks));
norma_sol_train_noisy = zeros(size(ks));

for i=1:length(ks)
    k = ks(i);
    A_train = zeros(n_train, k+1);

    for j=1:(k+1)
        for row=1:n_train
            f = @(x) x.^(j-1)./((x-t(row)).^2+1).^(3/2);
            A_train(row,j) = integral(f,0,1);
        end
    end
    %calculos de coeficientes sin ruido
    coefs_train = (A_train.'*A_train)\(A_train.'*F_train);
    p_inv_train = flip(coefs_train);
    Lamda_pred_train = polyval(p_inv_train,x_plot) ;
    cond_train = cond(A_train.'*A_train);
    fprintf('Para k = %d, el condicionamiento de la matriz es: %e\n', k, cond_train);
  %actualizar a partir de aqui
    norma_res_sinr(i) = norm(A*coefs - gravF);
    norma_sol_sinr(i) = norm(coefs);
    %calculo de coeficientes con ruido
    coefs_ruido = (A.'*A)\(A.'*gravF_ruido);
    p_invruido = flip(coefs_ruido);
    Lamda_predruido = polyval(p_invruido,x_plot) ;

    norma_res_r(i) = norm(A*coefs_ruido - gravF_ruido);
    norma_sol(i) = norm(coefs_ruido);

    % Gráfica 1: Datos sin perturbar 
    subplot(2, 1, 1); % 2 filas, 1 columna, posición 1
    plot(x_plot, Lamda_pred, 'DisplayName', sprintf('masa para k=%d',k))
    hold on
    
    % Gráfica 2: Datos perturbados
    subplot(2, 1, 2); % 2 filas, 1 columna, posición 2
    plot(x_plot, Lamda_predruido, 'DisplayName', sprintf('masa para k=%d',k))
    hold on
    
    % Calculo de otros datos
    residuo = norm(A*coefs - gravF);
    fprintf('Para k = %d, el residuo es: %e\n', k, residuo);
   % residuo_ruido = norm(A*coefs_ruido - gravF_ruido);
    fprintf('Para k = %d, el residuo con ruido es: %e\n', k, norma_res_r(i));
    
    %calculo del condicionamiento de A
    condici = cond(A);
    fprintf('Para k = %d, el número de condición es: %e\n', k, condici);
    fprintf('Los coeficientes encontrados son: \n');
    disp(coefs)
    fprintf('Los coeficientes encontrados con ruido son: \n');
    disp(coefs_ruido)
end

%plot(t,gravF)
subplot(2, 1, 1);
title('Reconstrucción: Datos Originales');
ylabel('\lambda(x)');
legend('Location', 'best');
grid on;

subplot(2, 1, 2);
title(['Reconstrucción: Datos con Ruido (\epsilon = ', num2str(eps), ')']);
ylabel('\lambda(x)');
legend('Location', 'best');
grid on;
hold off;

figure
plot(t,gravF, 'DisplayName', sprintf('valores de Fz'));
hold on 
plot(t, A*coefs, 'DisplayName', sprintf('valores aproximados'));
plot(t, A*coefs_ruido, 'DisplayName', sprintf('valores aproximados para el ruido'));
legend('Location', 'best');
hold off

disp('-------------------')
[min_err_sinr, idx_sinr] = min(norma_res_sinr);
fprintf('El k óptimo para la estabilidad sin ruido es: %d\n', ks(idx_sinr));

[min_err, idx] = min(norma_res_r);
fprintf('El k óptimo basado en estabilidad frente al ruido es: %d\n', ks(idx));

% Graficar la Curva L
figure('Name', 'Curva L');
loglog(norma_res_r, norma_sol, '-s', 'LineWidth', 2);
%loglog(norma_res_sinr, norma_sol_sinr, '-o', 'LineWidth', 2);
text(norma_res_r, norma_sol, arrayfun(@(x) num2str(x), ks, 'UniformOutput', false), 'VerticalAlignment','bottom');
xlabel('Norma del Residuo');
ylabel('Norma de la Solución');
title('Método de la Curva L (Escala Logarítmica)');
grid on;

% Graficar la Curva L datos sin ruido
figure('Name', 'Curva L (datos sin ruido)');
loglog(norma_res_sinr, norma_sol_sinr, '-o', 'LineWidth', 2);
text(norma_res_sinr, norma_sol_sinr, arrayfun(@(x) num2str(x), ks, 'UniformOutput', false), 'VerticalAlignment','bottom');
xlabel('Norma del Residuo');
ylabel('Norma de la Solución');
title('Método de la Curva L (datos sin ruido)');
grid on;

%------------------------------------------------
% FIN DE VALIDACION CRUZADA
%------------------------------------------------

% EJERCICIO 2, SEGUNDO SET DE DATOS
fprintf('Ejercicio 2 \n');
F_z2 = importdata("CampoGravitacional2.mat");
t2 = F_z2(1,:);
gravF2 = F_z2(2,:)';
n2 = numel(t2);
res_sol = zeros(size(ks));
sol= zeros(size(ks));
res_sol_ruido = zeros(size(ks));
sol_ruido = zeros(size(ks));

%perturbando los datos de gravF
ruido2 = (2*rand(size(gravF2)) - 1) * eps;
gravF_ruido2 = gravF2 + ruido2;

f2 = figure('Name', 'Comparacion de datos sin y con ruido 2do set de datos');

%for k=[2,4,5,7,10]
for m=1:length(ks)
    A2 = zeros(n2, k+1);
    k=m;

    for j=1:(k+1)
        for i=1:n2
            f2 = @(x) x.^(j-1)./((x-t2(i)).^2+1).^(3/2);
            A2(i,j) = integral(f2,0,1);
        end
    end
%calculos de coeficientes sin ruido
coefs2 = (A2.'*A2)\(A2.'*gravF2);
p_inv2 = flip(coefs2);
Lamda_pred2 = polyval(p_inv2,x_plot) ;
res_sol(m) = norm(A2*coefs2-gravF2);
sol(m) = norm(coefs2);

%calculo de coeficientes con ruido
coefs_ruido2 = (A2.'*A2)\(A2.'*gravF_ruido2);
p_invruido2 = flip(coefs_ruido2);
Lamda_predruido2 = polyval(p_invruido2,x_plot) ;
res_sol_ruido(m) = norm(A2*coefs_ruido2-gravF_ruido2);
sol_ruido(m)=norm(coefs_ruido2);

% Gráfica 1: Datos sin perturbar 
subplot(2, 1, 1); % 2 filas, 1 columna, posición 1
plot(x_plot, Lamda_pred2, 'DisplayName', sprintf('masa para k=%d',k))
hold on

% Gráfica 2: Datos perturbados
subplot(2, 1, 2); % 2 filas, 1 columna, posición 2
plot(x_plot, Lamda_predruido2, 'DisplayName', sprintf('masa para k=%d',k))
hold on

% Calculo de otros datos
residuo2 = norm(A2*coefs2 - gravF2);
fprintf('Para k = %d, el residuo es: %e\n', k, residuo2);
residuo_ruido2 = norm(A2*coefs_ruido2 - gravF_ruido2);
fprintf('Para k = %d, el residuo con ruido es: %e\n', k, residuo_ruido2);

%calculo del condicionamiento de A
condici2 = cond(A2);
fprintf('Para k = %d, el número de condición es: %e\n', k, condici2);
fprintf('Los coeficientes encontrados son: \n');
disp(coefs2)
fprintf('Los coeficientes encontrados con ruido son: \n');
disp(coefs_ruido2)
end

subplot(2, 1, 1);
title('Reconstrucción: Datos Originales 2do set de datos');
ylabel('\lambda(x)');
legend('Location', 'best');
grid on;

subplot(2, 1, 2);
title(['Reconstrucción: Datos con Ruido (\epsilon = ', num2str(eps), ')']);
ylabel('\lambda(x)');
legend('Location', 'best');
grid on;
hold off;

% Graficar la Curva L
figure('Name', 'Curva L (datos con ruido) 2do set de datos');
loglog(res_sol_ruido, sol_ruido, '-s', 'LineWidth', 2);
%loglog(norma_res_sinr, norma_sol_sinr, '-o', 'LineWidth', 2);
text(res_sol_ruido, sol_ruido, arrayfun(@(x) num2str(x), ks, 'UniformOutput', false), 'VerticalAlignment','bottom');
xlabel('Norma del Residuo');
ylabel('Norma de la Solución');
title('Método de la Curva L (Escala Logarítmica)');
grid on;

% Graficar la Curva L datos sin ruido
figure('Name', 'Curva L (datos sin ruido) 2do set de datos');
loglog(res_sol, sol, '-o', 'LineWidth', 2);
text(res_sol, sol, arrayfun(@(x) num2str(x), ks, 'UniformOutput', false), 'VerticalAlignment','bottom');
xlabel('Norma del Residuo');
ylabel('Norma de la Solución');
title('Método de la Curva L (datos sin ruido)');
grid on;

% figure
% plot(t2,gravF2, 'DisplayName', sprintf('2do set de datos valores de Fz'));
% hold on 
% plot(t2, A2*coefs2, 'DisplayName', sprintf('valores aproximados'));
% legend('Location', 'best');
% hold off
