F_z = importdata("CampoGravitacional1.mat");
t = F_z(1,:);
gravF = F_z(2,:)';
n = numel(t);

%para graficar
x_plot = linspace (0,1,100);


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
    fprintf('Para k = %d, el residuo con ruido es: %e\n', k, norma_res_r(i));
    
    %calculo del condicionamiento de A
    condici = cond(A);
    fprintf('Para k = %d, el número de condición es: %e\n', k, condici);
    fprintf('Los coeficientes encontrados son: \n');
    disp(coefs)
    fprintf('Los coeficientes encontrados con ruido son: \n');
    disp(coefs_ruido)
end

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
fprintf('pues presenta un error mínimo de: %d\n', min_err_sinr);

[min_err, idx] = min(norma_res_r);
fprintf('El k óptimo basado en estabilidad frente al ruido es: %d\n', ks(idx));
fprintf('pues presenta un error mínimo de: %d\n', min_err);

% Graficar la Curva L
figure('Name', 'Curva L');
loglog(norma_res_r, norma_sol, '-s', 'LineWidth', 2);
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

disp('------*******************----------')
k_codo = calc_codoCurvaL(norma_res_sinr,norma_sol_sinr);
fprintf('El k óptimo según la curva L (automático) es: %d\n', k_codo);
disp('------*******************----------')
disp('------*******************----------')
k_codo_ruido = calc_codoCurvaL(norma_res_r,norma_sol);
fprintf('El k óptimo según la curva L (automático) para datos con ruido es es: %d\n', k_codo_ruido);
disp('------*******************----------')

%-------------------------------------------------
% VALIDACION CRUZADA
%------------------------------------------------
fprintf('Iniciando validación cruzada...\n');

porcent = 0.8; %80 porciento para entrenamiento
data_train = round(porcent*numel(gravF));

%para graficar
x_plot = linspace (0,1,100);

ks = 1:12;

rep = 20;
k_opt_vc = zeros(1,rep);
k_codo_vc = zeros(1,rep);

for r=1:rep 
idx = randperm(numel(gravF));
idx_train = idx(1:data_train);
idx_test = idx(data_train+1:end);

F_train = gravF(idx_train);
F_test = gravF(idx_test);

t_train = t(idx_train);
t_test = t(idx_test);

n_train = numel(t_train);
n_test = numel(t_test);
res_validacion = zeros(size(ks));
norm_validacion = zeros(size(ks));
for i=1:length(ks)
    k = ks(i);
    A_train = zeros(n_train, k+1);
    A_test = zeros(n_test, k+1);

    for j=1:(k+1)
        for row=1:n_train
            f = @(x) x.^(j-1)./((x-t_train(row)).^2+1).^(3/2);
            A_train(row,j) = integral(f,0,1);
        end
    end
    
    coefs_train = (A_train.'*A_train)\(A_train.'*F_train);

    for j=1:(k+1)
        for row=1:n_test
            f = @(x) x.^(j-1)./((x-t_test(row)).^2+1).^(3/2);
            A_test(row,j) = integral(f,0,1);
        end
    end
        F_predi = A_test*coefs_train;
        
    res_validacion(i) = norm(F_predi - F_test);
    norm_validacion(i) = norm(coefs_train);
end


disp('-------------------')
[min_err_opt, idx_opt_val] = min(res_validacion);
k_opt_vc(r) = ks(idx_opt_val);
k_codo_vc(r) = calc_codoCurvaL(res_validacion,norm_validacion);

end
disp('----------------------')
fprintf('La mediana para k óptimo en la validación cruzada es: %d\n', median(k_opt_vc));
disp('----------------------')
fprintf('La mediana para el codo óptimo de la curva L en la validación cruzada es: %d\n', median(k_codo_vc));

%------------------------------------------------
% FIN DE VALIDACION CRUZADA
%------------------------------------------------

