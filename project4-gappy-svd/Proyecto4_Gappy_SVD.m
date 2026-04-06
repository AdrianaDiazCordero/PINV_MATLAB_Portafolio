f = @(x,y,z) x.^2.*(sin(5*pi*y+3*log(x.^3+y.^2+z+pi^2))-1).^2-4*x.^2.*y.^3.*(1-z).^(3/2)+(x+z-1).*(2*y-z).*cos(30*(x+z)).*log(6+x.^2.*y.^2+z.^3);
N = 100;
x = linspace(0, 1, N); % Define the range and number of points for x
y = linspace(0, 1, N); % Set a value for y
z = [0,0.5,1]; % Set a value for z
[X_mesh, Y_mesh] = meshgrid(x, y);

z_full = linspace(0,1,N);
[x_3d, y_3d, z_3d] = meshgrid(x,y,z_full);
vol_f = f(x_3d,y_3d,z_3d);
f_min_global = min(vol_f, [], 'all'); 
f_max_global = max(vol_f, [], 'all');

% Limpiamos memoria borrando las matrices 3D gigantes que ya no necesitamos
clear x_3d y_3d z_3d vol_f z_full;

z_tilde = f(X_mesh, Y_mesh, z(2));
surf(X_mesh, Y_mesh, z_tilde);
shading interp; % Para suavizar los colores
view(3);
title(['Funcion original para z=',num2str(z(2))]);

% normalizando la funcion
F = @(x,y,z) 2*((f(x,y,z) - f_min_global)./(f_max_global-f_min_global))-1;
x = linspace(0, 1, N); % Define the range and number of points for x
y = linspace(0, 1, N); % Set a value for y
[X_mesh, Y_mesh] = meshgrid(x, y);

A1 = F(X_mesh, Y_mesh, z(1));
A2 = F(X_mesh, Y_mesh, z(2));
A3 = F(X_mesh, Y_mesh, z(3));
surf(X_mesh, Y_mesh, A2);
shading interp; % Para suavizar los colores
view(3);
title(['Funcion con valores normalizados para z=',num2str(z(2))]);

porcent_remove = 0.40; %rehacer para 0.15 despues al hacer el informe

A1_gappy = A1;
A2_gappy = A2;
A3_gappy = A3;
n = size(A2,1);

elementos_totales = numel(A2_gappy); 
datos_a_borrar = round(porcent_remove*elementos_totales);

% Saca posiciones al azar
posiciones_al_azar = randperm(elementos_totales, datos_a_borrar);

% Borramos los datos que estén en esas posiciones
A2_gappy(posiciones_al_azar) = NaN;
%A1_gappy(posiciones_al_azar) = NaN;

figure;
surf(X_mesh, Y_mesh, A2_gappy);
shading interp;
view(3);
title(['Funcion gappy, con ',num2str(porcent_remove*100),' porciento de datos removidos']);

% Inicializando con otra matriz usando el fillmissing
A2_recons3 = fillmissing(A2_gappy,"makima",2,"EndValues","nearest");
figure;
surf(X_mesh, Y_mesh, A2_recons3);
shading interp;
view(3);
title('Reconstruccion usando metodo makima');

MaxE_noIter2 = max(abs(A2-A2_recons3),[],"all");
fprintf('El error maximo para la reconstruccion solo con makima: %.10f\n', MaxE_noIter2);
RMSE_noIter2 = (1/n)*norm(A2-A2_recons3);
fprintf('El RMSE para la versión reconstruccion solo con makima es: %.10f\n', RMSE_noIter2);



% Reconstruir la matriz gappy
A2_gappy(posiciones_al_azar) = 0;
A2_recons = A2_gappy;
r = rank(A2_recons);

MaxE_noIter0=max(abs(A2-A2_gappy),[],"all");
RMSE_Iter0 = (1/n)*norm(A2-A2_gappy);

[U_rec,S_rec,V_rec] = svd(A2_recons,"econ");
singularValues = diag(S_rec);
figure;
plot(1:length(singularValues), singularValues);
grid on;
xlabel('Index');
ylabel('Singular Values');
title('Singular Values of Matrix A2 Inicial');


m_op = sum(singularValues > 1e-10);
fprintf('La cantidad de modos optimos es: %.10f\n', m_op);

A2_recons2 = A2_gappy;
m_op = 20;
Vector_vals_gappy = zeros(size(A2_recons2));
for i = 1:m_op
     Ei = U_rec(:, i) * V_rec(:, i)';
     Vector_vals_gappy = Vector_vals_gappy + singularValues(i)*Ei;
     A2_recons2(posiciones_al_azar) = Vector_vals_gappy(posiciones_al_azar);
end
figure;
surf(X_mesh, Y_mesh, A2_recons2);
shading interp;
view(3);
title(['Funcion reconstruida para ',num2str(m_op),' modos']);


s = 20; %numero de iteraciones
A2_recons = A2_gappy;
m = 14; %numero de modos
for k = 1:s
    [U,S,V] = svd(A2_recons,"econ");
    U_m = U(:,1:m);
    S_m = S(1:m,1:m);
    V_m = V(:,1:m);
    Ak = U_m*S_m*V_m';
    A2_recons(posiciones_al_azar) = Ak(posiciones_al_azar);
end
figure;
surf(X_mesh, Y_mesh, A2_recons);
shading interp;
view(3);
title(['Reconstruccion de inicializacion con ceros y, ',num2str(m),' modos y ', num2str(s),' iteraciones']);

% Calculo de errores
MaxE_noIter = max(abs(A2-A2_recons2),[],"all");
fprintf('El error maximo para la versión sin iterar es: %.10f\n', MaxE_noIter);
MaxE_Iter = max(abs(A2-A2_recons), [],"all");
fprintf('El error maximo para la versión iterada es: %.10f\n', MaxE_Iter);
RMSE_noIter = (1/n)*norm(A2-A2_recons2);
fprintf('El RMSE para la versión sin iterar es: %.10f\n', RMSE_noIter);
RMSE_Iter = (1/n)*norm(A2-A2_recons);
fprintf('El RMSE para la versión iterada es: %.10f\n', RMSE_Iter);

for k = 1:s
    [U,S,V] = svd(A2_recons3,"econ");
    U_m = U(:,1:m);
    S_m = S(1:m,1:m);
    V_m = V(:,1:m);
    Ak = U_m*S_m*V_m';
    A2_recons3(posiciones_al_azar) = Ak(posiciones_al_azar);
end
figure;
surf(X_mesh, Y_mesh, A2_recons3);
shading interp;
view(3);
title(['Reconstrucion inicializada makima con, ',num2str(m),' modos y ', num2str(s),' iteraciones']);

MaxE_Iter3 = max(abs(A2-A2_recons3), [],"all");
fprintf('El error maximo para la versión inicializacion makima iterada es: %.10f\n', MaxE_Iter3);
RMSE_Iter3 = (1/n)*norm(A2-A2_recons3);
fprintf('El RMSE para la versión makima iterada es: %.10f\n', RMSE_Iter3);



img_orig = imread("mansion.jpg");
if size(img_orig,3) == 3
    img_2d = rgb2gray(img_orig);
else
    img_2d = img_orig;
end

A_img = double(img_2d);

figure;
imagesc(A_img);
colormap gray;
axis image;
title('Imagen original');


porcent_remove = 0.40; %rehacer para 0.15 despues al hacer el informe

A_img_gappy = A_img;

elementos_totales = numel(A_img_gappy); 
datos_a_borrar = round(porcent_remove*elementos_totales);

% Saca posiciones al azar
posiciones_al_azar = randperm(elementos_totales, datos_a_borrar);

% Borramos los datos que estén en esas posiciones
A_img_gappy(posiciones_al_azar) = NaN;
%A1_gappy(posiciones_al_azar) = NaN;

figure;
imagesc(A_img_gappy);
colormap gray;
axis image;
title(['Imagen con ',num2str(porcent_remove*100),' porciento de datos removidos']);

A_img_recon = fillmissing(A_img_gappy,"makima",2,"EndValues","nearest");
figure;
imagesc(A_img_recon);
colormap gray;
axis image;
title('Inicializacion usando metodo makima');

tamano = size(A_img_recon,1);

MaxE_Ini = max(abs(A_img-A_img_recon),[],"all");
fprintf('Error maximo para inicializacion con makima: %.10f\n', MaxE_Ini);
RMSE_Ini = (1/tamano)*norm(A_img-A_img_recon);
fprintf('El RMSE para inicializacion con makima es: %.10f\n', RMSE_Ini);

[~,S,~] = svd(A_img_recon,"econ");
singularValues = diag(S);
m_op = sum(singularValues > 1e-10);
fprintf('La cantidad de modos optimos es: %.10f\n', m_op);

s = 50; %numero de iteraciones
m = floor(0.1*m_op); %numero de modos
for k = 1:s
    [U,S,V] = svd(A_img_recon,"econ");
    U_m = U(:,1:m);
    S_m = S(1:m,1:m);
    V_m = V(:,1:m);
    Ak_img = U_m*S_m*V_m';
    A_img_recon(posiciones_al_azar) = Ak_img(posiciones_al_azar);
end
figure;
imagesc(A_img_recon);
colormap gray;
axis image;
title(['Reconstruccion, inicializacion makima y, ',num2str(m),' modos y ', num2str(s),' iteraciones']);

% Calculo de errores

MaxE_Iter_img = max(abs(A_img-A_img_recon), [],"all");
fprintf('El error maximo para la versión iterada es: %.10f\n', MaxE_Iter_img);

RMSE_Iter_img = (1/tamano)*norm(A_img-A_img_recon);
fprintf('El RMSE para la versión iterada es: %.10f\n', RMSE_Iter_img);


A_img_recon2 = fillmissing(A_img_gappy,'constant',1);
figure;
imagesc(A_img_recon2);
colormap gray;
axis image;
title('Inicializacion usando blancos');

%tamano = size(A_img_recon2,1);

MaxE_Ini = max(abs(A_img-A_img_recon2),[],"all");
fprintf('Error maximo para inicializacion con blanco: %.10f\n', MaxE_Ini);
RMSE_Ini = (1/tamano)*norm(A_img-A_img_recon2);
fprintf('El RMSE para inicializacion con blanco es: %.10f\n', RMSE_Ini);

[~,S2,~] = svd(A_img_recon2,"econ");
singularValues2 = diag(S2);
m_op = sum(singularValues2 > 1e-10);
fprintf('La cantidad de modos optimos es: %.10f\n', m_op);

s = 50; %numero de iteraciones
m = floor(0.1*m_op); %numero de modos
for k = 1:s
    [U,S,V] = svd(A_img_recon2,"econ");
    U_m = U(:,1:m);
    S_m = S(1:m,1:m);
    V_m = V(:,1:m);
    Ak_img = U_m*S_m*V_m';
    A_img_recon2(posiciones_al_azar) = Ak_img(posiciones_al_azar);
end
figure;
imagesc(A_img_recon2);
colormap gray;
axis image;
title(['Reconstruccion, ini con blancos y, ',num2str(m),' modos y ', num2str(s),' iteraciones']);

% Calculo de errores

MaxE_Iter_img = max(abs(A_img-A_img_recon2), [],"all");
fprintf('El error maximo para la versión iterada es: %.10f\n', MaxE_Iter_img);

RMSE_Iter_img = (1/tamano)*norm(A_img-A_img_recon2);
fprintf('El RMSE para la versión iterada es: %.10f\n', RMSE_Iter_img);


img_orig2 = imread("palm.jpg");
if size(img_orig2,3) == 3
    disp('--');
    img_2d_2 = rgb2gray(img_orig2);
else
    img_2d_2 = img_orig2;
end

A_img2 = double(img_2d_2);

figure;
imagesc(A_img2);
colormap gray;
axis image;
title('Imagen original');

porcent_remove = 0.40; %rehacer para 0.15 despues al hacer el informe

A_img_gappy2 = A_img2;

elementos_totales = numel(A_img_gappy2); 
datos_a_borrar = round(porcent_remove*elementos_totales);

% Saca posiciones al azar
posiciones_al_azar = randperm(elementos_totales, datos_a_borrar);

% Borramos los datos que estén en esas posiciones
A_img_gappy2(posiciones_al_azar) = NaN;
%A1_gappy(posiciones_al_azar) = NaN;

figure;
imagesc(A_img_gappy2);
colormap gray;
axis image;
title(['Imagen con ',num2str(porcent_remove*100),' porciento de datos removidos']);

A_img_recon3 = fillmissing(A_img_gappy2,'constant',1);
figure;
imagesc(A_img_recon3);
colormap gray;
axis image;
title('Inicializacion usando blancos');

tamano = size(A_img_recon3,1);

MaxE_Ini3 = max(abs(A_img2-A_img_recon3),[],"all");
fprintf('Error maximo para inicializacion con blanco: %.10f\n', MaxE_Ini3);
RMSE_Ini3 = (1/tamano)*norm(A_img2-A_img_recon3);
fprintf('El RMSE para inicializacion con blanco es: %.10f\n', RMSE_Ini3);

[~,S3,~] = svd(A_img_recon3,"econ");
singularValues3 = diag(S3);
m_op = sum(singularValues3 > 1e-10);
fprintf('La cantidad de modos optimos es: %.10f\n', m_op);

s = 50; %numero de iteraciones
m = floor(0.1*m_op); %numero de modos
for k = 1:s
    [U,S,V] = svd(A_img_recon3,"econ");
    U_m = U(:,1:m);
    S_m = S(1:m,1:m);
    V_m = V(:,1:m);
    Ak_img = U_m*S_m*V_m';
    A_img_recon3(posiciones_al_azar) = Ak_img(posiciones_al_azar);
end
figure;
imagesc(A_img_recon3);
colormap gray;
axis image;
title(['Reconstruccion, ini con blancos y, ',num2str(m),' modos y ', num2str(s),' iteraciones']);

% Calculo de errores

MaxE_Iter_img = max(abs(A_img2-A_img_recon3), [],"all");
fprintf('El error maximo para la versión iterada es: %.10f\n', MaxE_Iter_img);

RMSE_Iter_img = (1/tamano)*norm(A_img2-A_img_recon3);
fprintf('El RMSE para la versión iterada es: %.10f\n', RMSE_Iter_img);