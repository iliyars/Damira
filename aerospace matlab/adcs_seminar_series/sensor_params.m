%%%Update Rate
nextSensorUpdate = 1.0;
fsensor=1;
%%%Bias and Noise
MagscaleBias = (4e-7)*fsensor; %%T
MagFieldBias = MagscaleBias*(2*rand()-1); %%% -1 to 1 (0 to 1)



AngscaleBias = 0.01*fsensor; %%rad/s
AngFieldBias = AngscaleBias*(2*rand()-1); %%% -1 to 1 (0 to 1) (0 to 2

%%
%Углы Эйлера (φ, θ, ψ) (roll, pitch, yaw) обычно не измеряются напрямую. Вместо этого:
%Акселерометр даёт вектор гравитации → позволяет оценить углы наклона (φ, θ).
%Магнитометр даёт вектор магнитного поля → позволяет оценить азимут (ψ).
%На основе этих векторов и текущей матрицы ориентации или фильтра (например, комплементарного) вычисляются углы Эйлера.
EulerBias = 2*pi/180*fsensor; %%rad
EulerBias = EulerBias*(2*rand()-1);