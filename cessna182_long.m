%% Controle de arfagem
cessna182

[num, den] = ss2tf(Along,Blong,Clong,Dlong);

Kc = 0.0741;
%Kc = 4.22;
rlocus(tf(num(3, :), den))

print('-smodel_cessna182_long', '-dpng', 'cessna_long.png');
print('-smodel_cessna182_long', '-dsvg', 'cessna_long.svg');