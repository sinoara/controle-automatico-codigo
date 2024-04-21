%% Controle latero-direcional
cessna182

[num, den] = ss2tf(Ald,Bld,Cld,Dld, 1);

rlocus(tf(num(3, :), den))

Kp = 12;
Kr = 10;

print('-smodel_cessna182_ld', '-dpng', 'cessna_ld.png');
print('-smodel_cessna182_ld', '-dsvg', 'cessna_ldg.svg');