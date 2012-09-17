function g = rho_dct(I)

temp1=dct2(I);
temp2=temp1(:);
temp3=temp2(2:end);

%g=kurtosis(temp3);
g=bliinds2.rho_gen_gauss(temp3);