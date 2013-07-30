function g = gama_dct(I)

temp1=dct2(I);
temp2=temp1(:);
temp3=temp2(2:end);

%g=kurtosis(temp3);
g=bliinds2.gama_gen_gauss(temp3);