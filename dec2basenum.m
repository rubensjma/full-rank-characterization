function repBase = dec2basenum(D,B,N)
% D = Number to be expressed in base B
% N = Number of places in the representation

Daux = D;
for i = 1:N
    a(i) = rem(Daux,B);
    Daux = (Daux - a(i))/B;
end

repBase = a(end:-1:1);