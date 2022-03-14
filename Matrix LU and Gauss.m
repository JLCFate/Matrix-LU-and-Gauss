%LU Grzegorz Gra¿ewicz s183208
A=[7.290, -1.600, 0.730, -1.180;
   -0.220, 4.000, -9.070, -5.870;
   4.480, 0.640, 2.480, 6.730;
   6.490, -9.390, -4.400, 6.180];
   
B=[8.840; -75.160; 55.9600; -20.930];

%Metoda - Eliminacja Gaussa

N = rank(A);
% Rozwi¹zanie za pomoc¹ metody Gaussa
MA=A;
MB=B;
for ur=1:N-1
for r=1+ur:N
k=GA(r,ur)/MA(ur,ur);
for c=ur:N
MA(r,c)=GA(r,c)-k*MA(ur,c);
end
MB(r)=MB(r)-k*MB(ur);
end
end
for r=N:-1:1
p=0;
for c=r+1:N
q=MA(r,c)*MX(c);
p=p+q;
end
MX(r)=((MB(r)-p)/MA(r,r));
end
% Koniec Metody Gaussa
% Rozwiazanie metoda LU poczatek
% Tworzenie macierzy L i macierzy U
L=eye(N);
U=zeros(N);
U=A;
for s=1:(N-1)
n=s+1;
for i=n:N
L(i,s)=U(i,s)/U(s,s);
for k=s:N
U(i,k)=U(i,k)-L(i,s)*U(s,k);
end
end
end
y=L\B;
LUX=U\y;
% LU KONIEC
disp('Gauss:');
MX
disp('Sprawdzenie (A*X-B):');
(A*MX'-B)'
disp('LU:');
LUX'
disp('Sprawdzenie (A*X-B):');
(A*LUX-B)'
disp('LU:');
(A\B)'
disp('Sprawdzenie (A*X-B):');
(A*(A\B)-B)'