

A=dlmread('inputint.txt');
NumAng=A(1);




% finput=fopen('Sinograma.txt', 'r');
R=dlmread('Sinograma.txt');
deltatheta=360/NumAng;
B=transpose(R);

% subplot(1,2,1), imshow(R);
% subplot(1,2,2), imshow(B);
thetabun=180/pi*transpose(dlmread('AngulosBundle.txt'));
theta=0:deltatheta:(360-deltatheta);


I=iradon(B, theta);
str=sprintf('Reconstruccion con %d angulos y %d rayos cada ángulo',NumAng, A(2));
imshow(I),title(str);
savefig(str);
