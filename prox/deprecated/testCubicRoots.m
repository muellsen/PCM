% Test cubic roots

n=1e4
const = 100

pPoints = -const+2*const*rand(n,1);
qPoints = const+rand(n,1);

rootsMat = zeros(n,3);
posRealMat = zeros(n,3);
negRealMat = zeros(n,3);

for i=1:n
    p = pPoints(i);
    q = qPoints(i);
    rootfacs = [1 0 p q];
    cubroots = roots(rootfacs)
    rootsMat(i,:) = cubroots;
    for j=1:3
        posRealMat(i,j) = (real(cubroots(j))>0 && imag(cubroots(j))==0);
        negRealMat(i,j) = (real(cubroots(j))<0 && imag(cubroots(j))==0);
    end
end

figure;
plot(1:n,sum(posRealMat'),'.','MarkerSize',20)
xlabel('index i')
ylabel('# root types')
grid on
hold on
plot(1:n,sum(negRealMat'),'.','MarkerSize',20)
legend({'#positive real roots','#negative real roots'})


nPosRealRoots = sum(sum(posRealMat'))



% French Wikipedia example
p = -3;
q = 1;
rootfacs = [1 0 p q]
cubroots = roots(rootfacs)






