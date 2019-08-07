function beq = projCeq(b,PCeq,Ceq,rhsvec)
% Projection operator for equality constraints of the 
% type Ceq b = rhsvec 

% Standard formula for a given Ceq in R^kxp with k<p
%beq = b-Ceq'*inv(Ceq*Ceq')*(Ceq*b-rhsvec);

% Formula for a given Ceq in R^kxp with k>p
%beq = b-Ceq'*pinv(Ceq*Ceq')*(Ceq*b-rhsvec);

beq = b-PCeq*(Ceq*b-rhsvec);


