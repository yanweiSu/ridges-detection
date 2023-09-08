T=1000;
K=3;
C = ones(K,1);
C=C*C' + eye(K);
parameters = sqrt([.05 .93]);
p = 1;
q = 1;
type='Scalar';
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);

parameters = sqrt([.05 .05 .88]);
p=2;
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);
parameters = sqrt([.05 .05 .4 .4]);
q=2;
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);
q = 0;
parameters = sqrt([.2 .2]);
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);


p = 1;
q = 1;
type='Diagonal';
parameters = sqrt([.05 .93]);
parameters = repmat(parameters,K,1);
parameters = parameters(:);
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);

parameters = sqrt([.05 .05 .88]);
parameters = repmat(parameters,K,1);
parameters = parameters(:);
p=2;
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);
parameters = sqrt([.05 .05 .4 .4]);
parameters = repmat(parameters,K,1);
parameters = parameters(:);
q=2;
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);
q = 0;
parameters = sqrt([.2 .2]);
parameters = repmat(parameters,K,1);
parameters = parameters(:);
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);



p = 1;
q = 1;
type='CP';
parameters = sqrt([.05*ones(1,K) .98]);
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);
parameters = sqrt([.05*ones(1,K) .05*ones(1,K) .88]);
p=2;
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);
q = 0;
parameters = sqrt([.2*ones(1,K) .2*ones(1,K)]);
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);





K=3;
T=randn(1000,K);
C = ones(K,1);
C=C*C' + eye(K);
parameters = sqrt([.05 .93]);
p = 1;
q = 1;
type='Scalar';
[data, Ht] = rarch_simulate(T,C,parameters,p,q,type);
backCast = eye(K);
type = 1;
isJoint = false;
T = 1000;
data2 = zeros(K,K,T);
for i=1:T
    data2(:,:,i) = data(i,:)'*data(i,:);
end
isCChol = false;
[ll, lls, Ht2] = rarch_likelihood(parameters,data2,p,q,C,backCast,type,isJoint,isCChol);


[c,ceq] = rarch_constraint(parameters,data2,p,q,C,backCast,type,isJoint,isCChol)


sv = parameters;
LB = -ones(size(sv));
UB = ones(size(sv));
options = optimset('fmincon')
fmincon(@rarch_likelihood,sv,[],[],[],[],LB,UB,@rarch_constraint,options,data2,p,q,C,backCast,type,isJoint,isCChol);


parameters = sqrt([.05*ones(1,K) .93*ones(1,K)]);
sv = parameters;
UB = ones(size(sv)) * .99998;
LB = -BB;
options = optimset('fmincon');
type = 3;
options.Display = 'iter';
fmincon(@rarch_likelihood,sv,[],[],[],[],LB,UB,@rarch_constraint,options,data2,p,q,C,backCast,type,isJoint,isCChol);

profile on
rarch_likelihood(sv,data2,p,q,C,backCast,type,isJoint,isCChol);
profile report
profile off

method = '2-stage';
type = 'Scalar';
[parameters, ll, Ht, VCV, scores] = rarch(data,p,q,type,method);