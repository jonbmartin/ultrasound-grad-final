sig = 1;
N=20;
M=20;
X = ones([M N]);

u = optimvar('u');
w = optimvar('w',M);
r = optimvar('r',M);
prob = optimproblem('Objective',u,'ObjectiveSense','min');
prob.Constraints.c1 = -u*ones([N 1])<= X'*w;
prob.Constraints.c2 = X'*w <= u*ones([N 1]);
prob.Constraints.c3 = w >= 1 + sig*r;
prob.Constraints.c4 = -r <= w;
prob.Constraints.c5 = w <= r;

problem = prob2struct(prob);

[sol,fval,exitflag,output] = linprog(problem)