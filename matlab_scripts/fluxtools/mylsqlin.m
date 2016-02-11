function [X,LAMBDA,EXITFLAG] = MYLSQLIN(C,D,A,B,X0,varargin)
% MYLSQLIN Customized implementation of MATLAB's LSQLIN function.
%    X = MYLSQLIN(C,D,A,B,X0) solves the linear least-squares problem 
%   
%              min  0.5*(NORM(C*X-D)).^2       subject to    A*X <= B
%               X
%   
%    starting from the initial guess X0.
% 
%    X = MYLSQLIN(C,D,A,B,X0,AEQ,BEQ) applies the equality constraints 
%    AEQ*X = BEQ.
% 
%    [X,LAMBDA,EXITFLAG] = ... also returns the Lagrange multipliers and
%    exit flag. 
% 
%    See also LSQLIN, MYLINPROG.

EXITFLAG = -2;
X = [];
LAMBDA = [];

AEQ = [];
BEQ = [];
if nargin >5
    AEQ = varargin{1};
    BEQ = varargin{2};
end

nOpt = 10;
% nOpt = 100;

[beq_rows,beq_columns] = size(BEQ);
[b_rows,b_columns] = size(B);
[nMets_eq,nRxns_eq] = size(AEQ);
[nMets,nRxns] = size(A);
A = [AEQ;A];
b = [BEQ;B];
lb_eq(1:nRxns_eq,1) = -1000;
lb_in(1:nRxns,1) = -1000;
lb = [lb_eq;lb_in];
ub_eq(1:nRxns_eq,1) = 1000;
ub_in(1:nRxns,1) = 1000;
ub = [ub_eq;ub_in];
csense_eq(1:nMets_eq,1) = 'E';
csense_in(1:nMets,1) = 'L';
csense = [csense_eq; csense_in];
osense = 1;

NLPproblem.A = A;
NLPproblem.b = b;
NLPproblem.lb =lb;
NLPproblem.ub = ub;
NLPproblem.objFunction = 'mylsqlin_objFunction';
NLPproblem.gradFunction = 'mylsqlin_gradFunction';
NLPproblem.csense = csense;

% Current best solution (set for minimization)
% if (osense < 0)
%     currentSol.f = osense*inf;
% else
%     currentSol.f = osense*-inf;
% end
currentSol.f = osense*inf;
allObjValues = zeros(nOpt,1);
allSolutions = zeros(nRxns,nOpt);

%Define additional options for solveCobraNLP
majorIterationLimit = 1000;
printLevel = 0; %3 prints every iteration.  1 does a summary.  0 = silent.
NLPproblem.userParams.C = C; %pass the model into the problem for access by the nonlinear objective function
NLPproblem.userParams.D = D; %pass the model into the problem for access by the nonlinear objective function

NLPproblem.x0 = X0; %x0 now a cell within the NLP problem structure
for i = 1:nOpt
    solNLP = solveCobraNLP(NLPproblem, 'printLevel', printLevel, 'iterationLimit', majorIterationLimit);
    %solNLP = solveCobraNLP(NLPproblem, 'printLevel', printLevel, 'intTol', 1e-7, 'iterationLimit', majorIterationLimit); %New function call
    %solNLP = solveCobraNLP(NLPproblem,[],objArgs); Old Code
    fprintf('%d\t%f\n',i,osense*solNLP.obj);
    allObjValues(i) = osense*solNLP.obj;
    allSolutions(:,i) = solNLP.full;
    %if osense*solNLP.obj > currentSol.f
    if osense*solNLP.obj < currentSol.f
       currentSol.f = osense*solNLP.obj;
       currentSol.x = solNLP.full;
       currentSol.stat = solNLP.stat;
       X = solNLP.full;
       EXITFLAG = solNLP.stat;
       NLPproblem.x0 = solNLP.full;
    end
end