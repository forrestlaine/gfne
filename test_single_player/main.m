% Setup and Solve Game

% Define Primal Variables
varList = {};
for t = 0:9
    name = join(['x',string(t)],'');
    isPrimal = true;
    playerID = 0;
    dim = 4;
    varList = [varList; {name,isPrimal, playerID, dim}];
    
    name = join(['u',string(t)],'');
    isPrimal = true;
    playerID = 0;
    dim = 2;
    varList = [varList; {name,isPrimal, playerID, dim}];
end
name = join(['x',string(10)],'');
isPrimal = true;
playerID = 0;
dim = 4;
varList = [varList; {name,isPrimal, playerID, dim}];

playerID = 0;
cost = gameCostFunction(varList,varList,playerID,@genCost);
constraintDim = 44;
stageID = 0;
eqConstraints = gameConstraintFunction(varList, ...
                                       stageID, ...
                                       playerID, ...
                                       constraintDim, ...
                                       @genEqualityConstraints);
                                  
constraintDim = 0;
ineqConstraints = gameConstraintFunction(varList, ...
                                         stageID, ...
                                         playerID, ...
                                         constraintDim, ...
                                         @genInequalityConstraints);
                                    
playerFunctions = {cost, eqConstraints, ineqConstraints};                                     
varTable = variableTable(playerFunctions);

lag = lagrangian(0,0,cost,eqConstraints,ineqConstraints);
% global evalLag;
evalLag = lag.genEval(varTable);


% func = @cpfunjac;
% gradEvalLag = lag.genGradEval(varTable);
z = rand(108,1);

lagValue = evalLag(z)
% za = amatinit(z,2);
% grad = ajac(evalLag(za),0);
% hess = ahess(evalLag(za),0);

% lagGrad = gradEvalLag(z)';

% eval = eqConstraints.genEval(varTable);
% grad2 = eval(z);

% tic
% for t = 1:100
% evalLag(z);
% end
% toc / 100

% function [f,J,domerr] = cpfunjac(z,jacflag)
%     global evalLag
%     za = amatinit(z,2);
%     f = ajac(evalLag(za),0);
%     if jacflag
%         J = sparse(ahess(evalLag(za),0));
%     end
%     domerr = 0;
% end
                                    


