classdef lagrangian
    %UNTITLED13 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        playerIndex
        cost
        equalityConstraints
        inequalityConstraints
        stage
    end
    
    methods
        function obj = lagrangian(stage, ...
                                  playerIndex, ...
                                  cost, ...
                                  equalityConstraints,...
                                  inequalityConstraints)
            obj.stage = stage;
            obj.playerIndex = playerIndex;
            obj.cost = cost;
            obj.equalityConstraints = equalityConstraints;
            obj.inequalityConstraints = inequalityConstraints;
            
        end
        
        function eval = genEval(obj,Index)
            function val = evaluate(z)
                costFn = obj.cost.genEval(Index);
                eqConst = obj.equalityConstraints.genEval(Index);
                ineqConst = obj.inequalityConstraints.genEval(Index);
                lam = Index.getVar(join(['lambda',string(obj.stage)],''),obj.playerIndex,z);
                mu = Index.getVar(join(['mu',string(obj.stage)],''),obj.playerIndex,z);
                val = costFn(z) - lam'*eqConst(z) - mu'*ineqConst(z);
            end
            eval = @evaluate;
        end
        
        
        function gradEval = genGradEval(obj,Index)
            costFn = obj.cost.genEval(Index);
            eqConst = obj.equalityConstraints.genEval(Index);
            ineqConst = obj.inequalityConstraints.genEval(Index);
            function grad = evaluateGrad(z)
                
                function val = evaluate(z)
                    lam = Index.getVar(join(['lambda',string(obj.stage)],''),obj.playerIndex,z);
                    mu = Index.getVar(join(['mu',string(obj.stage)],''),obj.playerIndex,z);
                    val = costFn(z) - lam'*eqConst(z) - mu'*ineqConst(z);
                end
                
                za = amatinit(z);
                grad = ajac(evaluate(za),0);
            end
            gradEval = @evaluateGrad;
        end
    end
end

