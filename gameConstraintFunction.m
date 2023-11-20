classdef gameConstraintFunction
    %GAMEFUNCTION represents a function constituting some part 
    % of a differential game.
    
    properties 
        varList
        stageIndex
        playerIndex
        constraintDim
        genEvalFunc
    end
    
    methods
        function obj = gameConstraintFunction(varList, ...
                                    stageIndex, ...
                                    playerIndex, ...
                                    constraintDim, ...
                                    genEvalFunc)
            obj.varList = varList;
            obj.stageIndex = stageIndex;
            obj.playerIndex = playerIndex;
            obj.constraintDim = constraintDim;
            obj.genEvalFunc = genEvalFunc;
        end
        
        function eval = genEval(obj, Index)
            eval = obj.genEvalFunc(Index);
        end
        
    end
    
end

