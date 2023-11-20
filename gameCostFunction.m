classdef gameCostFunction
    %GAMEFUNCTION represents a function constituting some part 
    % of a differential game.
    
    properties 
        varList
        decisionVarList
        playerIndex
        genEvalFunc
    end
    
    methods
        function obj = gameCostFunction(varList, ...
                                    decisionVarList, ...
                                    playerIndex, ...
                                    genEvalFunc)
            obj.varList = varList;
            obj.decisionVarList = decisionVarList;
            obj.playerIndex = playerIndex;
            obj.genEvalFunc = genEvalFunc;
        end
        
        function eval = genEval(obj, Index)
            eval = obj.genEvalFunc(Index);
        end
        
    end
    
end

