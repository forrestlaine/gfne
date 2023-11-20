classdef variableTable
    
    properties
        varTable
        indexedTable
    end
    
    methods
        function obj = variableTable(playerFunctions)
            obj.varTable = table('Size',[0,4], ...
                'VariableTypes', {'string','logical','uint8','uint8'},...
                'VariableNames',{'name','isPrimal','player','dimension'});
            
            for i = 1:size(playerFunctions,1)
                cost = playerFunctions{i,1};
                eqConstraints = playerFunctions{i,2};
                ineqConstraints = playerFunctions{i,3};
                obj = obj.updateTable(cost.varList);
                obj = obj.updateTable(eqConstraints.varList);
                obj = obj.updateTable(ineqConstraints.varList);
                eqDuals = {join(['lambda',string(eqConstraints.stageIndex)],''), ...
                           false, ...
                           eqConstraints.playerIndex, ...
                           eqConstraints.constraintDim};
                ineqDuals = {join(['mu',string(ineqConstraints.stageIndex)],''), ...
                             false, ...
                             ineqConstraints.playerIndex, ...
                             ineqConstraints.constraintDim};
                obj = obj.updateTable([eqDuals; ineqDuals]);
            end
            obj = obj.generateIndexedTable();
        end
        
        function obj = updateTable(obj, gameVariableList)
            % Assumes gameVariableList is of the form: 
            % {name0, bool0, id0, dim0; 
            %  name1, bool1, id1, dim1;
            %  ... };
            obj.varTable = unique([obj.varTable; gameVariableList]);
        end
        
        function obj = generateIndexedTable(obj)
            dims = obj.varTable.dimension;
            i2 = cumsum(dims);
            i1 = [0;i2(1:end-1)]+1;
            
            idTable = table(i1,i2,'VariableNames',{'startIndex','endIndex'});
            obj.indexedTable = [idTable, obj.varTable];
        end
        
        function var = getVar(obj, name, player, z)
            row = obj.indexedTable(and(obj.indexedTable.name == name, ...
                                       obj.indexedTable.player == player),:);
            if size(row,1) > 1
                error('Duplicate variables exist');
            elseif size(row,1) < 1
                error('No variable matches name and player id');
            else
                var = z(row.startIndex:row.endIndex);
            end
        end
        
        function zc = combine(obj, primals, z, player)
            zc = z;
            rows = and(obj.indexedTable.player == player, ...
                       obj.indexedTable.isPrimal == true);
            primalTable = obj.indexedTable(rows,:);
            starts = primalTable.startIndex;
            ends = primalTable.endIndex;
            dims = primalTable.dimension;
            total = 1;
            for i = 1:size(primalTable,1)
                dim = dims(i);
                zc(starts(i):ends(i)) = primals(total:total+dim-1);
                total = total+dim;
            end
        end
        
        function primals = extract(obj, z, player)
            rows = and(obj.indexedTable.player == player, ...
                       obj.indexedTable.isPrimal == true);
            primalTable = obj.indexedTable(rows,:);
            starts = primalTable.startIndex;
            ends = primalTable.endIndex;
            dims = primalTable.dimension;
            primals = zeros(sum(dims),1);
            total = 1;
            for i = 1:size(primalTable,1)
                dim = dims(i);
                primals(total:total+dim-1) = z(starts(i):ends(i));
                total = total+dim;
            end
        end
    end
end

