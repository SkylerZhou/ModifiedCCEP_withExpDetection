function T = csv2tbl(csv_dir)


% Read the CSV into a table. 'TextType' set to 'string' ensures text is loaded as string arrays.
T = readtable(csv_dir, 'TextType', 'string');
    
% Get the list of variable names (columns)
varNames = T.Properties.VariableNames;
    
% Loop over each column
for i = 1:numel(varNames)
    colData = T.(varNames{i});
        
    % If the column is a string array or a char array, convert it to a cell array.
    if isstring(colData)
        colData = cellstr(colData);
    elseif ischar(colData)
        colData = cellstr(colData);
    elseif ~iscell(colData)
        continue;
    end
        
    % Update the table column to be a cell array.
    T.(varNames{i}) = colData;
        
    % Process each cell in the column
    for j = 1:height(T)
        value = T.(varNames{i}){j};
            
        % Check if the cell is non-empty and starts with '[' and ends with ']'
        if ~isempty(value) && value(1)=='[' && value(end)==']'
                
            value_inner = value(2:end-1); % remove []
            value_inner = strrep(value_inner, '''', ''); % remove ''
                
            % Split the remaining string by commas (and trim extra whitespace)
            parts = strsplit(value_inner, ',');
            parts = strtrim(parts);

            % If there is only one part that is empty, assign an empty cell array
            if numel(parts)==1 && isempty(parts{1}) 
                parts = {};
            end
            
            % Replace the cell content with the cell array of strings 
            T.(varNames{i}){j} = parts;
        end
    end
end

end