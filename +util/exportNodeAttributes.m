function exportNodeAttributes(filename, attributeName, specID, attValues)
%{
    Writes a cytoscape node attributes (species) file. Only takes
    attributes that are numerical values (no strings, class=Double)
%}
    attributeName2 = attributeName;
    attributeName = {attributeName, '(class=Double)', ''};
    attribute = attValues;

    equals = {};
    for i = 1:length(attribute)
        equals{end+1} = '=';
        attribute2{i} = num2str(attribute(i));
    end
    attribute = attribute2';

    temp = [specID, equals', attribute];
    towrite = vertcat(attributeName, temp);
    fname = sprintf(filename, attributeName2);
    textwrite(fname,towrite); 
    disp(['Wrote ',cd,'\',fname]);

end

    