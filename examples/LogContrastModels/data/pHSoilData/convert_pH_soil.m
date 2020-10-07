% Script that converts the csv files from R into MATLAB .mat files

%% Read the initial pH data
Y_csv = readtable('pHData.csv')

% Convert the second column to numerical data
Y = table2array(Y_csv(:,2));

%% Read and convert the initial OTU count data
OTU_csv = readtable('soilOTUData.csv')

% OTU identifiers
OTUlabels = table2cell(OTU_csv(1,2:end));

% Sample identifiers
samplelabels = table2cell(OTU_csv(2:end,1));

% Convert table entries OTU count data 

% Convert it to cell first
X_cell = table2array(OTU_csv(2:end,2:end));

% Convert strings to double
X_count = str2double(X_cell);

%% Read the taxonomic table
tax_csv = readtable('pHTaxTable.csv','Format','%s%s%s%s%s%s%s%s')

% Transform chars back to strings
[p,r] = size(tax_csv);

tax_cell = cell(p,r);

for i=1:p
    i
    for j=1:r
        j
        currEntry = tax_csv{i,j};
        bla = currEntry{1};
        newEntry = bla(2:end-1);
        if j==1
            tax_cell{i,j} = str2num(newEntry);
        else
            tax_cell{i,j} = newEntry;
        end
    end
end

tax_table = cell2table(tax_cell);

tax_table.Properties.VariableNames{'tax_cell1'} = 'OTU_ID';
tax_table.Properties.VariableNames{'tax_cell2'} = 'kingdom';
tax_table.Properties.VariableNames{'tax_cell3'} = 'phylum';
tax_table.Properties.VariableNames{'tax_cell4'} = 'class';
tax_table.Properties.VariableNames{'tax_cell5'} = 'order';
tax_table.Properties.VariableNames{'tax_cell6'} = 'family';
tax_table.Properties.VariableNames{'tax_cell7'} = 'genus';
tax_table.Properties.VariableNames{'tax_cell8'} = 'species';


save('pHSoilData.mat','tax_table','X_count','Y','OTUlabels','samplelabels')





