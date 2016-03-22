
if ~exist('ign','var'); load ign; end;
in=ign; 

elems={'Kv';'SiO2';'TiO2';'Al2O3';'FeO';'MgO';'CaO';'Na2O';'K2O';'H2O_Total';'CO2';'tc1Crust'};


% Start with all Fe as FeO
in=feconversion(in);
in.FeO=in.FeOT;
in.Fe2O3=zeros(size(in.Fe2O3));

% % Set all H2O 4%
% in.H2O_Total=ones(size(in.H2O_Total))*4;

% Set undefined H2O to the average
in.H2O_Total(isnan(in.H2O_Total))=nanmean(in.H2O_Total);

% Set undefined CO2 to the average
in.CO2(isnan(in.CO2))=nanmean(in.CO2);

data=zeros(length(in.SiO2),length(elems));
for i=1:length(elems)
    data(:,i)=in.(elems{i});
end

% Reject samples with missing data
test=~any(isnan(data),2);

% Reject samples with suspicious anhydrous normalizations
test=test &~ (sum(data(:,2:9),2)>100) &~ (sum(data(:,2:9),2)<90);

data=data(test,:);

size(data)

dlmwrite('ignmajors.csv',data,',')
