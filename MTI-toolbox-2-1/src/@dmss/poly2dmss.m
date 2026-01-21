function sys = poly2dmss(poly_str_cell, Ts, symbols)
% Torben Warnecke - 11/06/2024
arguments
    poly_str_cell string
    Ts (1,1) double
    symbols (1,:) string = ["xp", "x", "u", "y", "z"]
    % Can also also a sym-array e.g. [sym("dx"), sym("x"), sym("in"), sym("out"), sym("binary")] 
    % and will be converted into a string array by the input parser.
    % Also just parsing 4 symbols is possible, if dmss-model is not hybrid.
end
% Attention: Use multiplication sign "*"!

%% convert to char arrays
for k = 1:length(poly_str_cell)
    poly_str_cell{k} = char(poly_str_cell{k});
end

%% convert fraction
inequalityIndex = zeros([length(poly_str_cell),1]);
for k = 1:length(poly_str_cell)
    poly_sym = str2sym(poly_str_cell{k});    % convert string to symbolic equation
    poly_rhs = rhs(poly_sym);   % left hand side of the equation
    poly_sym = poly_sym - poly_rhs; % transfer left hand side to the ride side --> implicit equation
    if isequal(poly_sym, lhs(poly_sym) >= 0)
        poly_sym = lhs(poly_sym) * (-1);    % only using rhs, because lhs is zero and should stay so (no transfer from variables to the lhs is wished)
        inequalityIndex(k) = 1;
    elseif isequal(poly_sym, lhs(poly_sym) <= 0)
        poly_sym = lhs(poly_sym);   % only using rhs, because lhs is zero and should stay so (no transfer from variables to the lhs is wished)
        inequalityIndex(k) = 1;
    else
        poly_sym = lhs(poly_sym);
    end
    poly_sym = simplifyFraction(poly_sym); % simplifying fraction leads to the smallest amount of fraction left in representation
    poly_sym = expand(poly_sym);

    % Search for denominators
    [~, denominator] = numden(poly_sym);
    while denominator ~= 1
        poly_sym = poly_sym * denominator;
        poly_sym = simplify(poly_sym);
        poly_sym = expand(poly_sym);
        [~, denominator] = numden(poly_sym);
    end
    
    if inequalityIndex(k) == 1
        %poly_sym = 0 >= poly_sym;
        poly_str_cell{k} = append('0 >= ', char(poly_sym));
    else
        %poly_sym = 0 == poly_sym;
        poly_str_cell{k} = append('0 == ', char(poly_sym));
    end
    %poly_str_cell{i} = char(poly_sym);
end

%% legal variable signs, everything else will be handled as constants
variable_signs = symbols;

%% convert powers to auxiliary variables and equations (only integer powers are allowed, thus parameter-identification of powers is not possible)

% find maximum index of existing auxiliary (algebraic) variables
max_alg = regexp(poly_str_cell, {char('['+symbols(4)+']\d+')},'match');
for k = 1:length(max_alg)
    if isempty(max_alg{k}) == 0
        max_alg{k} = char(join(max_alg{k},''));%cell2mat(max_alg{k});
        index = regexp(max_alg{k},{'\d+'},'match');
        max_alg{k} = index{1};
    else
        max_alg{k} = '0';
    end
    max_alg{k} = max(str2double(max_alg{k}));
end
max_alg = max(cell2mat(max_alg));

% find all powers, check what are variables and what are not. Delete coefficient with powers from storage
powers = regexp(poly_str_cell, {'\w+\^\d+'},'match');
power_storage = [];
for k = 1:length(poly_str_cell)
    n = 1;
    while n <= length(powers{k})
        character = regexp(powers{k}(n),{'[a-zA-Z]'},'match');
        if ismember(character{1},variable_signs) == 0
            powers{k}(n) = [];
        else
            n = n+1;
        end
    end
    power_storage = [power_storage; powers{k}'];
end

% only safe unique strings
power_storage = unique(power_storage);

% sort storage descending by power
%power_storage = flip(sort(power_storage));

% sort storage ascending by power
power_storage = sort(power_storage);

% find variables in power_storage
variables_with_powers = string(regexp(power_storage, '[a-zA-Z]\d+','match'));
variables_with_powers = unique(variables_with_powers);

% begin replacing with highest powers first (e.g. "y1^3" gets replaced by
% "z1*y1^2")
auxiliary_eqn = {};
l = 1;
for k = 1:length(variables_with_powers)
    exponents = regexp(power_storage, variables_with_powers(k)+'\^\d+','match');
    exponents(cellfun(@isempty, exponents)) = {char()};
    exponents = string(exponents);
    exponents = regexp(exponents, '\^\d+','match');
    exponents(cellfun(@isempty, exponents)) = {char()};
    exponents = string(exponents);
    exponents = str2double(erase(exponents, "^"));

    while ~isempty(exponents) && sum(~isnan(exponents))>0
        max_alg = max_alg+1;
        aux = symbols(4) + string(max_alg);
        auxiliary_eqn{l,1} = {'0 == ' + aux + ' - ' + variables_with_powers(k)};
        l = l+1;

        for n = find(~isnan(exponents))'
            if exponents(n) > 2
                new_exponent = "^" + string(exponents(n)-1);
                rep_arg = aux + '*' + variables_with_powers(k) + new_exponent;
            elseif exponents(n) == 2
                rep_arg = aux + '*' + variables_with_powers(k);
            end
            poly_str_cell = strrep(poly_str_cell, variables_with_powers(k)+"^"+string(exponents(n)), rep_arg);
        end
        exponents = exponents -1;
        exponents(exponents==1) = [];
    end
end
% for k = 1:length(power_storage)
%     max_alg = max_alg+1;
%     aux = symbols(4) + string(max_alg);
%     variable = string(regexp(power_storage{k}, '[a-zA-Z]\d+','match'));
%     auxiliary_eqn{k,1} = {'0 == ' + aux + ' - ' + variable};
%     exponent = regexp(power_storage{k}, '\^\d+','match');
%     exponent = str2double(exponent{1}(2:end));
%     if exponent > 2
%         new_exponent = "^" + string(exponent-1);
%         rep_arg = aux + '*' + variable + new_exponent;
%     elseif exponent == 2
%         rep_arg = aux + '*' + variable;
%     end
%     poly_str_cell = strrep(poly_str_cell, power_storage{k}, rep_arg);
%     exponent = new_exponent;
% end
poly_str_cell = [poly_str_cell; auxiliary_eqn];

%% convert to string-cell-array
for k = 1:length(poly_str_cell)
    str_cell{k,1} = poly_str_cell(k);
end

% convert string-cell-array into iMTI object with str2iMTI
iMTI_str = str_cell;

sys = dmss.str2dmss(iMTI_str, Ts, symbols);%, symbols_str);
end