function sys = str2dmss(str_cell, Ts, symbols)
% Torben Warnecke - 11/06/2024
arguments
    str_cell cell
    Ts (1,1) double
    symbols (1,:) string = ["xp", "x", "u", "y", "z"]
    % Can also also a sym-array e.g. [sym("dx"), sym("x"), sym("in"), sym("out"), sym("binary")] 
    % and will be converted into a string array by the input parser.
    % Also just parsing 4 symbols is possible, if dmss-model is not hybrid.
end


% Create a iMTI model out of a string-cell-array and timestep size
% (0 for continous).
% example: {"0 = f1xp1 + f2x1 + f3u1 + f4x1u1"; "0 = f5y1 + f6x1"}
number_eq = length(str_cell);

% delete " ", "0=" and "*" and split at "+" and "-":
inequalityIndex = zeros([number_eq,1]);
for k = 1:number_eq
    str_cell{k} = erase(str_cell{k}, " ");

    str_cell{k} = erase(str_cell{k}, "0=");
    str_cell{k} = erase(str_cell{k}, "="); % erase "=", if from logical statement ("==") an additional equal sign is left!

    char_array = convertStringsToChars(str_cell{k}{1});
    if matches(char_array(1:2), '0>')
        char_array = erase(char_array, '0>');
        inequalityIndex(k) = 1;
    end

    lb = 0; % left bracket counter
    rb = 0; % right bracket counter
    ind = []; % indices of plus and minus signs outside of brackets
    l = 1;
    n = length(char_array);
    while l<=n
        if char_array(l) =='('
            lb = lb + 1;
            l = l+1;
        elseif char_array(l) == ')'
            rb = rb + 1;
            l = l+1;
        elseif (char_array(l)=='+')&&(lb == rb)%&&char_array(l-1)~='e'
            char_array(l) = '';
            ind = [ind ,l-1];
            n = n-1;
        elseif (char_array(l)=='-')&&(lb == rb)&&l>1%&&char_array(l-1)~='e'
            ind = [ind ,l-1];
            l = l+1;
        else
            l = l+1;
        end
    end
    ind_distance = [ind(1:end), length(char_array)] - [0, ind(1:end)];
    str_cell{k} = convertCharsToStrings(mat2cell(char_array(1:sum(ind_distance)),1,ind_distance));
end

% search and delete variables and create H.factor
Hfactor = [];
Hfactorstates = [];

variable_signs = symbols;

numbers = zeros(1, length(variable_signs));
for n = 1:length(variable_signs)
    sumc = 1;
    l = 1;
    container = [];
    statecontainer = [];
    if n == 2
        Hfactor = [Hfactor; Hfactorstates];
    end
    while sumc > 0
        Hfactor = [Hfactor; container];
        Hfactorstates = [Hfactorstates; statecontainer];

        container = [];
        statecontainer = [];
        for k = 1:number_eq
            if n == 1
                % For states, searching for xp and x must happen the same time,
                % e.g. if xp1 exist there needs to be a line for x1 in Hfactor, even if its all zeros!
                pattern = regexpPattern(variable_signs(n)+sprintf("%d",l)+"(?!\d)"+"|"+variable_signs(n)+sprintf("%d",l)+"$");
                %pattern = variable_signs(n)+sprintf("%d",j);
                container = [container, contains(str_cell{k}, pattern)];
                %str_cell{i} = replace(str_cell{i}, '*'+pattern+'*', '*');
                str_cell{k} = erase(str_cell{k}, '*'+pattern); % delete variables from cell-array
                str_cell{k} = erase(str_cell{k}, pattern+'*'); % delete variables and unnecessary multiplication sign from cell-array
                str_cell{k} = erase(str_cell{k}, pattern); % delete variables from cell-array
                %str_cell{i} = erase(str_cell{i}, pattern); % delete variables from cell-array
                str_cell{k}(str_cell{k}=="") = "1"; % where no variable and parameter is left, add 1 to the term
                str_cell{k}(str_cell{k}=="-") = "-1";

                statepattern = regexpPattern(variable_signs(n+1)+sprintf("%d",l)+"(?!\d)"+"|"+variable_signs(n+1)+sprintf("%d",l)+"$");
                %statepattern = variable_signs(n+1)+sprintf("%d",j);
                statecontainer = [statecontainer, contains(str_cell{k}, statepattern)];
                %str_cell{i} = replace(str_cell{i}, '*'+statepattern+'*', '*');

                str_cell{k} = erase(str_cell{k}, '*'+statepattern); % delete variables and unnecessary multiplication sign from cell-array
                str_cell{k} = erase(str_cell{k}, statepattern+'*'); % delete variables from cell-array
                str_cell{k} = erase(str_cell{k}, statepattern); % delete variables from cell-array

                str_cell{k}(str_cell{k}=="") = "1"; % where no variable and parameter is left, add 1 to the term
                str_cell{k}(str_cell{k}=="-") = "-1";
            elseif n >= 3
                pattern = regexpPattern(variable_signs(n)+sprintf("%d",l)+"(?!\d)"+"|"+variable_signs(n)+sprintf("%d",l)+"$");
                %pattern = variable_signs(n)+sprintf("%d",j);
                container = [container, contains(str_cell{k}, pattern)];
                str_cell{k} = replace(str_cell{k}, '*'+pattern+'*', '*');

                str_cell{k} = erase(str_cell{k}, '*'+pattern); % delete variables and unnecessary multiplication sign from cell-array
                str_cell{k} = erase(str_cell{k}, pattern+'*'); % delete variables from cell-array
                str_cell{k} = erase(str_cell{k}, pattern); % delete variables from cell-array

                str_cell{k}(str_cell{k}=="") = "1"; % where no variable and parameter is left, add 1 to the term
                str_cell{k}(str_cell{k}=="-") = "-1";
            end
        end
        l = l+1;
        sumc = sum(container);
    end
    numbers(n) = l-2;
end
%H.factor = sparse(Hfactor) > 0;

% number of variables
n = max(numbers(1:2));      %n
m = numbers(3);      %m
p = numbers(4);    %p
if length(variable_signs)>4
    q = numbers(5);    %q
else
    q = 0;
end

H = hyCPN1();
H.F.stateDerivative = Hfactor(1:n,:);
H.F.state = Hfactor(n+1:2*n,:);
H.F.input = Hfactor(2*n+1:2*n+m,:);
H.F.algebraic = Hfactor(2*n+m+1:2*n+m+p,:);
H.F.boolean = Hfactor(2*n+m+p+1:2*n+m+p+q,:);

% create symbloic H.phi from left over string cell
Hphi = [];
for k = 1:number_eq
    Hphi = [ Hphi ,[zeros( (k-1)*((k-1)>0) , length(str_cell{k}));...
        str2sym(str_cell{k});...
        zeros(number_eq-k, length(str_cell{k}))] ];
end

if isempty(symvar(Hphi)) == 1
    % check if H.phi has symbolic arguments or only consist of
    % numbers
    Hphi = double(Hphi);
end

H.phi.equality = Hphi(~inequalityIndex,:);
H.phi.inequality = Hphi(logical(inequalityIndex),:);

sys = dmss(H, Ts);
end