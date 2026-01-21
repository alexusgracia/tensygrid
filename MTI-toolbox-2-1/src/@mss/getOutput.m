function y = getOutput(sys,xu)
    if(isempty(sys.G))
        y = xu(1:sys.n); % This behaves like a "unity"tensor. Is not currently in use, since it moves data around unneccessarily
    else
        y = sys.data(2).processXU(xu);
    end
end