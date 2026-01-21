function x = getNextState(sys,xu)
    x = sys.data(1).processXU(xu);
end