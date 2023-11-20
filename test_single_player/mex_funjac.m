function [F, J, domerr] = mex_funjac(z, jacflag)
    mex gen.c -largeArrayDims
    F = gen('gamegrad',z);
    J = gen('gamehess',z);
    domerr = 0;
end