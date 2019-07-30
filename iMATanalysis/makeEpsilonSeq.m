function [epsilon_f, epsilon_r] = makeEpsilonSeq(model, reactions, epsilon0, k)
%default is take epsilon_f || epsilob_r = min(epsilon, Vmax_f || -Vmax_r)
epsilon_f = epsilon0 * ones(length(reactions),1);
epsilon_r = epsilon0 * ones(length(reactions),1);
for i = 1: length(reactions)
    testModel = model;
    testModel = changeObjective(testModel,reactions(i));
    Vm_f = optimizeCbModel(testModel, 'max');
    if Vm_f.obj < epsilon0 && Vm_f.obj > 0
        epsilon_f(i) = Vm_f.obj * k;
    end
    Vm_r = optimizeCbModel(testModel,'min');
    if -Vm_r.obj < epsilon0 && Vm_r.obj < 0
        epsilon_r(i) = -Vm_r.obj * k;
    end
end
