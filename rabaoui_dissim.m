% implements dissimilarity metric from Rabaoui et al 2008 (DOI 10.1109/TIFS.2008.2008216)
% model: from libsvm svmtrain()
% vec: 1xN normalized feature vector
function prob = rabaoui_dissim(model, vec)
    total = 0;
    for i=1:length(model.sv_coef)
        total = total + model.sv_coef(i) * rbf(vec, model.SVs(i,:), model.Parameters(4));
    end

    prob = -log(total) + log(model.rho);
end

function val = rbf(x, y, gamma)
    val = exp(-gamma*sum((x - y).^2));
end
