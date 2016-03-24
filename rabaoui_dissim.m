% implements dissimilarity metric from Rabaoui et al 2008 (DOI 10.1109/TIFS.2008.2008216)
% model: from libsvm svmtrain()
% vec: MxN normalized feature vectors (M vectors, N features)
function prob = rabaoui_dissim(model, vecs)
    total = zeros(size(vecs,1),1);
    for i=1:length(model.sv_coef)
        total = total + model.sv_coef(i) .* rbf(vecs, model.SVs(i,:), model.Parameters(4));
    end

    prob = -log(total) + log(model.rho);
end

function val = rbf(xs, y, gamma)
    val = exp(-gamma*sum(bsxfun(@minus, xs, y).^2, 2));
end
