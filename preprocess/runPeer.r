library(peer)

args = commandArgs(trailingOnly = TRUE)
# Load expression
foldernm = args[1]
tissue = args[2]
# num factors
k = args[3]

exprnm = paste(foldernm, "/", tissue, "-normfilt.csv", sep = "")
factorsnm = paste(foldernm, "/", tissue, "_", k, "_peer_factors.tsv", sep = "")
weightsnm = paste(foldernm, "/", tissue, "_", k, "_peer_weights.tsv", sep = "")
precisionnm = paste(foldernm, "/", tissue, "_", k, "_peer_precision.tsv", sep = "")
residualsnm = paste(foldernm, "/", tissue, "_", k, "_peer_residuals.tsv", sep = "")
expr = read.csv(exprnm, header=TRUE, row.names=1)
expr = t(expr)


#k = 5

alphaprior_a=0.001
alphaprior_b=0.01
epsprior_a=0.1
epsprior_b=10
max_iter=1000

model = PEER()
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setAdd_mean(model, TRUE)
PEER_setPriorAlpha(model, alphaprior_a, alphaprior_b)
PEER_setPriorEps(model,epsprior_a, epsprior_b)
PEER_setNmax_iterations(model, max_iter)
PEER_setNk(model, k)
PEER_update(model)
factors = PEER_getX(model)
weights = PEER_getW(model)
precision = PEER_getAlpha(model)
residuals = PEER_getResiduals(model)

rownames(residuals) = rownames(expr)
colnames(residuals) = colnames(expr)
residuals = t(residuals)

write.table(factors, factorsnm, sep="\t")
write.table(weights, weightsnm, sep="\t")
write.table(precision, precisionnm, sep="\t")
write.table(residuals, residualsnm, sep="\t")
