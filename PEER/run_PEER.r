library(peer)

args = commandArgs(trailingOnly = TRUE)
foldernm = args[1]

filenm = paste(foldernm, "/expression.csv", sep = "")
factorsnm = paste(foldernm, "/peer_factors.tsv", sep = "")
weightsnm = paste(foldernm, "/peer_weights.tsv", sep = "")
precisionnm = paste(foldernm, "/peer_precision.tsv", sep = "")
residualsnm = paste(foldernm, "/peer_residuals.tsv", sep = "")

expr = read.csv(filenm, header = TRUE, row.names = 1, sep='\t')
expr = t(expr)

k <- 10

model = PEER()
PEER_setPhenoMean(model, as.matrix(expr))
PEER_setAdd_mean(model, TRUE)
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


