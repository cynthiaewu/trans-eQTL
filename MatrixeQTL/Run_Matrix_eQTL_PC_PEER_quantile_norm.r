library("MatrixEQTL");

suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser()
parser$add_argument("-g", "--genotype", help="Input genotype file", required=TRUE)
parser$add_argument("-e", "--expression", help="Input expression file", required=TRUE)
parser$add_argument("-o", "--output", help="Output file", required=TRUE)
parser$add_argument("-c", "--covariates", help="Input covariates file")
args <- parser$parse_args()
# args <- parser$parse_args(c('--genotype', '--expression', '--output'))
#accumulate_fn <- get(args$accumulate)
#print(accumulate_fn(args$integers))
SNP_file_name <- args$genotype
expression_file_name <- args$expression
covariates_file_name <- args$covariates
output_file_name <- args$output
#print(c(SNP_file_name, expression_file_name, output_file_name))

#args = commandArgs(trailingOnly = TRUE)
#SNP_file_name = args[1]
#expression_file_name = args[2]
#covariates_file_name = args[3]
#output_file_name = args[4]

useModel = modelLINEAR;
#SNP_file_name = "./GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv_chr1_vcftools_addedfilter.intersect.table";
#expression_file_name = "./Nerve-Tibial_Euro_PEER_residual.intersect.tsv";
#covariates_file_name = character();
#covariates_file_name = "./gtex_v8_covariates_dropTRIDTH.intersect.tsv";
#output_file_name = "./matrixeqtltest.out";
pvOutputThreshold = 1;
errorCovariance = numeric();

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
#snps$fileOmitCharacters = "-1"; # denote missing values;
snps$fileOmitCharacters = "."; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in pieces of 2,000 rows
snps$LoadFile( SNP_file_name );

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

# Quantile normalization
for( sl in 1:length(gene) ) {
    mat = gene[[sl]];
    mat = t(apply(mat, 1, rank, ties.method = "average"));
    mat = qnorm(mat / (ncol(gene) + 1));
    gene[[sl]] = mat;
}
rm(sl, mat);

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$fileSliceSize = 2000;
if(length(covariates_file_name) > 0) {
  cvrt$LoadFile(covariates_file_name);
}

me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = TRUE);


