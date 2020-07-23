#Test all gene-SNP pairs and plot a histogram of all p-values
# Matrix eQTL by Andrey A. Shabalin
# http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
#
# Be sure to use an up to date version of R and Matrix eQTL.

options(error = quote({
      dump.frames(to.file=T, dumpto='last.dump')
        load('last.dump.rda')
        print(last.dump)
          q()
}))

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


## test if there is at least one argument: if not, return an error
#if (length(args)<3) {
#      stop("At least one argument must be supplied (input file).n", call.=FALSE)
#} else if (length(args)==3) {
#      # default output file
#      args[2] = "out.txt"
#}

# source("Matrix_eQTL_R/Matrix_eQTL_engine.r");
library(MatrixEQTL)

## Location of the package with the data files.
base.dir = find.package('MatrixEQTL');

## Settings

# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS

# Only associations significant at this level will be saved
#pvOutputThreshold = 1e-2;
pvOutputThreshold = 1;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");


## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "None"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);
print("Genotype file read")

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);
print("Expression file read")

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
#cvrt$fileSliceSize = 2000;
if(!is.null(covariates_file_name)) {
  cvrt$LoadFile(covariates_file_name);
  print("Covariates file read")
} else {
  print("No covariates file provided - skip reading")
}

## Run the analysis

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
  noFDRsaveMemory = FALSE);

#unlink(output_file_name);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
#show(me$all$eqtls)

## Plot the histogram of all p-values

plot(me)
