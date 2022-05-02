#' Function for splitting jobs
#' @param x joblist
#' @param n number of total file
#' @keywords split job
#' @export
#' @examples
#' file_helper()
file_helper <- function(x,n){
  split(x, cut(seq_along(x), n, labels = FALSE))
}


#' Function for filling in missing genotype
#' @param geno genotype dataframe
#' @param type filling missing value
#' @keywords constant method 
#' @export
#' @examples
#' geno_fill()
geno_fill <- function(geno, type){

  ##Filling in missing data##
  geno = data.matrix(geno)
  class(geno) = "numeric"
  k = which(is.na(geno), arr.ind=TRUE)
  geno[k] <- rowMeans(geno, na.rm=TRUE)[k[,1]]
  as.data.frame(geno)
}




#' Function for filling in snp coordinates
#' @param geno genotype dataframe
#' @keywords constant method 
#' @export
#' @examples
#' geno_coord()
geno_coord <- function(geno){
  geno$chromosome <- sapply(strsplit(rownames(geno)[1:nrow(geno)], "_"), function(x){as.character(x[1])})
  geno$start <- sapply(strsplit(rownames(geno)[1:nrow(geno)], "_"), function(x){as.numeric(x[2])}) - 1
  start_length = as.data.frame(sapply(sapply(strsplit(rownames(geno)[1:nrow(geno)], "_"), function(x){as.character(x[3])}), nchar))
  end_length = as.data.frame(sapply(sapply(strsplit(rownames(geno)[1:nrow(geno)], "_"), function(x){as.character(x[4])}), nchar))
  df = cbind(start_length, end_length)
  colnames(df)[1] = "start_length"
  colnames(df)[2] = "end_length"
  df$max = apply(df, 1, max)
  rm(start_length, end_length)
  geno$end <- geno$start + df$max
  rm(df)
  geno
}


#' Function for mapping snps to genes in 3D method (e.g. Domain, Loop, TAD)
#' @param geno genotype dataframe
#' @param expression expression dataframe
#' @param chr chromosome number
#' @param window +/- (window) kb to consider around each gene
#' @keywords constant method 
#' @export
#' @examples
#' constant_process()
constant_process <- function(geno, expression, chr, window){
  snp_range = makeGRangesFromDataFrame(geno, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
  gene_range = makeGRangesFromDataFrame(expression, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
  gene_range = resize(gene_range, width(gene_range) + window*1000, fix = 'start')
  gene_range = resize(gene_range, width(gene_range) + window*1000, fix = 'end')
  range = as.data.frame(ranges(gene_range))
  name = as.data.frame(cbind(range$names, paste(chr, "_", range$start, "_", range$end, sep="")))
  colnames(name) = c("gene_id", "origin_id")

  fo = findOverlaps(query=snp_range, subject=gene_range, type = "within")
  df = as.data.frame(cbind(rownames(expression[subjectHits(fo),]), rownames(geno[queryHits(fo),])))
  colnames(df) = c("gene_id", "snp_id")
  df$gene_id = sapply(strsplit(as.character(df$gene_id), "\\."), function(x){as.character(x[1])})
  df$snp_id = sapply(strsplit(as.character(df$snp_id), "\\."), function(x){as.character(x[1])})
  df = df %>% group_by(gene_id) %>% mutate(snp_list = paste0(snp_id, collapse = ";")) %>% dplyr::select(gene_id, snp_list) %>% distinct(gene_id, .keep_all =TRUE)

  merge(x = name, y = df, by.x = "gene_id", by.y = "gene_id") %>% dplyr::select(origin_id, gene_id, snp_list)

}



#' Function for mapping snps to genes in 3D method (e.g. Domain, Loop, TAD)
#' @param geno genotype dataframe
#' @param expression expression dataframe
#' @param func 3D map dataframe
#' @keywords 3D method 
#' @export
#' @examples
#' func_process()
func_process <- function(geno, expression, func){
  snp_range = makeGRangesFromDataFrame(geno, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
  gene_range = makeGRangesFromDataFrame(expression, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
  func_range = makeGRangesFromDataFrame(func, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)

  fo_func_snp = findOverlaps(query=snp_range, subject=func_range, type = "within")
  df_func_snp = as.data.frame(cbind( rownames(func[subjectHits(fo_func_snp),]), rownames(geno[queryHits(fo_func_snp),]) ))
  colnames(df_func_snp) = c("func_id", "snp_id")
  df_func_snp$func_id = sapply(strsplit(as.character(df_func_snp$func_id), "\\."), function(x){as.character(x[1])})
  df_func_snp$snp_id = sapply(strsplit(as.character(df_func_snp$snp_id), "\\."), function(x){as.character(x[1])})
  df_func_snp = df_func_snp %>% dplyr::group_by(func_id) %>% dplyr::mutate(snp_list = paste0(snp_id, collapse = ";")) %>% dplyr::distinct(func_id, .keep_all = TRUE) %>% dplyr::select(func_id, snp_list)

  fo_func_gene = findOverlaps(query=gene_range, subject=func_range, type = "within")
  df_func_gene = as.data.frame(cbind( rownames(func[subjectHits(fo_func_gene),]), rownames(expression[queryHits(fo_func_gene),]) ))
  colnames(df_func_gene) = c("func_id", "gene_id")
  df_func_gene$func_id = sapply(strsplit(as.character(df_func_gene$func_id), "\\."), function(x){as.character(x[1])})
  df_func_gene$gene_id = sapply(strsplit(as.character(df_func_gene$gene_id), "\\."), function(x){as.character(x[1])})
  df_func_gene = df_func_gene %>% dplyr::group_by(func_id) %>% dplyr::mutate(gene_list = paste0(gene_id, collapse = ";")) %>% dplyr::distinct(func_id, .keep_all = TRUE) %>% dplyr::select(func_id, gene_list)

  df_gene_snp = merge(x = df_func_gene, y = df_func_snp, by.x = "func_id", by.y = "func_id")
  s = strsplit(df_gene_snp$gene_list, split = ";")

  unique(data.frame(origin_id = rep(df_gene_snp$func_id, sapply(s, length)), gene_id = unlist(s), snp_list = rep(df_gene_snp$snp_list, sapply(s, length)) ))

}


#' Function for mapping snps to genes in 3D method pchic
#' @param geno genotype dataframe
#' @param expression expression dataframe
#' @param func 3D map dataframe
#' @keywords 3D method pchic
#' @export
#' @examples
#' pchic_process()
pchic_process <- function(geno, expression, func){
  snp_range = makeGRangesFromDataFrame(geno, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
  func_range = makeGRangesFromDataFrame(func, seqnames.field=c("chromosome"),start.field="start", end.field=c("end"), ignore.strand=TRUE)
  fo_func_snp = findOverlaps(query=snp_range, subject=func_range, type = "within")
  df_gene_snp = as.data.frame(cbind( func[subjectHits(fo_func_snp),], rownames(geno[queryHits(fo_func_snp),]) )) %>% dplyr::select(-c("chromosome", "start", "end"))
  colnames(df_gene_snp) = c("origin_id", "gene_id", "snp_id")
  df_gene_snp$snp_id = sapply(strsplit(as.character(df_gene_snp$snp_id), "\\."), function(x){as.character(x[1])})
  df_gene_snp$gene_id = sapply(strsplit(as.character(df_gene_snp$gene_id), "\\."), function(x){as.character(x[1])})
  df_gene_snp = df_gene_snp %>% group_by(gene_id) %>% mutate(snp_list = paste0(snp_id, collapse = ";")) %>% distinct(gene_id, .keep_all = TRUE) %>% dplyr::select(origin_id, gene_id, snp_list)
  subset(df_gene_snp, gene_id %in% rownames(expression))
}

#' Function for scaling genotype dataframe to mean of zero and sd of one
#' @param geno_subset input genotype matrix
#' @param sample_tuning index of tunning sample
#' @param sample_testing index of testing sample
#' @keywords Scaling genotype
#' @export
#' @examples
#' geno_scale()
geno_scale <- function(geno_subset, sample_tuning, sample_testing){
  ##Normalized genotype according to the tuning set##
  geno_tuning = geno_subset %>% dplyr::select(all_of(sample_tuning))
  geno_tuning = geno_fill(geno_tuning)
  center_x_tune = rowMeans(geno_tuning, na.rm = T)
  sd_x_tune = apply(geno_tuning, 1, sd, na.rm = T)
  geno_tuning = as.data.frame(t(scale(t(geno_tuning), center = center_x_tune, scale = sd_x_tune)))
  geno_tuning = geno_tuning[which(!(rowSums(is.na(geno_tuning)) == ncol(geno_tuning))),]

  geno_testing = geno_subset %>% dplyr::select(all_of(sample_testing))
  geno_testing = geno_testing[rownames(geno_testing) %in% rownames(geno_tuning),]
  center_x_tune = center_x_tune[rownames(geno_tuning)]
  sd_x_tune = sd_x_tune[rownames(geno_tuning)]
  geno_testing = as.data.frame(t(scale(t(data.matrix(geno_testing)), center = center_x_tune, scale = sd_x_tune)))
  geno_testing[is.na(geno_testing)] <- 0

  list(geno_tuning, geno_testing)
}

#' Function for scaling expression dataframe to mean of zero
#' @param expression_subset input expression
#' @param sample_tuning index of tunning sample
#' @param sample_testing index of testing sample
#' @keywords Scaling expression
#' @export
#' @examples
#' expression_scale()
expression_scale <- function(expression_subset, sample_tuning, sample_testing){
  ##Normalize expression according to the training set##
  expression_tuning = expression_subset %>% dplyr::select(all_of(sample_tuning))
  center_y_tune = rowMeans(expression_tuning, na.rm = T)
  expression_tuning = as.data.frame(t(scale(t(expression_tuning), center = center_y_tune, scale = FALSE)))

  expression_testing = expression_subset %>% dplyr::select(all_of(sample_testing))
  expression_testing = as.data.frame(t(scale(t(expression_testing), center = center_y_tune, scale = FALSE)))

  list(expression_tuning, expression_testing)
}



#' Cross-validation
#'
#' @param geno_tuning tuning genotype set
#' @param expression_tuning tuning expression set
#' @param geno_testing testing genotype set
#' @param expression_testing testing expression set
#' @param penalty_k penalty factor
#' @param lambda_list list of l1 penalty value
#' @keywords cross-validation
#' @export
#' @examples
#' elnet()
elnet <- function(geno_tuning, expression_tuning, geno_testing, expression_testing, penalty_k, lambda_list) {

  ##Fit tuning set and predict in testing sample
  testing <- tryCatch({
    alpha.fit <- glmnet::glmnet(t(geno_tuning), t(expression_tuning), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, lambda = lambda_list)
    alpha.predicted_testing <- predict(alpha.fit, newx = t(geno_testing))
    alpha.predicted_testing
  },
  error = function(cond) {
    alpha.predicted_testing <- as.matrix(rep(mean(unlist(expression_tuning)), ncol(geno_testing)))
    alpha.predicted_testing
  })

  rownames(testing) = colnames(geno_testing)
  colnames(testing) = 1:length(lambda_list)
  testing
}

#' Function for calculating coefficient of determination#' Read plink file through R
#'
#' @param y observation vector
#' @param y_pred prediction vector
#' @keywords calcaulate R^2
#' @export
#' @examples
#' calc_R2()
calc_R2 <- function(y, y_pred) {
  tss <- sum((y-mean(y))**2)
  rss <- sum((y - y_pred)**2)
  1 - rss/tss
}





#' Read plink file through R
#'
#' @param root plink file
#' @param impute How to impute missing value
#' @keywords read plink
#' @export
#' @examples
#' read_plink_custom()
read_plink_custom <- function(root, impute = c('none', 'avg', 'random')) {
  if(impute == 'random') {
    stop("The 'impute' random option has not been implemented.", call. = FALSE)
  }

  ## structure from https://github.com/gabraham/plink2R/blob/master/plink2R/R/plink2R.R
  proot <- path.expand(root)

  bedfile <- paste(proot, ".bed", sep="")
  famfile <- paste(proot, ".fam", sep="")
  bimfile <- paste(proot, ".bim", sep="")

  ## Could change this code to use data.table
  bim <- read.table(bimfile, header=FALSE, sep="", stringsAsFactors=FALSE)
  fam <- read.table(famfile, header=FALSE, sep="", stringsAsFactors=FALSE)
  ## Set the dimensions
  geno <- BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))

  ## Convert to a matrix
  geno <- as.matrix(geno)
  if(impute == 'avg') {
    ## Check if any are missing
    geno_na <- is.na(geno)
    if(any(geno_na)) {
      means <- colMeans(geno, na.rm = TRUE)
      geno[geno_na] <- rep(means, colSums(geno_na))
    }
  }
  colnames(geno) <- bim[,2]
  rownames(geno) <- paste(fam[,1], fam[, 2], sep=":")

  list(bed=geno, fam=fam, bim=bim)
}




#' Create gene expression prediction model using window type and penalty factor derived from the CV step
#'
#' @param geno_path genotype path, plink format
#' @param exp_path expression path, txt format
#' @param encode_path epigenomic annotation path, txt format
#' @param output_path dir of your output file
#' @param window_path_list A list of 3d genomic winodows you need
#' @keywords weights
#' @export
#' @examples
#' PUMICE.compute_weight()
PUMICE.compute_weight<-function(geno_path,exp_path,encode_path,output_path,window_path_list){

opt = list()

opt$geno = geno_path
opt$exp = exp_path
opt$out = output_path
opt$encode_path = encode_path

opt$Loop_path = window_path_list[["Loop"]]
opt$TAD_path = window_path_list[["TAD"]]
opt$Domain_path = window_path_list[["Domain"]]
opt$pcHiC_path = window_path_list[["pcHiC"]]

opt$fold = 5

#########################################################################################################################
##Processing genotype##
#########################################################################################################################
cat( "UPDATE: Processing genotype data\n" , file=stderr())


genoc = read_plink_custom(opt$geno, impute = 'avg')
geno = genoc$bed
geno = as.data.frame( t(geno) )
opt$chr<-genoc$bim$V1[1]

colnames(geno) = sapply(strsplit(colnames(geno), ":"), function(x){as.character(x[1])})
pos<-sapply(strsplit(rownames(geno), "_"), function(x){as.character(x[2])})


##Create subjectID-sample ID conversion dataframe for cross-validation step##
sample_conversion = as.data.frame(cbind( "sample_ID" = seq(1, ncol(geno)), "subject_ID" = colnames(geno)))
names(geno) = sample_conversion$sample_ID[match(names(geno), sample_conversion$subject_ID)]

# geno = mutate_all(geno, function(x) as.numeric(as.character(x)))

#head(geno[,1:10])


#########################################################################################################################
##Processing expression##
#########################################################################################################################

cat( "UPDATE: Processing expression data\n" , file=stderr())
expression = as.data.frame(fread(opt$exp, header = TRUE))
ID_exp<-(colnames(expression[,5:ncol(expression)]))
colnames(expression)<-c(c("gene_id","chromosome","start","end"),ID_exp)
index_exp<-which(match(ID_exp,sample_conversion$subject_ID)>0)
index<-c(c(1,2,3,4),index_exp+4)
expression = expression[,index]
#length(index_exp_male)
#length(ID_male)
sample_overlap = intersect(sample_conversion$subject_ID, ID_exp)
expression = expression %>% dplyr::select(gene_id, chromosome, start, end, all_of(sample_overlap))
names(expression)[5:ncol(expression)] = sample_conversion$sample_ID[match(names(expression)[5:ncol(expression)], sample_conversion$subject_ID)]

#########################################################################################################################
##Processing epigenomic/3D genomics data##
#########################################################################################################################
##Find sample overlap between genotype and expression data##
sample_overlap = intersect(colnames(geno), colnames(expression))

##Update the sample list in both genotype and expression data##
geno = geno %>% dplyr::select(all_of(sample_overlap))
expression = expression %>% dplyr::select(gene_id, chromosome, start, end, all_of(sample_overlap)) %>% remove_rownames %>% column_to_rownames(var="gene_id")
expression = expression[which(expression$chromosome==opt$chr),]



##Filling in SNP coordinates##
geno = geno_coord(geno)

##Filter out indel/cnv##
geno = geno %>% filter(end - start == 1)



#########################################################################################################################
##Import result from nested cross-validation##
#########################################################################################################################
##Import result from nested cv##

filename = dir(opt$out, pattern = "^result_cv")

filename = filename[grep(paste0("chr",opt$chr,".txt"),filename)]


pumice_int_total = data.frame()
for (i in 1:length(filename)){
  print(i)
  try({
    pumice_int_temp_total = read.table(paste(opt$out, "/", filename[i], sep = ""), header = TRUE)
    pumice_int_total = rbind(pumice_int_total, pumice_int_temp_total)
  })
}
pumice_int_total = as.data.frame(pumice_int_total %>% filter(penalty != "NA") %>% group_by(gene_id) %>% dplyr::slice(which.min(cvm_v)))
nrow(pumice_int_total)
iv_path = paste(opt$out, paste0("/total_cv_pumice_chr",opt$chr,".txt"), sep = "")
write.table(pumice_int_total, iv_path, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



geno_bed = data.frame("chr" = paste("chr", geno$chromosome, sep = ""), "start" = geno$start, "end" = geno$end, "snpid" = rownames(geno))


##Import ENCODE file##
penalties <- tryCatch({

  encode <- as.data.frame(fread(opt$encode_path, header = T))
  encode<-encode[which(encode$overlap==1),]

  geno_bed$overlap <- 0
  geno_bed[geno_bed$snpid %in% encode$snp, "overlap"] <- 1
  geno_bed <- geno_bed[match(rownames(geno), geno_bed$snpid ),]


  ##Penalty1 is based on cut off of number of annotation and no penalization on established predictors##
  penalty1 <- c(rep(1, nrow(geno)))
  penalty1[which(geno_bed$overlap > 0)] <- 0

  ##Penalty2 is based on cut off of number of annotation and set fixed penalty at 1/6##
  penalty2 <- c(rep(1, nrow(geno)))
  penalty2[which(geno_bed$overlap > 0)] <- 1/6

  ##Penalty3 is based on cut off of number of annotation and set fixed penalty at 2/6##
  penalty3 <- c(rep(1, nrow(geno)))
  penalty3[which(geno_bed$overlap > 0)] <- 2/6

  ##Penalty4 is based on cut off of number of annotation and set fixed penalty at 3/6##
  penalty4 <- c(rep(1, nrow(geno)))
  penalty4[which(geno_bed$overlap > 0)] <- 3/6

  ##Penalty5 is based on cut off of number of annotation and set fixed penalty at 4/6##
  penalty5 <- c(rep(1, nrow(geno)))
  penalty5[which(geno_bed$overlap > 0)] <- 4/6

  ##Penalty6 is based on cut off of number of annotation and set fixed penalty at 5/6##
  penalty6 <- c(rep(1, nrow(geno)))
  penalty6[which(geno_bed$overlap > 0)] <- 5/6

  ##Penalty7 is elastic net##
  penalty7 <- c(rep(1, nrow(geno)))
  penalty7[which(geno_bed$overlap > 0)] <- 1

  rbind(penalty1, penalty2, penalty3, penalty4, penalty5, penalty6, penalty7)
},
error = function(cond) {
  t(as.matrix(c(rep(1,nrow(geno)))))
})
colnames(penalties) = rownames(geno)
n_penalties = nrow(penalties)

#########################################################################################################################
##Create cross-validation folds##
#########################################################################################################################
set.seed(123)
cv_fold_testing = createFolds(sample_conversion$sample_ID, k = opt$fold, list = TRUE)

set.seed(123)
cv_fold_tuning = createFolds(sample_conversion$sample_ID, k = opt$fold, returnTrain = TRUE)

#########################################################################################################################
##Running PUMICE##
#########################################################################################################################
cat( "UPDATE: Start running cross-validation to create gene expression prediction models\n" , file=stderr())

coef = data.frame()
summary = data.frame()
for (i in 1:nrow(pumice_int_total)) {
  tryCatch({
    func_type = pumice_int_total$type[i]
    expression_subset = expression[rownames(expression) %in% pumice_int_total$gene_id[i],]

    cat( "UPDATE: Start running gene", rownames(expression_subset), "\n" , file=stderr())

    if (func_type == "pcHiC") {
      func_path = opt$pcHiC_path
      func = read.table(func_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
      func = func %>% filter(ensembl_gene_id == rownames(expression_subset)) %>% dplyr::select(-c("length", "MinusLog10Pval"))
      map = pchic_process(geno, expression_subset, func)
    } else if (func_type == "Loop") {
      func_path = opt$Loop_path
      func = read.table(func_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
      func$func_id = paste(func$chromosome, func$start, func$end, sep="_")
      func = func %>% filter(func_id == pumice_int_total$region[i]) %>% remove_rownames %>% column_to_rownames(var="func_id")
      map = func_process(geno, expression_subset, func)
    } else if (func_type == "TAD") {
      func_path = opt$TAD_path
      func = read.table(func_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
      func$func_id = paste(func$chromosome, func$start, func$end, sep="_")
      func = func %>% filter(func_id == pumice_int_total$region[i]) %>% remove_rownames %>% column_to_rownames(var="func_id")
      map = func_process(geno, expression_subset, func)
    } else if (func_type == "Domain") {
      func_path = opt$Domain_path
      func = read.table(func_path, header = TRUE, na.string ='.', as.is=TRUE,check.names=FALSE)
      func$func_id = paste(func$chromosome, func$start, func$end, sep="_")
      func = func %>% filter(func_id == pumice_int_total$region[i]) %>% remove_rownames %>% column_to_rownames(var="func_id")
      map = func_process(geno, expression_subset, func)
    } else {
      func_type = as.numeric(as.character(func_type))
      map = constant_process(geno, expression_subset, opt$chr, func_type)
    }

    geno_subset = geno[unlist(strsplit(as.character(map[1,3]), ";")),] %>% dplyr::select(-c("chromosome", "start", "end"))
    expression_subset = expression_subset %>% dplyr::select(-c("chromosome", "start", "end"))

    expression_subset<-t(scale(t(expression_subset)))
    expression_subset<-as.data.frame(expression_subset)

    origin_id = as.character(map[1,1])

    ##Need to create a list of predefined lambdas##
    ##Standardized genotype matrix##
    geno_subset_all = geno_fill(geno_subset)
    center_x_train_all = rowMeans(geno_subset_all, na.rm = T)
    sd_x_train_all = apply(geno_subset_all, 1, sd, na.rm = T)
    geno_subset_all = as.data.frame(t(scale(t(geno_subset_all), center = center_x_train_all, scale = sd_x_train_all)))
    geno_subset_all = geno_subset_all[which(!(rowSums(is.na(geno_subset_all)) == ncol(geno_subset_all))),]

    ##Center expression matrix##
    center_y_train = rowMeans(expression_subset, na.rm = T)
    expression_subset_all = as.data.frame(t(scale(t(expression_subset), center = center_y_train, scale = FALSE)))

    ##Retrieve penalty coefficients##
    penalty_k = penalties[pumice_int_total$penalty[i], rownames(geno_subset_all)]

    ##Find a list of lambda##
    lambda_list =  (glmnet::glmnet(t(geno_subset_all), t(expression_subset_all), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F))$lambda

    prediction_testing = data.frame()
    for (f in 1:opt$fold) {
      print(f)
      ##Identify individual in each folds##
      sample_testing = intersect(sample_overlap, cv_fold_testing[[f]])
      sample_tuning = intersect(sample_overlap, cv_fold_tuning[[f]])

      ##Scale the predictors accordingly##
      geno_output = geno_scale(geno_subset, as.character(sample_tuning), as.character(sample_testing))
      geno_tuning = geno_output[[1]]
      geno_testing = geno_output[[2]]

      ##Scale the expression accordingly##
      expression_output = expression_scale(expression_subset, as.character(sample_tuning), as.character(sample_testing))
      expression_tuning = expression_output[[1]]
      expression_testing = expression_output[[2]]

      if (nrow(geno_testing) > 1) {
        alpha.predicted_testing = elnet(geno_tuning, expression_tuning, geno_testing, expression_testing, penalty_k, lambda_list)
        prediction_testing = rbind(prediction_testing, cbind( t(expression_testing), alpha.predicted_testing ) )
      }
    }
    prediction_testing = prediction_testing[match( colnames(expression_subset), rownames(prediction_testing) ),]

    ##Identify best lambda value##
    cvm_testing = sapply(2:ncol(prediction_testing), function(x) { mean( (unlist( prediction_testing[,x] ) - unlist( prediction_testing[,1] ))^2 ) } )
    best_lambda = lambda_list[which.min(cvm_testing)]

    ##Create final model##
    alpha.fit = glmnet::glmnet(t(geno_subset_all), t(expression_subset_all), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, lambda = best_lambda)

    if ( all(alpha.fit$beta == 0) ) {
      summary_temp = data.frame( "gene_id" = rownames(expression_subset_all),
                                 "func_type" = func_type, "region" = origin_id, fit = "elnet",
                                 "alpha" = 0.5, "penalty" = pumice_int_total$penalty[i],
                                 "input_snps" = nrow(geno_subset_all), "nonzero_snps" = 0)
      summary = rbind(summary, summary_temp)
    } else {
      coef_temp = as.matrix(alpha.fit$beta[which(alpha.fit$beta[,1] != 0),])
      rownames(coef_temp) = rownames(coef(alpha.fit, s = 'lambda.min'))[coef(alpha.fit, s = 'lambda.min')[,1]!= 0][-1]

      coef_output = data.frame( "gene_id" = rep(rownames(expression_subset_all), nrow(coef_temp)),
                                "snp" = rownames(coef_temp), "weight" = coef_temp[,1]) %>% remove_rownames
      coef = rbind(coef, coef_output)

      summary_temp = data.frame( "gene_id" = rownames(expression_subset_all),
                                 "func_type" = func_type, "region" = origin_id, fit = "elnet",
                                 "alpha" = 0.5, "penalty" = pumice_int_total$penalty[i],
                                 "input_snps" = nrow(geno_subset_all), "nonzero_snps" = length(alpha.fit$beta[which(alpha.fit$beta[,1] != 0),]))
      summary = rbind(summary, summary_temp)
    }
  }, error=function(e){})

  filename_summ_result <- paste(opt$out, "/total_summary_pumice_chr",opt$chr,".txt", sep="")
  write.table(summary, filename_summ_result, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

  filename_coef_result <- paste(opt$out, "/total_coef_pumice_chr",opt$chr,".txt", sep="")
  write.table(coef, filename_coef_result, quote=FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

}

return(list("cv_list"=pumice_int_total,"summary"=summary,"coef"=coef))
}
