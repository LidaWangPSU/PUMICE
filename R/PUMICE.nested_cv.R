#' Nested Cross-validation
#'
#' @param geno_traning training genotype set
#' @param expression_training training expression set
#' @param fold_training index for traning for fold
#' @param geno_validating validating genotype set
#' @param expression_validating validating expression set
#' @param geno_tuning tuning genotype set
#' @param expression_tuning tuning expression set
#' @param fold_tuning index for tuning for fold
#' @param geno_testing testing genotype set
#' @param expression_testing testing expression set
#' @param penalty_k penalty factor
#' @keywords Nested cross-validation
#' @export
#' @examples
#' elnet.cv()
elnet.cv <- function(geno_training, expression_training, fold_training, geno_validating, expression_validating, geno_tuning, expression_tuning, fold_tuning, geno_testing, expression_testing, penalty_k) {

  ##Fit best lambda in training set and predict in validating sample
  alpha.predicted_validating <- tryCatch({
    alpha.fit <- glmnet::cv.glmnet(t(geno_training), t(expression_training), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, nfolds = opt$fold, type.measure='mse', foldid = fold_training)
    predict(alpha.fit, newx = t(geno_validating), s = "lambda.min")
  },
  error = function(cond) {
    as.matrix(rep(mean(unlist(expression_training)), ncol(geno_validating)))
  })
  rownames(alpha.predicted_validating) = colnames(geno_validating)

  ##Fit tuning set and predict in testing sample
  testing <- tryCatch({
    alpha.fit <- glmnet::cv.glmnet(t(geno_tuning), t(expression_tuning), alpha = 0.5, family="gaussian", penalty.factor = penalty_k, standardize = F, nfolds = opt$fold, type.measure='mse', foldid = fold_tuning)
    alpha.predicted_testing <- predict(alpha.fit, newx = t(geno_testing), s = "lambda.min")
    alpha_coef = as.matrix(alpha.fit$glmnet.fit$beta[, which.min(alpha.fit$cvm)])
    list(alpha.predicted_testing, alpha_coef)
  },
  error = function(cond) {
    alpha.predicted_testing <- as.matrix(rep(mean(unlist(expression_tuning)), ncol(geno_testing)))
    alpha_coef = as.matrix(rep(0, nrow(geno_tuning)))
    list(alpha.predicted_testing, alpha_coef)
  })

  alpha.predicted_testing = testing[[1]]
  alpha_coef = testing[[2]]
  rownames(alpha.predicted_testing) = colnames(geno_testing)
  colnames(alpha.predicted_testing) = "test"
  rownames(alpha_coef) = rownames(geno_tuning)
  colnames(alpha_coef) = "weight"

  ##Calculate encode importance##
  enc_imp <- tryCatch({
    coef_nonzero = as.matrix(abs(alpha_coef[which(alpha_coef[,1] != 0),]))
    encode_snp = geno_bed %>% filter(overlap > 0) %>% dplyr::select(snpid)
    coef_nonzero_encode = as.matrix(coef_nonzero[rownames(coef_nonzero) %in% encode_snp$snpid,])
    ifelse( !(is.na((sum(coef_nonzero_encode))/(sum(coef_nonzero)))), (sum(coef_nonzero_encode))/(sum(coef_nonzero)), 0)
  },
  error = function(cond) {
    return(0)
  })

  ##Calculate Pearson correlation, zscore, and pval of the validating fold##
  R2_v = calc_R2(unlist(expression_validating), unlist(alpha.predicted_validating))
  pear_folds_v = ifelse(sd(alpha.predicted_validating) != 0, cor( unlist(alpha.predicted_validating), unlist(expression_validating), method = "pearson"), 0)
  spear_folds_v = ifelse(sd(alpha.predicted_validating) != 0, cor( unlist(alpha.predicted_validating), unlist(expression_validating), method = "spearman"), 0)

  ##Calculate Pearson correlation, zscore, and pval of the test fold##
  R2_t = calc_R2(unlist(expression_testing), unlist(alpha.predicted_testing))
  pear_folds_t = ifelse(sd(alpha.predicted_testing) != 0, cor( unlist(alpha.predicted_testing), unlist(expression_testing), method = "pearson"), 0)
  spear_folds_t = ifelse(sd(alpha.predicted_testing) != 0, cor( unlist(alpha.predicted_testing), unlist(expression_testing), method = "spearman"), 0)
  pear_zscore_folds <- atanh(pear_folds_t)*sqrt(length(expression_testing) - 3) #Fisher transformation
  spear_zscore_folds <- atanh(spear_folds_t)*sqrt(length(expression_testing) - 3) #Fisher transformation

  list(alpha.predicted_validating, R2_v, pear_folds_v, spear_folds_v, enc_imp, R2_t, pear_folds_t, spear_folds_t, pear_zscore_folds, spear_zscore_folds)
}




#' Run nested cross-validation to determine which window type and penalty factor are optimal 
#'
#' @param geno_path genotype path, plink format
#' @param exp_path expression path, txt format
#' @param encode_path epigenomic annotation path, txt format
#' @param output_path dir of your output file
#' @param window specify the constant or 3D window you need
#' @param window_path_list A list of 3d genomic winodows you need
#' @keywords nested.cv
#' @export
#' @examples
#' PUMICE.nested_cv()
PUMICE.nested_cv<-function(geno_path,exp_path,encode_path,output_path,window,window_path_list=NULL){

opt = list()
opt$geno = geno_path
opt$exp = exp_path
opt$out = output_path
opt$encode_path = encode_path
opt$fold = 5


if(window=="1000"|window=="250"){
  opt$method = "constant"
  opt$type = window
  opt$window_path=NULL
}else if(window=="TAD"|window=="Loop"|window=="pcHiC"|window=="Domain"){
  opt$method = "3d"
  opt$type = window
  opt$window_path=window_path_list[[window]]
}else{
  error=paste0("Error: there is no window called ",window,", Please check the input window name.")
  return(error)
}







#########################################################################################################################
##Processing genotype##
#########################################################################################################################

genoc = read_plink_custom(opt$geno, impute = 'avg')
geno = genoc$bed
geno = as.data.frame( t(geno) )

opt$chr<-genoc$bim$V1[1]

colnames(geno) = sapply(strsplit(colnames(geno), ":"), function(x){as.character(x[1])})
pos<-sapply(strsplit(rownames(geno), "_"), function(x){as.character(x[2])})





##Create subjectID-sample ID conversion dataframe for cross-validation step##
sample_conversion = as.data.frame(cbind( "sample_ID" = seq(1, ncol(geno)), "subject_ID" = colnames(geno)))
names(geno) = sample_conversion$sample_ID[match(names(geno), sample_conversion$subject_ID)]



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


##Create map file##
cat( "UPDATE: Processing window", opt$method, "-", opt$type, "\n" , file=stderr())

if (opt$method == "3d"){
  if (opt$type == "pcHiC"){
    func = as.data.frame(read.table(opt$window_path, header = TRUE))
    func<-unique(func)
    func = func %>% dplyr::select(-c("length", "MinusLog10Pval"))
    map = pchic_process(geno, expression, func)
  } else if (opt$type == "Loop") {
    func = as.data.frame(fread(opt$window_path, header = TRUE))
    func<-unique(func)
    func$func_id = paste(func$chromosome, func$start, func$end, sep="_")
    func = func[which(func$chromosome==opt$chr),]
    func = func %>% remove_rownames %>% column_to_rownames(var="func_id")
    map = func_process(geno, expression, func)
  } else if (opt$type == "TAD") {
    func = as.data.frame(fread(opt$window_path, header = TRUE))
    func<-unique(func)
    func$func_id = paste(func$chromosome, func$start, func$end, sep="_")
    func = func %>% remove_rownames %>% column_to_rownames(var="func_id")
    func = func[which(func$chromosome==opt$chr),]
    map = func_process(geno, expression, func)
  }else {
    func = as.data.frame(fread(opt$window_path, header = TRUE))
    func<-unique(func)
    func$func_id = paste(func$chromosome, func$start, func$end, sep="_")
    func = func %>% remove_rownames %>% column_to_rownames(var="func_id")
    func = func[which(func$chromosome==opt$chr),]
    map = func_process(geno, expression, func)
  }

} else if (opt$method == "constant"){
  func_type = as.numeric(as.character(opt$type))
  map = constant_process(geno, expression, opt$chr, func_type)
}

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

  ##Penalty3 is based on cut off of number of annotation and set fixed penalty at	2/6##
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

##Update the gene list in expression dataframe according to map file##
geno = geno %>% dplyr::select(-c("chromosome", "start", "end"))
expression = expression[rownames(expression) %in% map$gene_id,] %>% dplyr::select(-c("chromosome", "start", "end"))
expression<-t(scale(t(expression)))
expression<-as.data.frame(expression)
#########################################################################################################################
##Create cross-validation folds##
#########################################################################################################################
set.seed(123)
cv_fold_testing = createFolds(sample_conversion$sample_ID, k = opt$fold, list = TRUE)

set.seed(123)
cv_fold_tuning = createFolds(sample_conversion$sample_ID, k = opt$fold, returnTrain = TRUE)

cv_fold_validating = list()
for (f in 1:opt$fold){
  cv_fold_validating[[f]] = cv_fold_testing[[(f%%opt$fold)+1]]
}

cv_fold_training = list()
for (f in 1:opt$fold){
  cv_fold_training[[f]] = setdiff(cv_fold_tuning[[f]], cv_fold_validating[[f]])
}

#########################################################################################################################
###Perform nested cross-validation##
#########################################################################################################################
cat( "UPDATE: Start running nested cross-validation\n" , file=stderr())

summ_result = data.frame()
for (i in 1:nrow(map)) {
  geno_subset = geno[unique(unlist(strsplit(as.character(map[i,3]), ";"))),]

  if (nrow(geno_subset) > 1) {

    expression_subset = expression[as.character(map[i,2]),]
    origin_id = as.character(map[i,1])

    cat( "UPDATE: Start running gene", rownames(expression_subset), "\n" , file=stderr())
    for (k in 1:nrow(penalties)){
      cat( "UPDATE: Start running gene", rownames(expression_subset), rownames(penalties)[k], "\n" , file=stderr())
      prediction_validating = data.frame()
      prediction_testing = data.frame()
      R2_folds_v = rep(0,opt$fold)
      pear_folds_v = rep(0,opt$fold)
      spear_folds_v = rep(0,opt$fold)
      R2_folds_t = rep(0,opt$fold)
      pear_folds_t = rep(0,opt$fold)
      spear_folds_t = rep(0,opt$fold)
      pear_zscore_folds = rep(0,opt$fold)
      spear_zscore_folds = rep(0,opt$fold)
      enc_imp = rep(0,opt$fold)
      penalty_k = penalties[k, rownames(geno_subset)]

      for (f in 1:opt$fold) {

        ##Identify individual in each folds##
        sample_testing = intersect(sample_overlap, cv_fold_testing[[f]])
        sample_tuning = intersect(sample_overlap, cv_fold_tuning[[f]])
        sample_validating = intersect(sample_overlap, cv_fold_validating[[f]])
        sample_training = intersect(sample_overlap, cv_fold_training[[f]])

        ##Ceate fold-ids for inner cv##
        set.seed(123)
        fold_tuning = sample(rep(seq(opt$fold), length.out = length(cv_fold_tuning[[f]])))
        fold_tuning = fold_tuning[match(sample_tuning, cv_fold_tuning[[f]])]
        set.seed(123)
        fold_training = sample(rep(seq(opt$fold), length.out = length(cv_fold_training[[f]])))
        fold_training = fold_training[match(sample_training, cv_fold_training[[f]])]

        ##Scale the predictors accordingly##
        geno_output = geno_scale(geno_subset, as.character(sample_training), as.character(sample_validating))
        geno_training = geno_output[[1]]
        geno_validating = geno_output[[2]]

        geno_output = geno_scale(geno_subset, as.character(sample_tuning), as.character(sample_testing))
        geno_tuning = geno_output[[1]]
        geno_testing = geno_output[[2]]

        ##Scale the expression accordingly##
        expression_output = expression_scale(expression_subset, as.character(sample_training), as.character(sample_validating))
        expression_training = expression_output[[1]]
        expression_validating = expression_output[[2]]

        expression_output = expression_scale(expression_subset, as.character(sample_tuning), as.character(sample_testing))
        expression_tuning = expression_output[[1]]
        expression_testing = expression_output[[2]]

        result = elnet.cv(geno_training, expression_training, fold_training, geno_validating, expression_validating, geno_tuning, expression_tuning, fold_tuning, geno_testing, expression_testing, penalty_k)

        alpha.predicted_validating = result[[1]]
        colnames(alpha.predicted_validating) = "predicted"
        prediction_validating = rbind(prediction_validating, cbind( t(expression_validating), alpha.predicted_validating ) )
        R2_folds_v[f] = result[[2]]
        pear_folds_v[f] = result[[3]]
        spear_folds_v[f] = result[[4]]
        enc_imp[f] = result[[5]]
        R2_folds_t[f] =	result[[6]]
        pear_folds_t[f] = result[[7]]
        spear_folds_t[f] = result[[8]]
        pear_zscore_folds[f] = result[[9]]
        spear_zscore_folds[f] = result[[10]]
      }

      prediction_validating = prediction_validating[match( colnames(expression_subset), rownames(prediction_validating) ),]

      ##Calculate validating cvm##
      cvm_validating = min(sapply(2:ncol(prediction_validating), function(x) { mean( (unlist( prediction_validating[,x] ) - unlist( prediction_validating[,1] ))^2 ) } ))

      ##Calculate average of correlation across folds##
      R2_avg_v = mean(R2_folds_v)
      pear_avg_v = mean(pear_folds_v)
      spear_avg_v = mean(spear_folds_v)

      R2_avg_t = mean(R2_folds_t)
      pear_avg_t = mean(pear_folds_t)
      spear_avg_t = mean(spear_folds_t)

      ##Combine Z-scores via Stouffer's method##
      pear_zscore_est = sum(pear_zscore_folds) / sqrt(opt$fold)
      pear_stouffer_pval <- 2*pnorm(abs(pear_zscore_est), lower.tail = FALSE)

      spear_zscore_est = sum(spear_zscore_folds) / sqrt(opt$fold)
      spear_stouffer_pval <- 2*pnorm(abs(spear_zscore_est), lower.tail = FALSE)

      ##Calculatre average of encode importance##
      enc_avg = mean(enc_imp)

      els_output = data.frame("gene_id" = rownames(expression_subset), "type" = opt$type, "region" = origin_id,
                              "penalty" = k, "cvm_v" = cvm_validating, "R2_avg_v" = R2_avg_v, "pear_avg_v" = pear_avg_v, "spear_avg_v" = spear_avg_v,
                              "R2_avg_t" = R2_avg_t, "pear_avg_t" = pear_avg_t, "spear_avg_t" = spear_avg_t,
                              "pear_stouffer_pval" = pear_stouffer_pval, "spear_stouffer_pval" = spear_stouffer_pval,
                              "enc_avg" = enc_avg) %>% remove_rownames
      summ_result = rbind(summ_result, els_output)
      filename_summ_result <- paste(opt$out, "/result_cv_", opt$type,"_chr", opt$chr,".txt", sep="")
      write.table(summ_result, filename_summ_result, quote=FALSE, sep = "\t", row.names = FALSE)
    }

  } else {
    next
  }
}

return(summ_result)
}
