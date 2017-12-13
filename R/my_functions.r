#' Reload a package.
#' @description Unload and reload a package and sets the namespace search order.  
#' @param package Unquoted package name.
#' @param pos Namespace search position.
#' @examples
#' \code{reload(trqwe)}
#' @export
reload <- function(package, pos=2) {
    pstring <- deparse(substitute(package))
    print(pstring)
    unloadNamespace(pstring)
    library(pstring, pos=pos, character.only=T)
}

#' Unload and reload trqwe.
#' @description Unload and reload trqwe.
#' Shortcut for \code{reload(trqwe)}
#' @export
reloadtrqwe <- function() {
    reload(trqwe)
}


#' Install Bioconductor package.
#' @description Utility function for installing packages from bioconductor easily.
#' @param package unquoted package name.
#' @examples
#' \code{bioc(DESeq2)}
#' @export
bioc <- function(package) {
    source("http://bioconductor.org/biocLite.R")
    pstring <- deparse(substitute(package))
    biocLite(pstring)
}

#' Install CRAN package.
#' @description Utility function for installing packages from CRAN easily.
#' @param package unquoted package name.
#' @examples
#' \code{install(Rcpp)}
#' @export
install <- function(package, repos="http://cran.us.r-project.org") {
    pstring <- deparse(substitute(package))
    install.packages(pstring, repos=repos)
}

#' Time profiling.
#' @description Reports time in seconds to the R prompt of the previous command.
#' This function adds a callback which saves the running of individual commands and reports the time in seconds on the next line.  
#' @examples
#' > timePrompt()
#' 0.000s> x <- sample(1:10, size=1e8, replace=T)
#' 1.240s> 
#' Note - this time is not accurate if child processes or multithreading is involved.  
#' @export
timePrompt <- function() {
    removeTaskCallback(".exec_time_prompt")
    updatePrompt <- function(...) {
        utime <- proc.time()
        utime_prev <- get0(".prev_time", envir=globalenv(), ifnotfound=utime)
        utime_diff <- sprintf("%.3f", sum(utime[c(1,2)] - utime_prev[c(1,2)]))
        options(prompt=paste0(utime_diff,"s> "))
        assign(".prev_time", utime, envir=globalenv())
        return(TRUE)
        }
    ret <- addTaskCallback(updatePrompt, name=".exec_time_prompt")
}

#' Variable information.  
#' @description Automatically stores basic information of variables in the previous command.
#' This function adds a callback which reports information on previous variables and stores this information in \code{.stats}.
#' @examples
#' statsCallback()
#' my_data <- VADeaths
#' .stats
#' > [1] "dim: 5 4, length: 20, class: matrix, typeof: double"
#' @export
statsCallback <- function() {
    removeTaskCallback(".stats_callback")
    statsUpdate <- function(...) {
        lv <- get0(".Last.value", envir=globalenv(), ifnotfound=NULL)
        dim <- paste0(dim(lv), collapse=" ")
        length <- length(lv)
        class <- class(lv)
        typeof <- typeof(lv)
        assign(".stats", paste0("dim: ",dim,", length: ", length, ", class: ", class, ", typeof: ",typeof), envir=globalenv())
        return(TRUE)
        }
    ret <- addTaskCallback(statsUpdate, name=".stats_callback")
}

#' All duplicates.
#' @description Finds all duplicates in a vector including first instances of duplicates.
#' @param vec A vector.
#' @return A boolean vector of the same length as \code{vec}.  TRUE if the element is duplicated, FALSE if not.  
#' @examples
#' allDups(sample(1:100, size=100, replace=T)
#' @export
allDups <- function(vec) {
	duplicated(vec) | rev(duplicated(rev(vec)))
}

#' Fast readLines.
#' @description Replacement for readLines, faster.  
#' @param fname Filename to read.
#' @param newlinechar The new line character in the file.  
#' @return A vector containing all the lines in the file.
#' @examples
#' library(microbenchmark)
#' writeLines(replicate(100000, sample(letters, size=100, replace=T)), con="/tmp/temp.txt")
#' microbenchmark(fastReadLines("/tmp/temp.txt"), 
#'                readLines("/tmp/temp.txt"), times=3)[,c("expr", "mean"), drop=F]
#' 
#' Unit: milliseconds
#'                           expr      min       lq     mean   median       uq
#' fastReadLines("/tmp/temp.txt") 335.3874 335.9188 336.5900 336.4502 337.1912
#'     readLines("/tmp/temp.txt") 724.9136 725.4523 727.2444 725.9911 728.4098
#' @export
fastReadLines <- function(fname, newlinechar="\n", nchars=NULL) {
    if(is.null(nchars)) {
        nchars = file.info(fname)$size
    }
    buf = readChar(fname, nchars, useBytes=T)
    strsplit(buf,newlinechar,fixed=T,useBytes=T)[[1]]
}

#' Upper-left corner of matrix.  
#' @description Shortcut function for previewing a matrix or data.frame by displaying the upper-left corner.  Similar to \code{head}. 
#' @param x A wide matrix or data.frame.  
#' @param n Number of lines to display.  Default 10. 
#' @param ncols Number of columns to display.  Default 10. 
#' @return \code{n} by \code{ncols} subset of the matrix taken from the upper-left corner.  
#' @export
head2 <- function(x, n = 10, ncols = 10) {
    head(x[,seq_len(ncols)])
}

#' Lower-left corner of matrix.  
#' @description Shortcut function for previewing the bottom of a matrix or data.frame by displaying the lower-left corner.  Similar to \code{tail}. 
#' @param x A wide matrix or data.frame.  
#' @param n Number of lines to display.  Default 10. 
#' @param ncols Number of columns to display.  Default 10. 
#' @return \code{n} by \code{ncols} subset of the matrix taken from the lower-left corner.  
#' @export
tail2 <- function(x, n = 10, ncols = 10) {
    tail(x[,seq_len(ncols)])
}

#' Parse TCGA barcode.
#' @description Taking in a full TCGA sample barcode, or any subset of the barcode, and return extracted values.  
#' @param x TCGA barcode, or a vector of barcodes.  
#' @param what Which information to return.  
#' @return The specified information contained in the barcode of the same length as \code{x}. 
#' @examples
#' TCGA_barcode(c("TCGA-02-0001-01C-01D-0182-01", "TCGA-02-0001-11C-01D-0182-01"), what="tissue")
#' 
#' [1] "01" "11"
#' @export
TCGA_barcode <- function(x, what = "patient") {
    reg <- "TCGA-(..)-(....)-?(..)?(.)?-?(..)?(.)?-?(....)?-?(..)?"
    group <- switch(what, 
        TSS = "\\1",
        participant = "\\2",
        patient = "TCGA-\\1-\\2",
        sample = "\\3",
        tissue = "\\3",
        vial = "\\4",
        portion = "\\5",
        analyte = "\\6",
        plate = "\\7",
        center = "\\8")
    gsub(reg, group, x, perl=T)
}

#' Multiple grepl.
#' @description Takes in a list of regex patterns and returns true if any pattern matches. 
#' @param pattern A vector of regex patterns.  
#' @param x A vector of strings to search.  
#' @return A vector of the same length as \code{x}, TRUE if any pattern matches.  
#' @examples
#' x <- fastReadLines("http://textfiles.com/ufo/ufobooks.ufo", newlinechar="\r\n", 1e5)
#' head(x[mgrepl(c("ALIENS", "UFO"), x)])
#' 
#' [1] "                         UFOLOGY BOOKS (REVISION 2.1 343 books)"
#' [2] "     H. S. Stewart on the  subject of UFO's. The list is alphabetic"
#' [3] "     Tom Mickus's most excellent board UFONET I.  (416-237-1204)"
#' [4] "     Bill Adler               *   LETTERS TO THE AIR FORCE ON UFOS  1967"
#' [5] "     Gordon W. Allen              OVERLORDS OLYMPIANS AND THE UFO   1974"
#' [6] "     Robert B. Beard              FLYING SAUCERS, UFO'S AND EXTRA"
#' @export
mgrepl <- function(patterns, x, ...) {
    Reduce("|",lapply(patterns, grepl, x, ...))
}

fileIgnoreCase <- function(x) {
    dir <- dirname(x)
    base <- basename(x)
    files <- list.files(dir)
    select <- grepl(base, files, ignore.case=T)
    if(any(select)) {
        paste0(dir,"/",files[select][1])
    } else {
        return(NA)
    }
}

tableVector <- function(x, factors=NULL) {
    if(!is.null(factors)) {
        x <- factor(x, levels = factors)
    }
    table(x) -> tab
    structure(as.vector(tab), .Names=names(tab))
}

#' Fast C++ tabulation.
#' @description Takes in a character, integer or factor vector and tabulates the number of times each element appears.  
#' @param x A character, integer or factor vector.  NAs are allowed.  
#' @param sort TRUE if the result names should be sorted alphanumerically.  
#' @return A integer vector of counts of each element.  
#' @examples
#'
#' x <- factor(sample(1e5, 1e8, replace=T))
#' microbenchmark(table(x), tablec(x), times=3)
#' 
#' Unit: milliseconds
#'      expr       min         lq       mean    median         uq        max
#'  table(x) 9777.6457 10479.0949 10717.3382 11180.544 11187.1844 11193.8246
#' tablec(x)  678.0364   685.9467   713.1181   693.857   730.6589   767.4608
#'
#'
#' x <- sample(letters, 1e8, replace=T)
#' microbenchmark(table(x), tablec(x), tablec(x,sort=T), times=3)
#' 
#'Unit: seconds
#'                expr      min       lq     mean   median       uq      max
#'            table(x) 5.778514 5.829125 5.855560 5.879737 5.894083 5.908430
#'           tablec(x) 1.589360 1.589381 1.589516 1.589402 1.589594 1.589786
#' tablec(x, sort = T) 1.589386 1.590520 1.591824 1.591655 1.593044 1.594432
#' @export
tablec <- function(x, sort=F) {
    cl <- class(x)
    if(cl == "character") {
        ret <- tablec_string(x)
    } else if(cl == "integer") {
        ret <- tablec_int(x)
    } else if(cl == "factor") {
        ret <- tablec_factor(x)
    } else {
        stop("x must be a character, integer or factor vector")
    }
    if(sort) {
        return(ret[sort(names(ret), na.last=F)])
    } else {
        return(ret)
    }
}

#' Size of R objects. 
#' @description Prints out the size of all R objects in the environment. 
#' @param env The environment to search (default global environment). 
#' @param units Units to print out for each variable.  
#' @return A data.frame containing the size of each object.  
#' @export
varSizes <- function(env=globalenv(), units="KB") {
    vars <- ls(envir=env)
    if(length(vars) == 0) return(NULL)
    sizes <- sapply(vars, function(v) {
        format(object.size(get(v, envir=env)), units=units)
    })
    data.frame(vars=vars, size=sizes, row.names=1:length(sizes), stringsAsFactors=F)
}

#' Calculate FPKM.
#' @description Calculates FPKM from a count matrix.
#' @param mat An integer matrix representing the counts of a RNA-Seq dataset.  
#' @param gene_lengths A named vector with the length of each gene.  
#' @param uq_norm If TRUE, will perform upper-quartile normalization.  
#' @return A matrix containing FPKM with the same dimensions as mat.  
#' @seealso
#' FPKM-UQ - A normalized read count in which gene expression values, in FPKM, are divided by the 75th percentile value.
#'
#' \url{https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/}
#'
#' \url{http://vinaykmittal.blogspot.com/2013/10/fpkmrpkm-normalization-caveat-and-upper.html}
#' @export
fpkmFromCounts<- function(mat, gene_lengths, uq_norm=T) {
    gene_lengths <- gene_lengths[rownames(mat)]
    gene_lengths <- gene_lengths / 1000
    lib_sizes <- colSums(mat) / 1e6
    mat <- diag(1 / gene_lengths) %*% mat
    mat <- subset(mat, apply(mat, 1, sum) != 0)
    if(uq_norm == F) {
        return(mat)
    } else {
        uqs <- apply(mat, 2, quantile, 0.75)
        print(summary(uqs))
        mat <- mat %*% diag(1 / uqs)
        return(mat)
    } 
}


# plotROC <- function(pred, class, showplot=T, main="ROC Curve", ...) {
    # require(pROC)
    # roc <- roc(class, pred)
    # auc <- auc(roc)
    
    # conf_mat <- data.matrix(table(pred >= 0.5, class))
    # acc <- (conf_mat[1,1] + conf_mat[2,2])/sum(conf_mat)
    
    # specificity <- roc$specificities
    # sensitivity <- roc$sensitivities
    # if(showplot) { 
        # plot(1-specificity, sensitivity, xlim = c(0,1), ylim = c(0,1), sub=paste0("auc=",round(auc,4), "; accuracy=",round(acc,4)), type="l", lwd=2, main=main,...)
        # abline(a=0,b=1, col="gray")
    # }
    # return(list(specificity, sensitivity))
# }

# needs work, kinda broken atm
# safemclapply <- function(arg_list, p_function, mc.cores=32, tasks_per_core_per_round=3, save_folder=".temp/") {
    # require(parallel)
    # require(digest)
    # # require(dplyr)
    
    # substr(
        # digest(list(capture.output(str(arg_list)),
            # p_function,
            # mc.cores,
            # tasks_per_core_per_round,
            # save_folder)),
    # 1,10) -> run_md5
    
    # function_ls <- ls()

    
    # print(paste("md5 hash: ", run_md5))
    # tasks_per_round = mc.cores*tasks_per_core_per_round
    # task_sequence = seq(1,ceiling(length(arg_list)/tasks_per_round))
    # dir.create(save_folder, showWarnings = FALSE)
    
    # # for(lsv in ls(envir=parent.frame(1))) {lockBinding(lsv,parent.frame(1))}
    # # for(lsv in ls(envir=parent.frame(2))) {lockBinding(lsv,parent.frame(2))}
    
    # for(i in task_sequence) {
        # savefile = paste0(save_folder,run_md5,"_",i,".Rds")
        # if(file.exists(savefile) && file.size(savefile) > 0) next
        # batch_indices = seq((i-1)*tasks_per_round+1,min(i*tasks_per_round,length(arg_list)))
        # print(paste("Processing ",batch_indices[1], " to ", tail(batch_indices,1)))
        # mclapply(arg_list[batch_indices], p_function, mc.cores=mc.cores) -> temp
        # saveRDS(temp, file=savefile)
        # gc()
    # }
    
    # # for(lsv in ls(envir=parent.frame(1))) {unlockBinding(lsv,parent.frame(1))}
    # # for(lsv in ls(envir=parent.frame(2))) {unlockBinding(lsv,parent.frame(2))}
    
    # ret_list <- list()
    # for(i in task_sequence) {
        # savefile = paste0(save_folder,run_md5,"_",i,".Rds")
        # readRDS(savefile) -> ret_list[[i]]
    # }
    # do.call(c, ret_list) -> ret
    
    
    # return(ret)
# }

# splitForMC <- function(arg_list, folder) {
    # if(is.null(names(arg_list))) {
        # names(arg_list) <- 1:length(arg_list)
    # }
    # for(i in 1:length(arg_list)) {
        # names(arg_list)[i] -> na
        # saveRDS(arg_list[[i]], file=sprintf("%s/%s.Rds",folder,na))
    # }
# }

#' Parallel split-matrix loop.
#' @description Splits a matrix into subsets based on a factor, and applies a function to each subset.  Typical use case: sum exon count data to gene count data.  
#' @param mat The matrix.
#' @param f A factor of length equal to nrow(mat).  The levels of this factor will split the matrix into subsets.  
#' @param func The function to apply to each subset.  
#' @param .combine The function to combine the results with.  Default is rbind.  Use NA to return a list. 
#' @param mc.cores The number of cores to use.  
#' @return A list or a combined object depending on the .combine parameter.  
#' @examples
#' library(pasilla)
#' library(DEXSeq)
#' library(trqwe)
#' data(pasillaDEXSeqDataSet)
#' exon_counts <- counts(dxd)
#' f <- rowData(dxd)$groupID
#' gene_counts <- mcsplitapply(exon_counts, f, colSums)
#' @export
mcsplitapply <- function(mat, f, func, mc.cores=4, .combine=rbind, ...) {
    if(mc.cores > 1) require(parallel)
    uf <- unique(f)
    ret <- parallel::mclapply(uf, function(fi) {
      func(mat[f == fi,,drop=F])
      }, mc.cores=mc.cores)
    names(ret) <- uf
    if(!is.function(.combine)) {
        return(ret)
    } else {
        ret <- do.call(.combine, ret)
        return(ret)
    }
}

#' Concatenate strings.  
#' @description Concatenates two strings.  
#' @param a First string.  
#' @param b Second string.  
#' @return The concatenated string.
#' @examples
#' 'Hello ' %Q% 'World'
#' 
#' [1] "Hello World"
#' @export
`%Q%` <- function(a, b) paste0(a, b)


#' Append to a vector.  
#' @description Appends the 2nd argument to the 1st.  
#' @param x A vector. 
#' @param value The element to append.  
#' @examples
#' x <- 1:5
#' append(x) <- 6
#' print(x)
#' 
#' [1] 1 2 3 4 5 6
#' @export
`append<-` <- function(x, value) {
    c(x, value)
}

#' Prepend to a vector.  
#' @description Prepends the 2nd argument to the 1st.  
#' @param x A vector. 
#' @param value The element to append.  
#' @examples
#' x <- 1:5
#' prepend(x) <- 6
#' print(x)
#' 
#' [1] 6 1 2 3 4 5
#' @export
`prepend<-` <- function(x, value) {
    c(value, x)
}

#' Highest elements in a vector.
#' @description Finds the top elements in a vector very quickly.  Equivalent ot \code{-sort(-x, partial=1:n)}
#' @param x A numeric vector.  
#' @param n The number of top elements to return.  
#' @param lowest If TRUE, returns the lowest elements instead of the highest. 
#' @param value If TRUE, returns the values of the top elements.  If FALSE, returns the indices.  
#' @return A vector containing the indices or the values of the top elements.  
#' @examples
#' naive_top <- function(x, n) {
#'    -sort(-x, partial=1:n)
#' }
#' x <- runif(1e7)
#' microbenchmark(naive_top(x,100), topn(x,100,value=T), times=10)
#'
#' Unit: milliseconds
#'                    expr       min        lq     mean    median        uq
#'       naive_top(x, 100) 1070.0180 1071.5951 1075.964 1072.3520 1073.9989
#' topn(x, 100, value = T)  433.6682  433.8882  434.771  434.4986  435.6029
#' @seealso
#' \url{http://stackoverflow.com/questions/18450778/}
#' @export
topn <- function(x, n=100, value=F, lowest=F) {
    if(lowest) x <- -x
    nx <- length(x)
    p <- nx-n
    xp <- sort(x, partial=p)[p]
    wh <- x > xp
    if(value) {
        return(x[wh])
    } else {
        return(which(wh))
    }
}

#' Standard error.  
#' @description Calculates the standard error of a sampling distribution.  
#' @param x A vector.  
#' @return The standard error of x.  
#' @examples
#' x <- rnorm(1e3)
#' se(x)
#' [1] 0.03192027
#' @export
se <- function(x) {x<-na.omit(x);sd(x)/sqrt(length(x))}

# rescale <- function(x, min=0, max=1) {
    # (x-min(x))/(max(x)-min(x)) * (max-min) + min
# }

colorInterpolate <- function(x, color1="white", color2="blue") {
    colramp <- colorRamp(c(color1, color2))
    cols <- colramp(rescale(x))
    cols <- cols/255
    ncol(cols) -> nc
    apply(cols, 1, function(r) {
        if(nc == 3) {
            rgb(r[1], r[2], r[3])
        } else {
            rgb(r[1], r[2], r[3], r[4])
        }
    })

}

#devel branch : https://github.com/GuangchuangYu/clusterProfiler/blob/master/R/GMT.R
read.gmt <- function(gmtfile) {
    require(GSEABase)
    gmt <- getGmt(con=gmtfile)
    ont2gene <- stack(geneIds(gmt))
    ont2gene <- ont2gene[, c("ind", "values")]
    colnames(ont2gene) <- c("ont", "gene")
    return(ont2gene)
}


#specific function for pulling down pathway descriptions from GSEA
# pathwayDescriptionPullDown <- function(gmt_file) {
    # require(curl)
    # readLines(gmt_file) -> gmt
    # strsplit(gmt, split ="\t") -> gmt
    # t(sapply(gmt, function(p) {
        # c(p[1], p[2])
    # })) -> pathway_urls
    # sapply(pathway_urls[,2], function(url) {
        # tryCatch({
        # con <- curl(url)
        # lines <- readLines(con)
        # desc <- lines[which(grepl("Brief", lines))+1]
        # desc <- gsub("^.*<td>","",desc)
        # desc <- gsub("</td>$","",desc)
        # desc <- gsub("<a href.+?</a>","",desc)
        # desc <- gsub("  "," ",desc)
        # desc <- gsub("[\\. ]+$","",desc)
        # print(desc)
        # Sys.sleep(.5)
        # return(desc)
        # }, error=function(e) e)
    # }) -> descs
    # data.frame(short=pathway_urls[,1], url=pathway_urls[,2], desc=descs)
# }

#' Cosine distance. 
#' @description Calculates the cosine distance of rows of a matrix.  
#' @param x A matrix.  
#' @return Cosine distance as a dist object.  
#' @examples
#' rbind(c(1,1,1,0,0,0), c(0,0,0,1,1,1), c(1,1,1,0,0,1)) -> x
#' cosineDist(x)
#'           1         2
#' 2 1.0000000
#' 3 0.1339746 0.7113249
#' @seealso
#' \url{http://stackoverflow.com/questions/2535234/find-cosine-similarity-in-r}
#' @export
cosineDist <- function(x){
  as.dist(1 - x%*%t(x)/(sqrt(rowSums(x^2) %*% t(rowSums(x^2))))) 
}


intRangeString <- function(vec) {
    if(length(vec) ==1) return(as.character(vec))
    vec <- sort(vec)
    res <- list()
    res[[1]] <- vec[1]
    for(i in 2:length(vec)) {
        if(vec[i] == vec[i-1] + 1 || vec[i] == vec[i-1]) {
            res[[length(res)]] <- c(res[[length(res)]],vec[i])
        } else {
            res[[length(res)+1]] <- vec[i]
        }
    }
    sapply(res, function(r) {
        if(min(r) == max(r)) {
            return(as.character(max(r)))
        } else {
            return(paste0(min(r), "-", max(r)))
        }
    }) -> res
    paste0(res, collapse=",")
}

# leadingEdgeGenes <- function(geneList, geneSet, exponent=1) {
    # require(DOSE)
    # df <- DOSE:::gseaScores(geneList=geneList, geneSet=geneSet, exponent=exponent,fortify=T)
    # names(geneList)[which(df$position==1)]
    # upreg <- max(df$runningScore) > abs(min(df$runningScore))
    # if(upreg) {
        # set <- 1:which.max(df$runningScore)
    # } else {
        # set <- which.min(df$runningScore):length(df$runningScore)
    # }
    # intersect(names(geneList)[set], geneSet)
# }

scaleVec <- function(x, ...) {
    scale(x, ...)[,1]
}

# getAuc <- function(true_Y, probs) {
    # probsSort = sort(probs, decreasing = TRUE, index.return = TRUE)
    # val = unlist(probsSort$x)
    # idx = unlist(probsSort$ix)  

    # roc_y = true_Y[idx];
    # stack_x = cumsum(roc_y == 0)/sum(roc_y == 0)
    # stack_y = cumsum(roc_y == 1)/sum(roc_y == 1)    

    # auc = sum((stack_x[2:length(roc_y)]-stack_x[1:length(roc_y)-1])*stack_y[2:length(roc_y)])
    # return(list(stack_x=stack_x, stack_y=stack_y, auc=auc))
# }

# getAuc <- function(class, pred) {
    # ord <- order(pred)
    # pred <- pred[ord]
    # class <- class[ord]
    # as.data.frame(table(pred, class)) -> tab
    # # tab$pred <- unique(pred)
    # npos <- sum(subset(tab, class==1)$Freq)
    # nneg <- sum(subset(tab, class==0)$Freq)
    # tpr <- npos - cumsum(subset(tab, class==1)$Freq)
    # tpr <- tpr/npos
    # tpr <- c(0, rev(tpr))
    # fpr <- cumsum(subset(tab, class==0)$Freq)
    # fpr <- fpr/fpr[length(fpr)]
    # fpr <- c(0, fpr)
    # auc <- sum(sapply(1:(nrow(tab)/2), function(i) {
        # (tpr[i] + tpr[i+1])/2 * (fpr[i+1]-fpr[i])
    # }))
    # return(list(auc=auc, tpr=tpr, fpr=fpr))
# }



#style converter for titles
# titleConvert <- function(title, style="camel", sep=" ", original_sep="_") {
    # title <- unlist(strsplit(title, split=original_sep))
    # title <- tolower(title)
    # if(style == "camel") {
        # title <- gsub("^(.)", "\\U\\1", title, perl=T)
    # } else if(style == "upper") {
        # title <- toupper(title)
    # }
    # paste0(title, collapse=sep)
# }

# do.call.args <- function(what, args, quote = FALSE, envir = parent.frame(), ...) {
    # extra_args <- list(...)
    # do.call(what, c(args, extra_args))
# }

#' @export
saveListFF <- function(list, return_ff_list=F, ...) {
    require(ff)
    if(is.null(names(list))) {
        names(list) <- 1:length(list)
    }
    if(any(is.na(names(list)))) {
        stop("list must have names")
    }
    ff_list <- sapply(list, as.ff, simplify=F)
    with(ff_list, {
        ffsave(list=ls(),...)
    })
    if(return_ff_list) return(ff_list)
}

#' @export
loadFFList <- function(file, ...) {
    require(ff)
    ff_env <- new.env()
    ffload(file, envir=ff_env, ...)
    return(as.list(ff_env))
}

#' Cleans leading and trailing whitespace.  
#' @description Removes leading and trailing whitespace in a vector of strings.  
#' @param x A character vector.  
#' @return A vector with leading and trailing whitespace removed.  
#' @examples
#' chop(c("    hello ", "  123  \t"))
#'
#' [1] "hello" "123"
#' @export
chop <- function(x) {
    gsub("\\s+$","",gsub("^\\s+","",x))
}

multipagePlot <- function(plot_list, nrow=2, ncol=2,...) {
  require(grid)
  require(gridExtra)
  while(length(plot_list) > 0) {
    page_plots <- plot_list[seq(min(nrow*ncol, length(plot_list)))]
    plot_list <- plot_list[-seq(nrow*ncol)]
    grid.arrange(grobs=page_plots, ncol=ncol, nrow=nrow,...)
  }
}



#mclapply is broken and has memory issues.  Workaround: save return value to temp folder and read at the end of loop
#to do: still needs work
# null.mclapply <- function(X,FUN,mc.cores,...) {
    # require(parallel)
    # base <- paste0(tempdir(),"/",paste0(sample(c(letters,LETTERS,0:9), replace=T,size=12),collapse=""))
    # message(paste0("Saving results to ",base))
    # mclapply(seq_along(X),function(i) {
        # saveRDS(FUN(X[i]),file=paste0(base,"_",i,".Rds"),compress=F)
        # return(NULL)
    # }, mc.cores=mc.cores,...)
    # message(paste0("Reading results from ",base))
    # mclapply(seq_along(X),function(i) {
        # readRDS(paste0(base,"_",i,".Rds"))
    # },mc.cores=mc.cores,...)  
# }





# testAUC <- function(probs, class) {
    # x <- probs
    # y <- class
    # x1 = x[y==1]; n1 = length(x1); 
    # x2 = x[y==0]; n2 = length(x2);
    # r = rank(c(x1,x2))
    # h <- ???
    # r_aug <- c(0,r[1:n1],n2)
    # w <- sapply(1:n1, function(ii) {
        # (r[ii]+r[ii+2])/2 - 1
    # })
    # auc <- sum(w*h)/(n1*n2)
# }


# which operations actually benefit?
mcReduce <- function(FUN,x,mc.cores=4) {
    require(parallel)
    while(length(x) != 1) {
        x_tail <- tail(x,length(x) %% 2)
        x <- mclapply(1:floor(length(x)/2), function(i) {
            FUN(x[[i*2-1]],x[[i*2]])
        }, mc.cores=mc.cores)
        if(length(x_tail) > 0) x <- c(x,x_tail)
        # print(length(x))
    }
    return(x[[1]])
}

#' Multi-threaded saveRDS
#' @description Uses the pigz utility to improve saving large R objects.  This is compatible with saveRDS and readRDS functions.  Requires pigz (sudo apt-get install pigz on Ubuntu).  
#' @param object An r object to save.  
#' @param file The filename to save to.
#' @param mc.cores How many cores to use in pigz.  The program does not seem to benefit after more than about 4 cores.  
#' @examples
#' x <- sample(1e4, 1e7, replace=T)
#' y <- sample(1e4, 1e7, replace=T)
#' microbenchmark(mcsaveRDS(x, file="temp.Rds"), saveRDS(y, file="temp2.Rds")
#' 
#' Unit: seconds
#'                             expr      min       lq     mean   median       uq
#'  mcsaveRDS(x, file = "temp.Rds") 1.908310 1.908310 1.908310 1.908310 1.908310
#'   saveRDS(y, file = "temp2.Rds") 6.271499 6.271499 6.271499 6.271499 6.271499
#' @seealso
#' \url{http://stackoverflow.com/questions/28927750/}
#' @export
mcsaveRDS <- function(object,file,mc.cores=min(parallel::detectCores(),4)) {
  con <- pipe(paste0("pigz -p",mc.cores," > ",file),"wb")
  saveRDS(object, file = con)
  close(con)
}

#' Multi-threaded readRDS
#' @description Uses the pigz utility to improve loading large R objects.  This is compatible with saveRDS and readRDS functions.  Requires pigz (sudo apt-get install pigz on Ubuntu).  
#' @param file The filename of the rds object. 
#' @param mc.cores How many cores to use in pigz.  The program does not seem to benefit after more than about 4 cores.  
#' @return The R object.  
#' @examples
#' x <- sample(1e4, 1e7, replace=T)
#' saveRDS(x, file="temp.Rds")
#' xmc <- mcreadRDS("temp.Rds")
#' @seealso
#' \url{http://stackoverflow.com/questions/28927750/}
#' @export
mcreadRDS <- function(file,mc.cores=min(parallel::detectCores(),4)) {
  con <- pipe(paste0("pigz -d -c -p",mc.cores," ",file))
  object <- readRDS(file = con)
  close(con)
  return(object)
}

#' Nelson Aelen estimator
#' @description Nelson–Aalen estimator is an estimator of the cumulative hazard function in survival data.  It can be used to compare the overall risks of two groups or used to estimate the number of deaths before a certain time.  
#' @param time The time of each patient.  
#' @param event The death of each patient: 1 for a patient death, 0 for censored.  
#' @return A list containing the cumulative hazard function.  
#' @examples
#' library(survival)
#' data(veteran)
#' with(nelson_aelen_surv(veteran$time, veteran$status), plot(ti, Hi, type="b"))
#' @seealso
#' #https://en.wikipedia.org/wiki/Nelson-Aalen_estimator
#' @export
nelson_aelen_surv <- function(time, event) {
    ord <- order(time)
    time <- time[ord]
    event <- event[ord]
    di <- table(time[event==1])
    ti <- as.numeric(names(di))
    di <- as.vector(di)
    ni <- length(time) - cumsum(di) + di
    list(Hi=cumsum(di/ni),ti=ti)
}

#' Wald-Wolfowitz Runs Tests for Randomness
#' @description This is the k-category asymptotic Z Test with continuity correction.  Imagine rolling a die multiple times to obtain a sequence of rolls.  This statistic tests whether there exists a "run" within the sequence where one particular number comes up more times in a row than expected randomly.  If the test is significant, it can be concluded that the die rolls are not independent.  
#' @param x A vector of items, coerced into a factor.  
#' @return A p-value for the statistical test.  
#' @examples
#' set.seed(1)
#' ww_test(sample(2, 100, replace=T))
#' ww_test(c(sample(6, 90, replace=T),rep(1,10)))
#' @seealso
#' Reference https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Analysis_of_Runs.pdf
#' @export
ww_test <- function(x) {
    x <- factor(x)
    n <- length(x)
    njs <- table(x)
    njs <- structure(.Data=as.vector(njs), .Names=names(njs))
    mu <- (n*(n+1) - sum(njs^2))/n
    sdev_l <- n^2*(n-1)
    sdev_u <- sum(njs^2*(sum(njs^2)+n*(n+1))) - 2*n*sum(njs^3)-n^3
    sdev <- sqrt(sdev_u/sdev_l)
    r <- sum(diff(as.numeric(x)) != 0) + 1
    z <- (r - mu) / sdev
    pvalue <- 2*pnorm(-abs(z))
    return(pvalue)
}

#' Fast AUC
#' @description This function calculates the Area Under the Reciever-Operator Curve from the results of a classifcation model.  
#' @param probs A numeric vector of probabilities or likelihoods for each data point.  
#' @param class A numeric vector where 1 is positive and 0 is negative.  
#' @return The AUC.  
#' @examples
#' library(AppliedPredictiveModeling)
#' data(abalone)
#' library(glmnet)
#' class <- factor(as.numeric(abalone$Type == "M"))
#' data <- data.matrix(abalone[,-1])
#' fit <- cv.glmnet(data, class, family="binomial")
#' probs <- predict(fit, newx=data)[,1]
#' fastAUC(probs, class)
#' @seealso
#' Reference https://stat.ethz.ch/pipermail/r-help/2005-September/079872.html
#' @export
fastAUC <- function(probs, class, method="trqwe") {
    #ROCR is now much faster
    if(method=="ROCR") {
      require(ROCR)
      pred <- ROCR::prediction(probs, class)
      perf <- performance(pred, "auc")
      return(perf@y.values[[1]])
    }
    x <- probs
    y <- class
    x1 = x[y==1]; n1 = length(x1); 
    x2 = x[y==0]; n2 = length(x2); 
    r = rank(c(x1,x2))  
    auc = (sum(r[1:n1]) - n1*(n1+1)/2) / n1 / n2
    return(auc)
}

#' Fast ROC
#' @description Calculates the points in a Reciever-Operator Curve from the results of a classifcation model.  
#' @param probs A numeric vector of probabilities or likelihoods for each data point.  
#' @param class A numeric vector where 1 is positive and 0 is negative.  
#' @return a list containing the ROC curve.  
#' @examples
#' library(AppliedPredictiveModeling)
#' data(abalone)
#' library(glmnet)
#' class <- factor(as.numeric(abalone$Type == "M"))
#' data <- data.matrix(abalone[,-1])
#' fit <- cv.glmnet(data, class, family="binomial")
#' probs <- predict(fit, newx=data)[,1]
#' with(fastROC(probs, class), plot(fpr, tpr, type="l"))
#' @export
fastROC <- function(probs, class) {
    if(is.factor(class)) {class = as.numeric(class) - 1}
    class_sorted <- class[order(probs, decreasing=T)]
    TPR <- cumsum(class_sorted) / sum(class)
    FPR <- cumsum(class_sorted == 0) / sum(class == 0)
    return(list(tpr=TPR, fpr=FPR))
}

#' Precision-Recall curve
#' @description Calculates the points in a Precision-Recall from the results of a classifcation model.  
#' @param probs A numeric vector of probabilities or likelihoods for each data point.  
#' @param class A numeric vector where 1 is positive and 0 is negative.  
#' @return a list containing the ROC curve.  
#' @examples
#' library(AppliedPredictiveModeling)
#' data(abalone)
#' library(glmnet)
#' class <- factor(as.numeric(abalone$Type == "M"))
#' data <- data.matrix(abalone[,-1])
#' fit <- cv.glmnet(data, class, family="binomial")
#' probs <- predict(fit, newx=data)[,1]
#' with(fastPR(probs, class), plot(recall, precision, type="l"))
#' @export
fastPR <- function(probs, class) {
    if(is.factor(class)) {class = as.numeric(class) - 1}
    class_sorted <- class[order(probs, decreasing=T)]
    recall <- cumsum(class_sorted) / sum(class)
    precision <- cumsum(class_sorted) / 1:length(class_sorted)
    return(list(recall=recall, precision=precision))
}

#' Matthew's correlation coefficient
#' @description Calculates Matthew's correlation coefficient.  Requires the Rmpfr package.  
#' @param probs A numeric vector where 1 is predicted positive and 0 is predicted negative.  
#' @param class A numeric vector where 1 is positive and 0 is negative.  
#' @return The MCC score.  
#' @seealso
#' \url{https://en.wikipedia.org/wiki/Matthews_correlation_coefficient}
#' @export
mccscore <- function(probs, class) {
    require(Rmpfr)
    m <- function(x) mpfr(x, precBits=20)
    tp <- m(sum((probs == 1) & (class == 1)))
    tn <- m(sum((probs == 0) & (class == 0)))
    fp <- m(sum((probs == 1) & (class == 0)))
    fn <- m(sum((probs == 0) & (class == 1)))
    as.numeric(((tp * tn) - (fp * fn)) / sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
}

#' Posterior probability adjustment.
#' @description Adjusts the posterior probability of a classifier based on unbalanced datasets.  In classification model where the negative data is randomly under-sampled and all the positive data is used, the adjustment factor (beta) is p(s=1|-) = p(+)/p(-).  I.e., the probability that a negative datapoint is selected in the classifier.  beta ~ N+/N-.  
#' @param probs The original posterior probability.  
#' @param beta The adjustment factor.  
#' @param Nplus The number of positive examples in the real data.  Only used if beta is NULL.  
#' @param Nminus The number of negative examples in the real data.  Only used if beta is NULL.  
#' @return The adjusted posterior probability.  
#' @seealso
#' Dal Pozzolo, Andrea, et al. "Calibrating probability with undersampling for unbalanced classification." Computational Intelligence, 2015 IEEE Symposium Series on. IEEE, 2015.
#' @export
posteriorBalance <- function(probs, beta=NULL, Nplus, Nminus) {
  if(is.null(beta)) {
    beta <- Nplus/Nminus
  }
  beta*probs / (beta*probs-probs+1)
}

#' F1 score.  
#' @description Calculates F1 score from the results of a classification model.  
#' @param probs A numeric vector where 1 is predicted positive and 0 is predicted negative.  
#' @param class A numeric vector where 1 is positive and 0 is negative.  
#' @return The F1 score.  
#' @seealso
#' \url{https://en.wikipedia.org/wiki/F1_score}
#' @export
f1score <- function(probs, class) {
    tp <- sum((probs == 1) & (class == 1))
    fp <- sum((probs == 1) & (class == 0))
    fn <- sum((probs == 0) & (class == 1))
    precision <- tp / (tp + fp)
    recall <- tp / (tp + fn)
    2 * precision * recall / (precision + recall)
}

#' Concordance Index.  
#' @description Calculates the concordance index from the results of a censored survival model.  Very fast compared to other packages.  
#' @param probs The prognostic score of each patient.  
#' @param time The time of each patient.  
#' @param event The death of each patient: 1 for a patient death, 0 for censored. 
#' @return The concordance index.  
#' @export
cindex <- function(probs, time, event)
{
    wh <- which(event==1)
    time_mat <- outer(time, time[wh], ">")
    pred_mat <- outer(probs, probs[wh], "<")
    pred_mat_eq <- outer(probs, probs[wh], "==")
    total <- sum(time_mat)
    concord <- sum(time_mat & pred_mat) + sum(time_mat & pred_mat_eq)*0.5
    concord/total
}

#' Sigmoid Function. 
#' @description Calculates the results of the sigmoid function.  
#' @param probs x Input to the function.  
#' @return Sigmoid(x)
#' @export
sigmoid <- function(x) {
    return(1/(1+exp(-x)))
}

#' Logit Transformation. 
#' @description Calculates the results of the Logit Transformation.  
#' @param x Input to the function (e.g., probabilities from 0 to 1).  
#' @return Logit(x)
#' @export
logit <- function(x) {
    return(log(x/(1-x)))
}

#' Matrix Factor Design.  
#' @description From a factor, returns a design matrix with a column for each level.  
#' @param x A factor.  
#' @param names The name of each instance in the resulting matrix (i.e., the rownames).  
#' @return The design matrix.  
#' @examples
#' matrixFactor(factor(letters))
#' @export
matrixFactor <- function(x, names=NULL) {
    x <- factor(x)
    xlevels <- levels(x)
    x_mat <- matrix(nrow = length(x), ncol = length(xlevels), 0)
    colnames(x_mat) <- xlevels
    if(!is.null(names)) {
        rownames(x_mat) <- names
    }
    col <- match(x, xlevels)
    for(i in 1:length(x)) {
        x_mat[i,col[i]] <- 1
    }
    return(x_mat)
}

#' Make percentage from number
#' @description Make percentage from number
#' @param x A number.
#' @param digits Number of decimal places, default 2.  
#' @return A percentage.
#' @examples
#' make_percent(0.424, 2)
#' @export
make_percent <- function(x, digits=2) {
  paste0(round(100*x,digits),"%")
}


drew_table_plot <- function(my_data_frame, title=NULL) {
    library(gridExtra)
    library(gtable)
    library(grid)
    library(ggplot2)
  tt3 = ttheme_minimal(  core=list(bg_params = list(fill = c("azure3","azure2", "azure1", "white"), col=NA),
                                   fg_params=list(fontface=3, cex = .55)), colhead=list(fg_params=list(col="black", fontface=4L,
                                                                                                       cex = .55)), rowhead=list(fg_params=list(col="white", fontface=3L, cex = .55)))
  g1 = tableGrob(my_data_frame, theme = tt3)
  title = textGrob(title, gp=gpar(fontsize=14))
  padding = unit(5,"mm")
  table = gtable_add_rows(g1, heights = grobHeight(title) + padding, pos = 0)
  table = gtable_add_grob(table, title, 1, 1, 1, ncol(table))
  grid.arrange(table, nrow=1, newpage = FALSE)
}

gg_center_title <- function() {
    theme(plot.title = element_text(hjust = 0.5))
}

gg_rotate_xlabels <- function(angle=90, hjust=1, vjust=0.5, ...) {
    theme(axis.text.x = element_text(angle = angle, hjust = hjust, vjust=vjust))
}

gg_legend_bottom <- function() {
  theme(legend.position="bottom")
}

#https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
gg_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)
}


#' Fast binary KNN classifier
#' @description Fast binary KNN classifier
#' @param distmat A NxN pre-computed distance matrix
#' @param train_idx Train indicies
#' @param test_idx Test indicies
#' @param classes vector length N, 1 or 0
#' @param K Number of nearest neighbors parameter
#' @param mc.cores Number of threads to use
#' @return A prediction vector for the test set based on the class labels of the train set.
#' @export
trqwe_KNN <- function(distmat, train_idx, test_idx, classes, K, mc.cores=1) {
  if(mc.cores==1) {
    sapply(test_idx, function(ti) {
      sum(classes[topn(distmat[ti,train_idx], n=K, lowest=T)])/K
    })
  } else {
    unlist(mclapply(test_idx, function(ti) {
      sum(classes[topn(distmat[ti,train_idx], n=K, lowest=T)])/K
    }, mc.cores=mc.cores))
  }
}

#' Set rownames of data.frame or matrix
#' @description Set rownames of data.frame or matrix and return it, for use with pipes
#' @param df data.frame or matrix
#' @param rownames rownames to add
#' @return df with rownames added
#' @export
add_rownames <- function(df, rownames) {
  rownames(df) <- rownames
  return(df)
}
#' Set colnames of data.frame or matrix
#' @description Set colnames of data.frame or matrix and return it, for use with pipes
#' @param df data.frame or matrix
#' @param colnames colnames to add
#' @return df with colnames added
#' @export
add_colnames <- function(df, colnames) {
  colnames(df) <- colnames
  return(df)
}
