library(dplyr)
library(fst)
library(compiler)
library(trqwe)
enableJIT(3)

purge_cache <- function() {
  system("sync")
  system("purge")
}

x <- vector(mode="character", length=1e6)
for(i in 1:length(x)) {
  x[i] <- paste0(sample(letters, replace=T, size=sample(100)), collapse="")
}
x <- as.list(x)

lz4write <- function(x,file, compressor="LZ4") {
  purge_cache()
  start <- Sys.time()
  writeBin(compress_fst(serialize(x, connection=NULL), compressor=compressor), con=file, useBytes=T)
  as.numeric(Sys.time() - start)
}

lz4read <- function(file) {
  purge_cache()
  start <- Sys.time()
  xu <- unserialize(decompress_fst(readBin(file, "raw", n=file.info(file)$size)))
  as.numeric(Sys.time() - start)
}

saveR <- function(x, file, compress=T) {
  purge_cache()
  start <- Sys.time()
  saveRDS(x, file="ctest.Rds", compress=compress)
  as.numeric(Sys.time() - start)
}

readR <- function(file) {
  purge_cache()
  start <- Sys.time()
  xu <- readRDS(file="ctest.Rds")
  as.numeric(Sys.time() - start)
}
mcsaveR <- function(x, file) {
  purge_cache()
  start <- Sys.time()
  mcsaveRDS(x, file="ctest.Rds")
  as.numeric(Sys.time() - start)
}

mcreadR <- function(file) {
  purge_cache()
  start <- Sys.time()
  xu <- mcreadRDS(file="ctest.Rds")
  as.numeric(Sys.time() - start)
}

summ <- function(times) {
  cat(mean(times), " +/- ", sd(times), "\n")
}

cat("lz4_write ")
replicate(5, {lz4write(x, "ctest.rs.lz4")}) %>% summ

cat("lz4_read ")
replicate(5, {lz4read("ctest.rs.lz4")}) %>% summ

cat("zstd_write ")
replicate(5, {lz4write(x, "ctest.rs.zstd", compressor = "ZSTD")}) %>% summ

cat("zstd_read ")
replicate(5, {lz4read("ctest.rs.zstd")}) %>% summ

cat("saveRDS ")
replicate(5, {saveR(x, "ctest.Rds")}) %>% summ

cat("readRDS ")
replicate(5, {readR("ctest.Rds")}) %>% summ

cat("saveRDS_uncompressed ")
replicate(5, {saveR(x, "ctest2.Rds", compress=F)}) %>% summ

cat("readRDS_uncompressed ")
replicate(5, {readR("ctest2.Rds")}) %>% summ

cat("mcsaveRDS_uncompressed ")
replicate(5, {mcsaveR(x, "ctest3.Rds")}) %>% summ

cat("mcreadRDS_uncompressed ")
replicate(5, {mcreadR("ctest3.Rds")}) %>% summ

cat("lz4_file_size ", file.info("ctest.rs.lz4")$size /1e6, "\n")
cat("zstd_file_size ", file.info("ctest.rs.zstd")$size /1e6, "\n")
cat("rds_file_size ", file.info("ctest.Rds")$size /1e6, "\n")

