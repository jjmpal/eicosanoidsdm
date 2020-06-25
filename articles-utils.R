mzrt2mz <- function(mzid) {
  mzid %>%
    gsub("mzid_", "", .) %>%
    gsub("_.*", "", .) %>%
    as.numeric()
}

mzrt2rt <- function(mzid) {
  mzid %>%
    gsub("mzid_.*_", "", .) %>%
    as.numeric()
}

bioproperty <- function(name, formatter = "%.4f / %.4f") {
    mzid.subset <- grepl("mzid_", name)
    name[mzid.subset] <- sprintf(formatter,
                                 mzrt2mz(name[mzid.subset]),
                                 mzrt2rt(name[mzid.subset]))
    name
}

gmean <- function(x, na.rm=TRUE){
  if (any(x < 0, na.rm = TRUE)) return(NaN)
  if (any(x == 0, na.rm = TRUE)) return(0)
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

pub.mzrt.naive <- function(mzid, formatter = "%.4f / %.4f") {
    sprintf(formatter,
            as.numeric(gsub("_.*", "", gsub("mzid_", "", mzid))),
            as.numeric(gsub(".*_", "", mzid)))
}

getmetabolistes <- function(dset) {
    dset %>% select(starts_with("mzid_")) %>% colnames
}

mkdir <- function(...) {
    vmkdir <- Vectorize(dir.create, vectorize.args = "path")
    vmkdir(path = c(...), showWarnings = FALSE) %>%
        invisible
}

rmdir <- function(...) {
    vrmdir <- Vectorize(unlink, vectorize.args = "x")
    vrmdir(x = c(...), recursive=TRUE) %>%
        invisible
}

getcategorical <- function(..., values = c("binomial", "categorical")) {
    c(...)[c(...) %in% values] %>% names
}

getcontinuous <- function(..., values = c("continuous")) {
    c(...)[c(...) %in% values] %>% names
}

getall <- function(...) {
    c(getcategorical(...), getcontinuous(...))
}

