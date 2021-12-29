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
    mzid_subset <- grepl("mzid_", name)
    name[mzid_subset] <- sprintf(formatter,
                                 mzrt2mz(name[mzid_subset]),
                                 mzrt2rt(name[mzid_subset]))
    name
}

pub.p <- function(p) {
  p <- as.numeric(p)
  fo <- case_when(p < 0.001 ~ "%.2e",
                  p < 0.01 ~ "%.3f",
                  TRUE ~ "%.2f")
  sprintf(fo, p)
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

getmetabolistes <- function(obj) {
    if (inherits(obj, "data.frame")) {
        return(obj %>% select(starts_with("mzid_")) %>% colnames)
    }
    if (inherits(obj, "coxph")) {
        return(obj %>% tidy %>% pull(term) %>% mygrep(word = "mzid_"))
    }
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

getcategorical <- function(..., values = c("binomial", "categorical"), exclude = c("Sample_ID", "plate")) {
    rbind(...) %>% filter(type %in% values) %>% pull(var) %difference% exclude
}

getcontinuous <- function(..., values = c("continuous")) {
    rbind(...) %>% filter(type %in% values) %>% pull(var)
}

getmissing <- function(...) {
    rbind(...) %>% filter(dropmissing == TRUE) %>% pull(var)
}

getfilters <- function(...) {
    rbind(...) %>% filter(useforfiltering == TRUE) %>% pull(var)
}

getall <- function(...) {
    rbind(...) %>% pull(var)
}

qsub <- function(chr, arg, script, grun = "/apps/general/grun.py", name = "eicosanoiddm") {
    stopifnot(!missing(chr), !missing(script))
    parameters <- list(grun = grun, script = script, name = name, chr = chr, arg = arg, cwd = getwd())
    "%(grun)s -n %(name)s-%(arg)s-%(chr)s -c '%(cwd)s/scripts/%(script)s %(cwd)s %(chr)s %(arg)s' -q havulinna.q --log-dir %(cwd)s/gwaslogs/" %format%
        parameters %>%
        system(ignore.stdout = TRUE, ignore.stderr = TRUE)
}

qsubwait <- function(pattern = "eicosanoiddm-*") {
    dummy <- sprintf("%s/scripts/dummy.sh", getwd())
    sprintf("qsub -sync y -hold_jid '%s'  -q havulinna.q  %s", pattern, dummy) %>% 
        system(ignore.stdout = TRUE, ignore.stderr = TRUE)
}

model.metabolites <- function(model) {
    model %>% tidy %>% filter(grepl("mzid", term)) %>% pull(term)
}

snptest.firstrow <- function(x) {
    firstrow <- colnames(x)
    names(firstrow) <- firstrow
    firstrow[grepl("ID|missing", firstrow)] <- 0
    firstrow[grepl("riskpersd|^mzid", firstrow)] <- "P"
    firstrow[grepl("riskclass", firstrow)] <- "B"
    firstrow[grepl("BATCH|female", firstrow)] <- "D"
    firstrow[grepl("AGE|^C[0-9]", firstrow)] <- "C"
    firstrow %>%
        as.list %>%
        as_tibble
}

snptest.generatefiles <- function(df, chr) {
    df.raw <- sprintf("gwasdata/fr0207_METAG_chr%s.samples", chr) %>%
        read_delim(delim = " ", col_types = cols(ID_1 = col_character(),
                                                 ID_2 = col_character(),
                                                 missing = col_double())) %>%
        filter(ID_1 != 0) %>%
        left_join(., df, by = c("ID_1" = "Sample_ID"))
    
    snptest.firstrow(df.raw) %>%
        rbind(df.raw) %>%
        write_delim(path = sprintf("gwascache/riskscore_chr%s.samples", chr),
                    delim = " ")

    df.raw %>%
        filter(is.na(BL_AGE)) %>%
        pull(ID_1) %>%
        as.data.frame %>%
        write_delim(path = sprintf("gwascache/exclusions_chr%s.exclude", chr),
                    delim = " ",
                    col_names = FALSE)
    TRUE
}

read.gwas <- function(chr) {
    parallel::mclapply(chr, function(x) {
        read_delim(sprintf("gwascache-manual/gwas_%s.out", x), " ", skip = 13, progress = FALSE) %>%
            mutate(chr = as.numeric(gsub("chr(.*):.*", "\\1", rsid)),
                   start = as.numeric(gsub(".*:([^_]*)_.*", "\\1", rsid))) %>%
            filter(!is.na(start), !is.na(chr)) %>%
            filter(cohort_1_hwe > 1e-10,
                   all_maf > 0.01,
                   nchar(alleleA) == 1,
                   nchar(alleleB) == 1)}, mc.cores = 4) %>%
        map_df(~.x) 
}

getrsid <- function(snpmart, chr, pos) {
    lapply(seq(length(chr)), function(i) {
        query <- getBM(attributes=c("refsnp_id", "allele", "chrom_start", "chrom_strand"),
              filters = c("chr_name", "start", "end"),
              values = list(chr[i], pos[i], pos[i]),
              mart = snpmart) %>%
            pull(refsnp_id)
        ret <- ifelse(length(query) == 1, query[1], NA)
        message(sprintf("(%i/%i) chr : %s, pos : %s, id : %s", i, length(chr), chr[i], pos[i], ret))
        ret
    }) %>% map_chr(~ifelse(is.na(.x), NA, paste0(.x)))
}

mybernoulli <- function(x) {
    rbinom(length(x), 1, 1/(1+exp(-x)))
}

colsuffix <- function(df, suffix, exclude = "term") {
    stopifnot(!is_missing(suffix))
    df %>%
        mutate_at(vars(contains("p.value")), .funs = format,  digits=3, scientific=TRUE) %>%
        mutate(mzrt = pub.mzrt.naive(term),
               mean_ci = sprintf("%.2f (%.2f to %.2f)", estimate, conf.low, conf.high)) %>%
        rename_at(vars(-one_of(exclude)), ~paste0(., suffix))
}
