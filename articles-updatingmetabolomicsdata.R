library(dplyr)
library(tibble)
library(broom)
library(tidyr)

mzidstandarizer <- function(mz, rt, formatter = "mzid_%0.6f_%0.4f") {
  sprintf(formatter, mz, rt)
}

mycolrenamer <- function(old, rename, leadingzero = FALSE) {
    ret <- gsub(".*Plate_([0-9][0-9]?)_([0-9][0-9]?)_.*driftCorrect.*", "\\1_\\2", old)
    if (leadingzero) {
        ret <- gsub(".*Plate_([0-9][0-9]?)_([0-9][0-9]?)_.*driftCorrect.*", "00\\1_00\\2", old)
    }
    
    if (!missing(rename)) {
        return(rename[match(ret, pull(rename, PLATEWELL)),  "Sample_ID"])
    }
    return(ret)
}

# FILES

files <- list("FHS" = list("source" = "eicdata/metabolites/elokuu/FHS_EIC_MAD_Norm_EIC_only.csv",
                           "target" =  "eicdata/metabolites/fhs-eicosanoids-MAD.rds"),
              "FR02" = list("source" = "eicdata/metabolites/elokuu/FINRISK_EIC_MAD_Norm_EIC_only.csv",
                            "target" =  "eicdata/metabolites/fr02-eicosanoids-MAD.rds"))

lapply(files, function(file) {
    dset <- read.csv(file = file$source, header=TRUE, sep=",", fileEncoding="latin1") %>%
        dplyr::mutate(mzid = mzidstandarizer(MZ, RT))

    df.data.abundances <- dset %>%
        select(-MZ, -RT, -local_Lab, -Identity) %>%
        gather(variable, value, -mzid) %>%
        mutate(variable = mycolrenamer(variable, leadingzero = grepl("FHS", file$source))) %>%
        spread(mzid, value) %>%
        rename(key = variable)

    df.mapping <- dset %>% 
        mutate(name = gsub("^\\s+|\\s+$", "", Identity)) %>%
        select(mzid, name, MZ, RT)

    saveRDS(list(data = df.data.abundances,
             mapping = df.mapping,
             desc = file$source,
             date = format(Sys.time(), format="%Y%m%d-%H%M%S")),
            file = file$target)
    })

