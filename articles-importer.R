importdata <- function(logtransform = FALSE,
                       normalize = FALSE,
                       replacenawithmin = FALSE) {
    importfinriskidata() %>%
        { if (replacenawithmin == TRUE) 
              mutate_at(., vars(starts_with('mzid_')), list(replacewithmin))
          else
              .
        } %>%
        { if (logtransform == TRUE)
            mutate_at(., vars(starts_with('mzid_')), list(log))
          else
              .
        } %>%
        { if (normalize == TRUE)
              mutate_at(., vars(starts_with('mzid_')), list(scale))
          else
              .
        }
}

importfinriskidata <- function(fmeta = "eicdata/metabolites/2015_60_Salomaa_Jain_dataFR02_FU17_2020-02-19.txt.gz",
                               fabund = "eicdata/metabolites/fr02-eicosanoids-MAD.tsv.gz") {
    abund <- read_tsv(fabund, col_types = cols(.default = col_double(),
                                               sample_run_id = col_character(),
                                               Sample_ID = col_character()))
    
    meta <- read_tsv(fmeta, na = "", col_types = cols(.default = col_double(),
                                                  Sample_ID = col_character(),
                                                  Barcode = col_character(),
                                                  WESID = col_character(),
                                                  WGSID = col_character(),
                                                  BATCH = col_character(),
                                                  FID = col_character(),
                                                  CRP = col_character(),
                                                  CDT = col_character(),
                                                  K_TPKS = col_character(),
                                                  K_VKS = col_character(),
                                                  K_M1 = col_character(),
                                                  K_M2 = col_character(),
                                                  APOE_BATCH = col_character())) 

    full_join(meta, abund, by = "Sample_ID") %>%
        dplyr::mutate(female = factor(ifelse(MEN == 0, 1, 0)),
                      hsCRP = case_when(CRP == "<0.1" ~ 0.05, !is.na(CRP) ~ as.numeric(CRP)),
                      pregnant = case_when(GRAVID == 2 ~ 1,
                                           GRAVID == 1 ~ 0,
                                           female == 0 ~ 0)) %>%
        mutate_at(vars(c("female", "BP_TREAT", "KOULGR", "CURR_SMOKE", "Q57X", "PREVAL_DIAB",
                         "DIAB_FAMILYHIST", "INCIDENT_DIAB_T2", "plate")),
                  as.factor)
        
}


replacewithmin <- function(...) {
    list <- c(...)
    min.value <- min(list, na.rm = TRUE)
    ifelse(is.na(list), min.value, list)
}

importinfo <- function(fnames = "eicdata/metabolites/2015_60_Salomaa_Jain_dataFR02_FU17_2020-02-19_names.txt") {
    read_tsv(fnames, col_types = cols(.default = col_character())) %>%
        filter(!is.na(LONGNAME)) %>%
        select(var = VARIABLE, desc = LONGNAME)
}


