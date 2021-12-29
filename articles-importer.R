importdata <- function(importmethod = importfinriskidata,
                       logtransform = FALSE,
                       normalize = TRUE,
                       replacenawithmin = TRUE) {
    importmethod() %>%
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

fakedata <- function(n = 2438) {
    data.frame(Sample_ID = 1:n,
               ALKI2_FR02=rnorm(n, 91.1, 134.0),
               BL_AGE=rnorm(n, 46.0, 12.7),
               BMI=rnorm(n, 26.7, 4.5),
               DIAB_T2_AGEDIFF=rnorm(n, 14.3, 4.5),
               DIASM=rnorm(n, 78.7, 11.3),
               HBA1C=rnorm(n, 36.3, 5.6),
               HDL=rnorm(n, 1.4, 0.4),
               hsCRP=rnorm(n, 2.3, 4.8),
               LANTIO=rnorm(n, 98.9, 8.8),
               LDL_DIRECT=rnorm(n, 3.3, 0.8),
               PAINO=rnorm(n, 76.6, 15.4),
               PITUUS=rnorm(n, 1.6, 0.1),
               SYSTM=rnorm(n, 133.8, 19.4),
               TRIG=rnorm(n, 1.4, 1.0),
               VYOTARO=rnorm(n, 88.9, 13.4),
               WHR=rnorm(n, 0.8, 0.1),
               BP_TREAT=rbinom(n, 1, 804/6548),
               CURR_SMOKE=rbinom(n, 1, 1767/6548),
               DIAB_FAMILYHIST=rbinom(n, 1, 1628/6548),
               EAST=rbinom(n, 1, 4450/6548),
               female=rbinom(n, 1, 3423/6548),
               INCIDENT_DIAB_T2=rbinom(n, 1, 586/6548),
               INCIDENT_CR_ANYCANC = 0,
               INCIDENT_DIAB_T1 = 0,
               pregnant = 0,
               PREVAL_CR_ANYCANC = 0,
               PREVAL_DIAB_T1 = 0,
               KOULGR = sample(c(1,2,3), n, replace=TRUE, prob=c(0.336, 0.321, 0.326)),
               Q57X = sample(c(1,2,3), n, replace=TRUE, prob=c(0.228, 0.531, 0.241)),
               plate = sample(seq(1, 20), n, replace=TRUE, prob=rep(0.05, 20)),
               mzid_359.218200_4.7483 = rnorm(n, 0, 1),
               mzid_319.228100_4.8285 = rnorm(n, 0, 4),
               mzid_309.207200_3.7555 = rnorm(n, 0, 2),
               mzid_297.243700_4.6404 = rnorm(n, 0, 2),
               mzid_335.295700_6.7525 = rnorm(n, 0, 3),
               mzid_349.202400_2.6023 = rnorm(n, 0, 2),
               mzid_377.233700_2.8922 = rnorm(n, 0, 4),
               mzid_375.221300_3.0278 = rnorm(n, 0, 3),
               mzid_293.212400_3.9343 = rnorm(n, 0, 2)) %>%
        mutate(PREVAL_DIAB_T2 = mybernoulli(1 
                                            + 2*mzid_359.218200_4.7483 
                                            + 3*mzid_319.228100_4.8285 
                                            - 1*mzid_309.207200_3.7555 
                                            + 0.2*BL_AGE + 2*female + 0.5*Q57X),
               INCIDENT_DIAB_T2 =  mybernoulli(2*(mzid_335.295700_6.7525 > 0)
                                               + 3*( mzid_377.233700_2.8922 > 0)
                                               - 2* mzid_293.212400_3.9343)) %>%
        mutate_at(vars(c("female", "DIAB_FAMILYHIST", "BP_TREAT", "Q57X", "plate")), as.factor)
}

simulate_association <- function(df, formula) {
    values <- eval(substitute(mutate(df, xbet = formula))) %>%
        pull(xbet) %>%
        ilogit()
    rbinom(nrow(df), 1, values)
}

importdilgomdata <- function(fmeta = "eicdata/dilgom/2015_60_Salomaa_Jain_dataFR07_FU17_2019-11-24.txt",
                             fabund = "eicdata/dilgom/BAL-DILGOM07/2020.10.06 DILOM 2007 BaL Data.csv",
                               fplate = "eicdata/metabolites/2020.10.06_DILOM_2007_BaL_Data.Key.txt") {

    abund <- read_csv(fabund) %>%
        mutate(term = case_when(Identity == "EIC_62" ~ "mzid_311.223100_2.9230",
                                Identity == "Eicosanoid_8-iso-PGA1 [M-H]" ~ "mzid_335.223200_2.6455",
                                MZ == "279.1967" & RT == "3.72" ~ "mzid_279.196600_3.7247")) %>%
        filter(!is.na(term)) %>%
        select(-Method, -Polarity, -MZ, -RT, -local_Lab, -Identity) %>%
        gather(Sample_ID, value, -term) %>%
        spread(term, value) 
    
    meta <- read_tsv(fmeta) %>%
        mutate(Sample_ID = as.character(PLASMA_DG07_ID))

    meta_with_plate <- read_delim(fplate, " ") %>%
        select(Barcode, plate = Plate)

    full_join(meta, meta_with_plate, by = c("PLASMA_DG07_ID" = "Barcode")) %>%
     full_join(., abund, by = "Sample_ID") %>%
        dplyr::mutate(female = factor(ifelse(MEN == 0, 1, 0)),
                      hsCRP = case_when(CRP == "<0.1" ~ 0.05, !is.na(CRP) ~ as.numeric(CRP)),
                      pregnant = case_when(GRAVID == 2 ~ 1,
                                           GRAVID == 1 ~ 0,
                                           female == 0 ~ 0)) %>%
        mutate_at(vars(c("female", "BP_TREAT", "KOULGR", "CURR_SMOKE", "Q57X", "PREVAL_DIAB",
                         "DIAB_FAMILYHIST", "INCIDENT_DIAB_T2", "plate")),
                  as.factor)
    
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
                                                  APOE_BATCH = col_character())) %>%
        mutate(Sample_ID = ifelse(is.na(Sample_ID), FID, Sample_ID))
    
    full_join(meta, abund, by = "Sample_ID") %>%
        dplyr::mutate(female = factor(ifelse(MEN == 0, 1, 0)),
                      hsCRP = case_when(CRP == "<0.1" ~ 0.05, !is.na(CRP) ~ as.numeric(CRP)),
                      pregnant = case_when(GRAVID == 2 ~ 1,
                                           GRAVID == 1 ~ 0,
                                           female == 0 ~ 0),
                      HOMAIR = FR02_GLUK_NOLLA*FR02_INS_0H/22.5,
                      HOMAB = 20*FR02_INS_0H/(FR02_GLUK_NOLLA-3.5)) %>%
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


