getformula <- function(df, str = "(%.3f)*%s") {
    df %>%
        tidy %>%
        filter(grepl("mzid", term)) %>%
        mutate(product = sprintf(str, estimate, term)) %>%
        pull(product) %>%
        paste(collapse = " + ")
}

getriskset <- function(df, model, nclass = 4) {
    df %>%
        mutate(risk := !!rlang::parse_expr(getformula(model)),
               riskclass = factor(dplyr::ntile(risk, nclass)),
               riskpersd = scale(risk))
                                           
}
