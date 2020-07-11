## Helper functions

c2l <- function(...) {
    l <- as.list(c(...))
    names(l) <- c(...)
    l
}

mygrep <- function(..., word, ignorecase = TRUE, complement = FALSE) {
    c(...)[xor(grepl(word, c(...), ignore.case = ignorecase), (complement == TRUE))]
}

dropna <- function(...) {
    c(...)[!is.na(c(...))]
}

pass <- function(x, fun) {
    fun(x)
    x
}

lastelement <- function(x) {
    x[[length(x)]]
}

## Operators

`%difference%` <- function(a, b) {
    a[!a %in% b]
}

`%intersect%` <- function(a, b) {
    intersect(a, b)
}


`%union%` <- function(a, b) {
    c(a, b)
}

## General functions

mytidy <- function(model) {
    dependant <- all.vars(model$call)[1]
    tidy(model, exponentiate = FALSE) %>%
        tibble::add_column(response = dependant, .before = 1)
}

myformula <- function(response, term, covariates) {
    terms <- term %union% covariates
    sprintf("%s ~ %s", response, paste(terms, collapse = " + "))
}

loop.results <- function(..., filterstr = "mzid_", exponentiate = FALSE, method = "fdr") {
    purrr::map_df(c(...), ~broom::tidy(.x, exponentiate = exponentiate)) %>%
        dplyr::filter(grepl(filterstr, term)) %>%
        dplyr::mutate(conf.low = estimate - qnorm(0.975) * std.error,
                      conf.high = estimate + qnorm(0.975) * std.error,
                      qval = p.adjust(p.value, method = method))
}

## Linear and logistic models

loop.lm <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    models <- parallel::mclapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates)
        ret <- stats::lm(formula = as.formula(fo), data = dset, na.action = na.omit)
        ret$call <- as.formula(fo)
        ret
    }, mc.cores = min(length(loops), 4))
}

loop.binomial <- function(dset,
                          response,
                          loops,
                          covariates = c()) {
    stopifnot(!missing(dset), !missing(response), !missing(loops))
    lapply(c2l(loops), function(loop) {
        fo <- myformula(response, loop, covariates)
        ret <- stats::glm(formula = as.formula(fo),
                          family=binomial(link='logit'),
                          data = dset,
                          na.action = na.omit)
        ret$call <- as.formula(fo)
        ret
    })
}

## COX

mysurvformula <- function(response, term, covariates) {
    terms <- term %union% covariates
    sprintf("Surv(%s_AGEDIFF, INCIDENT_%s) ~ %s", response, response, paste(terms, collapse = " + "))
}

loop.cox <- function(dset,
                    response,
                    loops,
                    covariates = c()) {
    parallel::mclapply(c2l(loops), function(loop) {
        fo <- mysurvformula(response, loop, covariates)
        ret <- survival::coxph(formula = as.formula(fo), ties = "breslow", data = dset)
        ret$call <- str2lang(sprintf('survival::coxph(formula = %s, ties = "breslow", data = dset)', fo))
        ret
    }, mc.cores = min(length(loops), 4)) %>%
        { if (length(.) == 1) .[[1]] else . }
}

loop.coxresiduals <- function(...) {
    lapply(c(...), function(model) {
        term <- all.vars(model$call)[3]
        survminer::ggcoxdiagnostics(model,
                         type = "deviance",
                         linear.predictions = FALSE,
                         point.shape = ".",
                         hline.col = "gray",
                         hline.lty = "solid",
                         hline.size = 0.5,
                         sline.lty = "solid",
                         sline.size = 0.5,
                         ggtheme = theme_bw()) +
            ggtitle(sprintf("%s", bioproperty(term))) +
            theme(plot.title = element_text(size = 11),
                  aspect.ratio = 1,
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  strip.background = element_blank(),
                  strip.text.x = element_blank())
    })
}

check.proportionality <- function(...) {
    lapply(c(...), function(model) cox.zph(model))
}

filter.significant <- function(..., filterstr = "mzid_") {
    lapply(c(...), function(model) {
        broom::tidy(model) %>%
            dplyr::filter(grepl(filterstr, term)) %>%
            dplyr::top_n(n = 1, wt = -p.value) %>%
            dplyr::pull(p.value) %>%
            { if (. < 0.05) model else NA }
    }) %>%
        dropna
}

loop.residuals <- function(...) {
    lapply(c(...), function(model) {
        term <- ifelse(isS4(model),
                       all.vars(model@call)[2],
                       all.vars(model$call)[2])
        ggplot(broom::augment(model), aes(x = .fitted, y = .resid)) +
            geom_hline(yintercept = 0, color = "gray") +
            geom_point(size = 1) +
            geom_smooth(method = 'gam', formula = y ~ s(x, bs = "cs"), size = 1) +
            ggtitle(sprintf("%s", bioproperty(term))) +
            theme_classic() +
            theme(plot.title = element_text(size = 11),
                  aspect.ratio = 1,
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  strip.background = element_blank(),
                  strip.text.x = element_blank())

    })
}

loop.qq <- function(...) {
    lapply(c(...), function(model) {
        term <- ifelse(isS4(model),
                       all.vars(model@call)[2],
                       all.vars(model$call)[2])
        ggplot(data = broom::augment(model), mapping = aes(sample = .std.resid)) +
            geom_abline(intercept = 0, size = 1, slope=1, color = "gray") +
            stat_qq(geom = "point", size = 0.1) +
            coord_fixed() +
            ggtitle(sprintf("%s", bioproperty(term))) +
            theme_classic() +
            theme(plot.title = element_text(size = 11),
                  aspect.ratio = 1,
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  strip.background = element_blank(),
                  strip.text.x = element_blank())
    })
}


results.table <- function(df, exponentiate = c("htn", "htn3", "htn.followup"), percent = FALSE, drop = FALSE) {
    dplyr::mutate(df,
                  term = bioproperty(term),
                  estimate = ifelse(response %in% exponentiate, exp(estimate), estimate),
                  conf.low = ifelse(response %in% exponentiate, exp(conf.low), conf.low),
                  conf.high = ifelse(response %in% exponentiate, exp(conf.high), conf.high)) %>%
        betacip(percent = percent) %>%
        myspread() %>%
        { if (drop) filter_at(., vars(contains("p.value")), any_vars(. < 0.05)) else .}
}

# Others

spearmancorrelation  <- function(dset, vars) {
    dat <- dset %>% dplyr::select(one_of(vars))
    colnames(dat) <- colnames(dat) %>% bioproperty()
    cor(dat, method = 'spearman')
}
