plot.manhattanplot <- function(df,
                               psignf = NULL,
                               nlabels = 10,
                               pointsize = 0.8,
                               xlab = "Eicosanoids Rank Ordered by Mass to Charge Ratio",
                               ylab = "FDR corrected negative Log P Value",
                               color_scheme = c("positive" = "red",
                                                "negative" = "blue",
                                                "insignificant" = "gray")) {
    dat <- df %>%
        dplyr::mutate(neg_log10_pvalue = -1 * log10(qval),
                      coloring = case_when(qval > psignf ~ "insignificant",
                                           estimate < 0 ~ "negative",
                                           estimate >= 0 ~ "positive"),
                      mzrt = pub.mzrt.naive(term),
                      mz = mzrt2mz(term))

  ggplot2::ggplot(dat, ggplot2::aes(x = mz, y = neg_log10_pvalue)) +
    ggplot2::geom_point(ggplot2::aes(color = coloring), size = pointsize) +
    ggrepel::geom_text_repel(data = top_n(dat, n = nlabels, wt = neg_log10_pvalue),
                             ggplot2::aes(label = mzrt),
                             segment.size = 0.05,
                             min.segment.length = unit(1, 'mm'),
                             size = 2) +
      ggplot2::geom_hline(yintercept = -1 * log10(psignf), linetype = 2) +
      ggplot2::scale_x_continuous(name = xlab, breaks = seq(200, 650, 50)) +
      ggplot2::scale_y_continuous(name = ylab, expand = c(0.02, 0)) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::scale_color_manual(values = color_scheme)
}

mycorplot <- function(dset, file, height, width) {
    png(height = height, width = width, pointsize = 14, res = 300, file = file)
    gplots::heatmap.2(dset,
                      notecex = 2,
                      notecol = "black",
                      distfun = dist,
                      hclustfun = hclust,
                      trace = "none",
                      col = colorRampPalette(c("blue", "white", "red"))(n = 300),
                      dendrogram = "column",
                      density.info = "none",
                      labRow = rownames(dset),
                      labCol = NA,
                      asp = 1,
                      lmat = rbind(c(0, 3), c(2, 1), c(0, 4)),
                      lhei = c(2, 10, 1),
                      lwid = c(0.001, 9),
                      margins = c(0, 10),
                      key.par = list(mar=c(2.0, 1.8, 0.5, 15.0), cex=0.5),
                      cexRow=0.25,
                      keysize=0.50,
                      key.title = NA,
                      key.xlab = "")
    dev.off()
}

replicationforestplot <- function(dset, file) {
    png(width = 1000, height = 520, res = 150, file = file)
    forestplot::forestplot(
                    labeltext = cbind(c("Eicosanoid", NA, dset$name),
                                      c("FINRISK\nHR (95% CI)", NA, dset$mean_ci)),
                    mean = cbind(c(NA, NA, dset$estimate)),
                    lower = cbind(c(NA, NA, dset$conf.low)),
                    upper = cbind(c(NA, NA, dset$conf.high)),
                    align = c("l", "l", "l"),
                    graph.pos = 3,
                    title = "",
                    xlog = FALSE,
                    xlab = "HR (95% CI)",
                    txt_gp = fpTxtGp(label = gpar(cex = 1.25),
                                     ticks = gpar(cex = 1.1),
                                     xlab = gpar(cex = 1.2),
                                     title = gpar(cex = 1.2)),
                    xticks = seq(0.0, 2.0, 0.5),
                    clip =exp(c(-1, 1)),
                    col = fpColors(box=c("#0433ff")),
                    zero = 1, 
                    lineheight = unit(16, "mm"),
                    boxsize = 0.15,
                    colgap = unit(4, "mm"),
                    lwd.ci = 1)
    dev.off()
}
