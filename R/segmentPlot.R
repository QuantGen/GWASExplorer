extractSegment <- function(data, chromosome, base_pair_position, padding = 1e6) {
    begin <- base_pair_position - padding
    end <- base_pair_position + padding
    data[
        data$chromosome == chromosome &
        data$base_pair_position >= begin &
        data$base_pair_position <= end,
    ]
}

# name: name of the grob
# title: title of the plot
# variants: data.frame with three columns: base_pair_position (int), SMR_y (double),
# and LBR_d (double)
segmentPlot <- function(name, title, variants) {
    variants$base_pair_position <- variants$base_pair_position / 1e6
    xScale <- c(
        min(variants$base_pair_position),
        max(variants$base_pair_position)
    )
    if (xScale[1] < 0) {
        xScale[1] <- 0
    }
    yScaleSMR <- c(
        0,
        ceiling(max(variants$SMR_y) + ((max(variants$SMR_y) - min(variants$SMR_y)) / 10))
    )
    yScaleLBR <- c(1, 0)
    xDividers <- round(
        seq(
            from = xScale[1],
            to = xScale[2],
            length.out = 5
        ),
        digits = 2
    )
    yDividersSMR <- round(
        seq(
            from = yScaleSMR[1],
            to = yScaleSMR[2],
            length.out = 5
        ),
        digits = 1
    )
    yDividersLBR <- round(
        seq(
            from = yScaleLBR[1],
            to = yScaleLBR[2],
            by = -0.05
        ),
        digits = 1
    )
    grob <- grid::gTree(
        children = grid::gList(
            grid::textGrob(
                label = title,
                vp = grid::vpPath(paste0(name, ".layout"), paste0(name, "_title.layout")),
            ),
            grid::gTree(
                data = variants,
                vp = grid::vpPath(paste0(name, ".layout"), paste0(name, "_smr.layout")),
                childrenvp = grid::viewport(
                    x = grid::unit(1, "npc") - grid::unit(2, "lines"), # leave padding and space for plot label
                    y = grid::unit(1, "npc") - grid::unit(1, "lines"), # leave padding
                    just = c("right", "top"), # makes it easier to position xaxis and yaxis
                    width = grid::unit(1, "npc") - grid::unit(6, "lines"), # leave space for padding (2), yaxis (1), ylab (1), and plot label (1)
                    height = grid::unit(1, "npc") - grid::unit(3, "lines"), # leave space for padding (2) and xaxis (1)
                    xscale = xScale,
                    yscale = yScaleSMR,
                    name = paste0(name, "_smr.vp")
                ),
                children = grid::gList(
                    grid::polylineGrob(
                        x = grid::unit(rep(c(0, 1), length(yDividersSMR)), "npc"),
                        y = grid::unit(rep(yDividersSMR, each = 2), "native"),
                        id = rep(1:length(yDividersSMR), each = 2),
                        gp = grid::gpar(
                            col = "grey",
                            lty = "dashed",
                            lwd = 0.7
                        ),
                        vp = paste0(name, "_smr.vp")
                    ), # horizontal dividers
                    grid::polylineGrob(
                        x = grid::unit(rep(xDividers, each = 2), "native"),
                        y = grid::unit(rep(c(0, 1), length(xDividers)), "npc"),
                        id = rep(1:length(xDividers), each = 2),
                        gp = grid::gpar(
                            col = "grey",
                            lty = "dashed",
                            lwd = 0.7
                        ),
                        vp = paste0(name, "_smr.vp")
                    ), # vertical dividers
                    grid::linesGrob(
                        x = grid::unit(c(0, 1), "npc"),
                        y = grid::unit(8, "native"),
                        gp = grid::gpar(
                            col = "blue",
                            lty = "dashed",
                            lwd = 1.5
                        ),
                        vp = paste0(name, "_smr.vp")
                    ), # pValue threshold
                    grid::pointsGrob(
                        x = variants$base_pair_position,
                        y = variants$SMR_y,
                        gp = grid::gpar(
                            cex = 0.5,
                            pch = 19
                        ),
                        vp = paste0(name, "_smr.vp")
                    ),
                    grid::xaxisGrob(
                        at = xDividers,
                        vp = paste0(name, "_smr.vp")
                    ),
                    grid::yaxisGrob(
                        at = yDividersSMR,
                        vp = paste0(name, "_smr.vp")
                    ),
                    grid::textGrob(
                        label = "-log10(p-value)",
                        x = grid::unit(-3, "lines"),
                        y = grid::unit(0.5, "npc"),
                        rot = 90,
                        vp = paste0(name, "_smr.vp")
                    ),
                    grid::textGrob(
                        label = "Single Marker Regression",
                        x = grid::unit(1, "npc") + grid::unit(1, "lines"),
                        y = grid::unit(0.5, "npc"),
                        rot = -90,
                        vp = paste0(name, "_smr.vp")
                    )
                )
            ),
            grid::gTree(
                data = variants,
                vp = grid::vpPath(paste0(name, ".layout"), paste0(name, "_bayes.layout")),
                childrenvp = grid::viewport(
                    x = grid::unit(1, "npc") - grid::unit(2, "lines"), # leave padding and space for plot label
                    y = grid::unit(1, "npc") - grid::unit(1, "lines"), # leave padding
                    just = c("right", "top"), # makes it easier to position xaxis and yaxis
                    width = grid::unit(1, "npc") - grid::unit(6, "lines"), # leave space for padding (2), yaxis (1), ylab (1), and plot label (1)
                    height = grid::unit(1, "npc") - grid::unit(3, "lines"), # leave space for padding (2) and xaxis (1)
                    xscale = xScale,
                    yscale = yScaleLBR,
                    name = paste0(name, "_bayes.vp")
                ),
                children = grid::gList(
                    grid::polylineGrob(
                        x = grid::unit(rep(c(0, 1), length(yDividersLBR)), "npc"),
                        y = grid::unit(rep(yDividersLBR, each = 2), "native"),
                        id = rep(1:length(yDividersLBR), each = 2),
                        gp = grid::gpar(
                            col = "grey",
                            lty = "dashed",
                            lwd = 0.7
                        ),
                        vp = paste0(name, "_bayes.vp")
                    ), # horizontal dividers
                    grid::polylineGrob(
                        x = grid::unit(rep(xDividers, each = 2), "native"),
                        y = grid::unit(rep(c(0, 1), length(xDividers)), "npc"),
                        id = rep(1:length(xDividers), each = 2),
                        gp = grid::gpar(
                            col = "grey",
                            lty = "dashed",
                            lwd = 0.7
                        ),
                        vp = paste0(name, "_bayes.vp")
                    ), # vertical dividers
                    grid::pointsGrob(
                        x = variants$base_pair_position,
                        y = variants$LBR_d,
                        gp = grid::gpar(
                            cex = 0.5,
                            pch = 19
                        ),
                        vp = paste0(name, "_bayes.vp")
                    ),
                    grid::xaxisGrob(
                        at = xDividers,
                        label = FALSE, # use same labels as top plot
                        main = FALSE, # draw at the top (what a stupid name)
                        vp = paste0(name, "_bayes.vp")
                    ),
                    grid::yaxisGrob(
                        at = yDividersLBR,
                        vp = paste0(name, "_bayes.vp")
                    ),
                    grid::textGrob(
                        label = expression(paste("P(", beta[j] != 0," | data)")),
                        x = grid::unit(-3, "lines"),
                        y = grid::unit(0.5, "npc"),
                        rot = 90,
                        vp = paste0(name, "_bayes.vp")
                    ),
                    grid::textGrob(
                        label = "Bayesian Regression",
                        x = grid::unit(1, "npc") + grid::unit(1, "lines"),
                        y = grid::unit(0.5, "npc"),
                        rot = -90,
                        vp = paste0(name, "_bayes.vp")
                    )
                )
            )
        ),
        childrenvp = grid::vpTree(
            parent = grid::viewport(
                layout = grid::grid.layout(
                    nrow = 4,
                    ncol = 1,
                    heights = grid::unit.c(
                        grid::unit(1, "lines"),
                        grid::unit(rep_len(1, 2), "null"),
                        grid::unit(1, "lines")
                    )
                ),
                name = paste0(name, ".layout")
            ), # no clip, no nothing allowed in layout viewport ...
            children = grid::vpList(
                grid::viewport(
                    layout.pos.row = 1,
                    layout.pos.col = 1,
                    name = paste0(name, "_title.layout"),
                    clip = TRUE
                ),
                grid::viewport(
                    layout.pos.row = 2,
                    layout.pos.col = 1,
                    name = paste0(name, "_smr.layout"),
                    clip = TRUE
                ),
                grid::viewport(
                    layout.pos.row = 3,
                    layout.pos.col = 1,
                    name = paste0(name, "_bayes.layout"),
                    clip = TRUE
                ),
                grid::viewport(
                    layout.pos.row = 4,
                    layout.pos.col = 1,
                    name = paste0(name, "_xlab.layout"),
                    clip = TRUE
                )
            )
        )
    )
    return(grob)
}
