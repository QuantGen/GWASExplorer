# name: name of the grob
# title: title of the plot
# chromosomes: data.frame with two columns: chromosome (factor), length (int)
# segments: data.frame with four columns: chromosome, start, length, pValue
# variants: data.frame with two columns: chromosome, position
chromosomePlot <- function(name, title, chromosomes, segments, variants, pValuePalette = grDevices::colorRamp(c("yellow", "red"))) {
    bpLimit <- shiftceil(chromosomes$length)
    pValueColors <- grDevices::rgb(
        pValuePalette(range01(segments$pValue)),
        maxColorValue = 255
    )
    grob <- grid::gTree(
        childrenvp = grid::viewport(
            x = grid::unit(1, "npc") - grid::unit(1, "lines"), # leave padding
            y = grid::unit(1, "npc") - grid::unit(3, "lines"), # leave padding
            just = c("right", "top"), # makes it easier to position xaxis and
                                      # yaxis
            width = grid::unit(1, "npc") - grid::unit(6, "lines"), # leave space for yaxis
            height = grid::unit(1, "npc") - grid::unit(5, "lines"), # leave padding and
                                                        # space for xaxis
            xscale = range(as.integer(chromosomes$chromosome)) + c(-0.5, 0.5),
            yscale = c(0, bpLimit),
            name = paste0(name, ".vp")
        ),
        children = grid::gList(
            # chromosomes
            grid::rectGrob(
                x = grid::unit(as.integer(chromosomes$chromosome), "native"),
                y = grid::unit(0, "native"),
                just = c("center", "bottom"),
                width = grid::unit(0.5, "native"),
                height = grid::unit(chromosomes$length, "native"),
                gp = grid::gpar(
                    lwd = 0,
                    fill = "grey",
                    col = "grey"
                ),
                vp = paste0(name, ".vp")
            ),
            # p-values
            grid::rectGrob(
                x = grid::unit(as.integer(segments$chromosome), "native"),
                y = grid::unit(segments$start, "native"),
                just = c("center", "bottom"),
                width = grid::unit(0.5, "native"),
                height = grid::unit(segments$length, "native"),
                gp = grid::gpar(
                    lwd = 0,
                    fill = pValueColors,
                    col = pValueColors # lwd = 0 by itself doesn't work in
                                       # pdf() to eliminate borders
                ),
                vp = paste0(name, ".vp")
            ),
            # Bayesian hits
            grid::rectGrob(
                x = grid::unit(as.integer(variants$chromosome), "native"),
                y = grid::unit(variants$base_pair_position, "native"),
                just = c("center", "bottom"),
                width = grid::unit(0.3, "native"),
                height = grid::unit(3e5, "native"), # pdf() produces a minimum height
                                              # just to be visible, for
                                              # everything else we need to be
                                              # specific
                gp = grid::gpar(
                    lwd = 0,
                    fill = "blue",
                    col = "blue" # lwd = 0 by itself doesn't work in pdf() to
                                 # eliminate borders
                ),
                vp = paste0(name, ".vp")
            ),
            grid::xaxisGrob(
                at = 1:nrow(chromosomes),
                label = levels(chromosomes$chromosome),
                vp = paste0(name, ".vp")
            ),
            grid::linesGrob(
                x = grid::unit(c(0, 1), "npc"),
                y = grid::unit(c(0, 0), "npc"),
                vp = paste0(name, ".vp")
            ), # extend range of xaxis to include padding
            grid::yaxisGrob(
                at = seq(0, bpLimit, by = 50e6),
                label = paste0(seq(0, bpLimit / 1e6, by = 50), " Mbp"),
                vp = paste0(name, ".vp")
            ),
            grid::textGrob(
                label = title,
                y = grid::unit(1, "npc") + grid::unit(1, "lines"),
                just = c("center", "top"),
                vp = paste0(name, ".vp")
            )
        )
    )
    return(grob)
}
