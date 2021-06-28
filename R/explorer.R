# Global variables
globals <- new.env(parent = emptyenv())
globals$chromosomeLengths <- loadChromosomeLengths()
globals$segments <- loadSegments()
globals$variants <- loadVariants()
globals$state <- "chromosomeViewer"
globals$currentPlot <- NULL
globals$traits <- c(
    glucose = "Glucose",
    log_urate = "Urate",
    log_creatinine = "Creatinine",
    log_hdl = "HDL",
    ldl = "LDL",
    cholesterol = "Cholesterol",
    log_triglycerides = "Triglycerides"
)
globals$trait <- "cholesterol"
globals$chromosome <- 18
globals$base_pair_position <- 50e6
globals$bfdrCutoff <- 0.05


createChromosomePlot <- function() {
    segments <- extractSegments(globals$segments, globals$trait)
    variants <- extractVariants(globals$variants, globals$trait, globals$bfdrCutoff)
    chromosomePlot(
        name = "chromosomeViewer",
        title = globals$traits[globals$trait],
        chromosomes = globals$chromosomeLengths,
        segments = segments,
        variants = variants
    )
}


createSegmentPlot <- function() {
    variants <- extractVariants(globals$variants, globals$trait)
    segment <- extractSegment(
        data = variants,
        chromosome = globals$chromosome,
        base_pair_position = globals$base_pair_position
    )
    segmentPlot(
        name = "segmentViewer",
        title = paste0(globals$traits[globals$trait], " (Chromosome ", globals$chromosome, ")"),
        variants = segment
    )
}


drawPlot <- function(grob) {
    grid::grid.newpage()
    grid::grid.draw(grob)
}


explore <- function() {
    if (!interactive()) {
        stop("R needs to be run interactively")
    }
    globals$currentPlot <- createChromosomePlot()
    drawPlot(globals$currentPlot)
    while (TRUE) {
        system("clear")
        message("=== GWAS Explorer ===")
        if (globals$state == "chromosomeViewer") {
            response <- readline(
                prompt = paste(
                    "",
                    "Options:",
                    "1) Refresh plot",
                    paste0("2) Change trait (current: ", globals$trait, ")"),
                    paste0("3) Change BFDR cutoff (current: ", globals$bfdrCutoff, ")"),
                    "4) Select segment",
                    "5) Save plot",
                    "6) Quit",
                    "",
                    "Select option: ",
                    sep = "\n"
                )
            )
            if (response == "1") {
                globals$currentPlot <- createChromosomePlot()
                drawPlot(globals$currentPlot)
            } else if (response == "2") {
                response <- readline(
                    prompt = paste(
                        "",
                        paste0("Current trait: ", globals$trait),
                        "",
                        "Available traits:",
                        paste(names(globals$traits), collapse = "\n"),
                        "",
                        "Select trait: ",
                        sep = "\n"
                    )
                )
                if (response %in% names(globals$traits)) {
                    globals$trait <- response
                    globals$currentPlot <- createChromosomePlot()
                    drawPlot(globals$currentPlot)
                } else {
                    message("\nError: Invalid trait ...\n")
                    Sys.sleep(1)
                }
            } else if (response == "3") {
                response <- readline(
                    prompt = paste(
                        "",
                        paste0("Current BFDR cutoff: ", globals$bfdrCutoff),
                        "",
                        "Available BFDR cutoffs:",
                        "0.05",
                        "0.1",
                        "",
                        "Select BFDR cutoff: ",
                        sep = "\n"
                    )
                )
                if (response %in% c("0.05", "0.1")) {
                    globals$bfdrCutoff <- as.double(response)
                    globals$currentPlot <- createChromosomePlot()
                    drawPlot(globals$currentPlot)
                } else {
                    message("\nError: Invalid BFDR cutoff ...\n")
                    Sys.sleep(1)
                }
            } else if (response == "4") {
                message("\nSelect a segment in the plot ...\n")
                grid::downViewport("chromosomeViewer.vp")
                locs <- grid::grid.locator()
                grid::upViewport()
                chromosome <- as.integer(round(locs$x))
                base_pair_position <- as.integer(round(locs$y))
                message("Selected coord: ", base_pair_position)
                if (!chromosome %in% 1:24) {
                    stop("Invalid chromosome selected: ", chromosome)
                } else if (base_pair_position < 0 || base_pair_position > max(globals$chromosomeLengths$length)) {
                    stop("Invalid position selected: ", base_pair_position)
                } else {
                    globals$chromosome <- chromosome
                    globals$base_pair_position <- base_pair_position
                    globals$currentPlot <- createSegmentPlot()
                    globals$state <- "segmentViewer"
                    drawPlot(globals$currentPlot)
                }
            } else if (response == "5") {
                response <- readline(
                    prompt = paste(
                        "",
                        "Choose path to save plot: ",
                        sep = "\n"
                    )
                )
                grDevices::pdf(file = response)
                drawPlot(globals$currentPlot)
                grDevices::dev.off()
            } else if (response == "6") {
                message("Bye bye!")
                quit(save = "no")
            } else {
                message("\nError: Invalid selection ...\n")
                Sys.sleep(1)
            }
        } else if (globals$state == "segmentViewer") {
            response <- readline(
                prompt = paste(
                    "",
                    "Options:",
                    "1) Refresh plot",
                    "2) Return to chromosome plot",
                    "3) Save plot",
                    "4) Quit",
                    "",
                    "Select option: ",
                    sep = "\n"
                )
            )
            if (response == "1") {
                globals$currentPlot <- createSegmentPlot()
                drawPlot(globals$currentPlot)
            } else if (response == "2") {
                globals$currentPlot <- createChromosomePlot()
                globals$state <- "chromosomeViewer"
                drawPlot(globals$currentPlot)
            } else if (response == "3") {
                response <- readline(
                    prompt = paste(
                        "",
                        "Choose path to save plot: ",
                        sep = "\n"
                    )
                )
                grDevices::pdf(file = response)
                drawPlot(globals$currentPlot)
                grDevices::dev.off()
            } else if (response == "4") {
                message("Bye bye!")
                quit(save = "no")
            } else {
                message("\nError: Invalid selection ...\n")
                Sys.sleep(1)
            }
        }
    }
}
