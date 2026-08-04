#!/usr/bin/env Rscript

# Plot CSV output from the Agilent GC/MS Python ingestion pipeline.
# Direct external dependency: ggplot2 only.

ensure_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 is not installed; installing from CRAN...")
    libs <- .libPaths()
    writable <- libs[file.access(libs, 2) == 0]
    target_lib <- if (length(writable)) writable[[1]] else Sys.getenv("R_LIBS_USER")
    if (!nzchar(target_lib)) {
      stop("no writable R library found for installing ggplot2", call. = FALSE)
    }
    if (!dir.exists(target_lib)) {
      dir.create(target_lib, recursive = TRUE, showWarnings = FALSE)
    }
    .libPaths(unique(c(target_lib, .libPaths())))
    ok <- tryCatch({
      install.packages("ggplot2", lib = target_lib,
                       repos = "https://cloud.r-project.org")
      TRUE
    }, error = function(e) {
      message("failed to install ggplot2: ", conditionMessage(e))
      FALSE
    })
    if (!ok || !requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 is required. Install it with install.packages('ggplot2').",
           call. = FALSE)
    }
  }
  suppressPackageStartupMessages(library(ggplot2))
}

usage <- function() {
  cat(paste0(
    "Usage: Rscript ingest/plot_r.R [options]\n\n",
    "Options:\n",
    "  --ingested PATH   ingested/ directory (default: ./ingested)\n",
    "  --out PATH        plot directory (default: <ingested>/plots_r)\n",
    "  --only RUN...     restrict to these run names\n",
    "  --dpi N           PNG resolution (default: 300)\n",
    "  --logy            log-like y axis on the multi-run spectrum grid\n",
    "  --skip-eic        skip EIC plots even if ingested/eic exists\n",
    "  --help            show this help\n"
  ))
}

parse_args <- function(argv) {
  args <- list(ingested = "ingested", out = NULL, only = character(),
               dpi = 300L, logy = FALSE, skip_eic = FALSE)
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (a == "--help" || a == "-h") {
      usage()
      quit(status = 0)
    } else if (a == "--ingested") {
      i <- i + 1L
      if (i > length(argv)) stop("--ingested needs a path", call. = FALSE)
      args$ingested <- argv[[i]]
    } else if (a == "--out") {
      i <- i + 1L
      if (i > length(argv)) stop("--out needs a path", call. = FALSE)
      args$out <- argv[[i]]
    } else if (a == "--dpi") {
      i <- i + 1L
      if (i > length(argv)) stop("--dpi needs an integer", call. = FALSE)
      args$dpi <- as.integer(argv[[i]])
      if (is.na(args$dpi) || args$dpi <= 0L) {
        stop("--dpi must be a positive integer", call. = FALSE)
      }
    } else if (a == "--logy") {
      args$logy <- TRUE
    } else if (a == "--skip-eic") {
      args$skip_eic <- TRUE
    } else if (a == "--only") {
      vals <- character()
      while (i + 1L <= length(argv) && !startsWith(argv[[i + 1L]], "--")) {
        i <- i + 1L
        vals <- c(vals, argv[[i]])
      }
      args$only <- vals
    } else {
      stop(sprintf("unknown option: %s", a), call. = FALSE)
    }
    i <- i + 1L
  }
  args
}

path_join <- function(...) file.path(..., fsep = .Platform$file.sep)

mkdir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

strip_d <- function(name) {
  sub("\\.[dD]$", "", name)
}

sort_key <- function(name) {
  if (toupper(name) == "STD") return(sprintf("%d:%08d:%s", 0L, 0L, name))
  if (grepl("^[0-9]+$", name)) {
    return(sprintf("%d:%08d:%s", 1L, as.integer(name), name))
  }
  digits <- gsub("[^0-9]", "", name)
  n <- if (nchar(digits)) as.integer(digits) else 0L
  sprintf("%d:%08d:%s", 2L, n, name)
}

ordered_runs <- function(data_dir, only = character()) {
  runs <- basename(list.dirs(data_dir, full.names = TRUE, recursive = FALSE))
  runs <- runs[file.exists(path_join(data_dir, runs, "tic.csv"))]
  runs <- runs[order(vapply(runs, sort_key, character(1)))]
  if (length(only)) {
    want <- tolower(strip_d(only))
    runs <- runs[tolower(runs) %in% want]
  }
  runs
}

read_csv <- function(path) {
  if (!file.exists(path)) stop(sprintf("missing input: %s", path), call. = FALSE)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

num <- function(x) as.numeric(x)

sci_label <- function(x) {
  ifelse(is.na(x), "", format(x, scientific = TRUE, digits = 3, trim = TRUE))
}

fmt_mz <- function(x) {
  sub("\\.?0+$", "", format(x, scientific = FALSE, trim = TRUE))
}

base_theme <- function(base_size = 9) {
  theme_classic(base_size = base_size) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(color = "black", face = "bold"),
      plot.title = element_text(color = "black", face = "bold"),
      plot.subtitle = element_text(color = "#333333"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
}

save_plot <- function(plot, filename, width, height, dpi) {
  ggplot2::ggsave(filename, plot = plot, width = width, height = height,
                  units = "in", dpi = dpi, bg = "white", limitsize = FALSE)
}

save_stacked_plots <- function(plots, filename, width, height, dpi,
                               heights = NULL) {
  n <- length(plots)
  if (is.null(heights)) heights <- rep(1, n)
  grDevices::png(filename, width = width, height = height, units = "in",
                 res = dpi, bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    layout = grid::grid.layout(nrow = n, ncol = 1,
                               heights = grid::unit(heights, "null"))
  ))
  for (i in seq_len(n)) {
    print(plots[[i]], vp = grid::viewport(layout.pos.row = i, layout.pos.col = 1))
  }
  invisible(TRUE)
}

find_peaks <- function(x, y, min_rel = 0.02, min_sep = 0.05, limit = 12L) {
  x <- num(x)
  y <- num(y)
  if (length(y) < 3L || max(y, na.rm = TRUE) <= 0) return(integer())
  idx <- 2L:(length(y) - 1L)
  cand <- idx[y[idx] >= y[idx - 1L] & y[idx] > y[idx + 1L]]
  cand <- cand[y[cand] >= min_rel * max(y, na.rm = TRUE)]
  cand <- cand[order(y[cand], decreasing = TRUE)]
  picked <- integer()
  for (i in cand) {
    if (!length(picked) || all(abs(x[i] - x[picked]) >= min_sep)) {
      picked <- c(picked, i)
    }
    if (length(picked) >= limit) break
  }
  picked[order(x[picked])]
}

spaced_rows <- function(df, x_col, y_col, max_labels, min_gap) {
  if (!nrow(df)) return(df[FALSE, , drop = FALSE])
  ord <- order(df[[y_col]], decreasing = TRUE, na.last = NA)
  chosen <- integer()
  for (i in ord) {
    if (length(chosen) >= max_labels) break
    if (!length(chosen) || all(abs(df[[x_col]][i] - df[[x_col]][chosen]) >= min_gap)) {
      chosen <- c(chosen, i)
    }
  }
  df[sort(chosen), , drop = FALSE]
}

read_tic <- function(ingested, run) {
  d <- read_csv(path_join(ingested, "data", run, "tic.csv"))
  d$rt_min <- num(d$rt_min)
  d$tic <- num(d$tic)
  d$run <- run
  d
}

read_spectrum <- function(ingested, run, unit = FALSE) {
  fname <- if (unit) "avg_spectrum_unit.csv" else "avg_spectrum.csv"
  d <- read_csv(path_join(ingested, "data", run, fname))
  d$mz <- num(d$mz)
  d$relative_pct <- num(d$relative_pct)
  d$run <- run
  d
}

plot_tic_run <- function(run, tic, out_path, dpi) {
  pos <- tic$tic[tic$tic > 0]
  floor <- if (length(pos)) max(min(pos), stats::quantile(pos, 0.02) / 5) else 1
  pk1 <- find_peaks(tic$rt_min, tic$tic, min_rel = 0.05, min_sep = 0.12, limit = 12L)
  pk2 <- find_peaks(tic$rt_min, tic$tic, min_rel = 0.004, min_sep = 0.3, limit = 18L)
  labels1 <- data.frame(rt_min = tic$rt_min[pk1], tic = tic$tic[pk1],
                        label = sprintf("%.2f", tic$rt_min[pk1]))
  labels2 <- data.frame(rt_min = tic$rt_min[pk2], tic = pmax(tic$tic[pk2], floor),
                        label = sprintf("%.2f", tic$rt_min[pk2]))

  p1 <- ggplot(tic, aes(rt_min, tic)) +
    geom_line(color = "#1f4e79", linewidth = 0.35) +
    geom_text(data = labels1, aes(label = label), color = "#b03030",
              size = 2.1, vjust = -0.45, check_overlap = TRUE) +
    labs(title = sprintf("TIC - %s", run),
         subtitle = sprintf("%d scans, %.2f-%.2f min", nrow(tic),
                            min(tic$rt_min), max(tic$rt_min)),
         x = NULL, y = "total ion current") +
    scale_y_continuous(labels = sci_label, expand = expansion(mult = c(0, 0.08))) +
    coord_cartesian(ylim = c(0, max(tic$tic) * 1.05)) +
    base_theme()

  log_data <- tic
  log_data$tic <- pmax(log_data$tic, floor)
  p2 <- ggplot(log_data, aes(rt_min, tic)) +
    geom_line(color = "#1f4e79", linewidth = 0.35) +
    geom_text(data = labels2, aes(label = label), color = "#b03030",
              size = 2.0, vjust = -0.45, check_overlap = TRUE) +
    labs(title = "log scale - reveals minor peaks",
         x = "retention time (min)", y = "total ion current (log10)") +
    scale_y_log10(labels = sci_label,
                  limits = c(floor, max(tic$tic) * 3),
                  expand = expansion(mult = c(0, 0.05))) +
    base_theme()

  save_stacked_plots(list(p1, p2), out_path, 11, 6.6, dpi)
}

plot_tic_grid <- function(runs, tics, out_path, dpi) {
  pdata <- do.call(rbind, tics)
  pdata$run <- factor(pdata$run, levels = runs)
  ncol <- if (length(runs) > 4L) 3L else max(1L, length(runs))
  nrow <- ceiling(length(runs) / ncol)
  p <- ggplot(pdata, aes(rt_min, tic)) +
    geom_line(color = "#1f4e79", linewidth = 0.25) +
    facet_wrap(~ run, ncol = ncol, axes = "all", axis.labels = "all") +
    labs(title = "TIC - all runs",
         subtitle = "Common intensity scale",
         x = "retention time (min)", y = "total ion current") +
    scale_y_continuous(labels = sci_label, limits = c(0, max(pdata$tic) * 1.05)) +
    base_theme(8) +
    theme(legend.position = "none")
  save_plot(p, out_path, 4.6 * ncol, 2.5 * nrow, dpi)
}

plot_tic_overlay <- function(runs, tics, out_path, dpi) {
  raw <- do.call(rbind, tics)
  stack <- raw
  offsets <- seq_along(runs) - 1L
  names(offsets) <- runs
  max_by_run <- tapply(stack$tic, stack$run, max)
  stack$value <- stack$tic / max_by_run[stack$run] + offsets[stack$run]
  raw$run <- factor(raw$run, levels = runs)
  stack$run <- factor(stack$run, levels = runs)
  label_data <- data.frame(run = factor(runs, levels = runs),
                           rt_min = min(raw$rt_min),
                           value = offsets + 0.35)

  p1 <- ggplot(raw, aes(rt_min, tic, color = run)) +
    geom_line(linewidth = 0.25, alpha = 0.85) +
    labs(title = "TIC overlay - absolute intensity",
         x = NULL, y = "total ion current") +
    scale_color_manual(values = grDevices::hcl.colors(length(runs), "Dark 3")) +
    scale_y_continuous(labels = sci_label, limits = c(0, max(raw$tic) * 1.05)) +
    base_theme()

  p2 <- ggplot(stack, aes(rt_min, value, color = run)) +
    geom_line(linewidth = 0.25, alpha = 0.85) +
    geom_text(data = label_data, aes(label = run), hjust = 1.05, size = 2.2,
              show.legend = FALSE) +
    labs(title = "TIC stack - each run scaled to its own base peak",
         x = "retention time (min)",
         y = "normalised TIC (offset per run)") +
    scale_color_manual(values = grDevices::hcl.colors(length(runs), "Dark 3")) +
    scale_y_continuous(breaks = NULL) +
    coord_cartesian(clip = "off") +
    base_theme() +
    theme(plot.margin = margin(5.5, 5.5, 5.5, 36),
          legend.position = "none")

  save_stacked_plots(list(p1, p2), out_path, 12, 8, dpi)
}

top_spectrum_labels <- function(spec, n = 15L, min_gap = 7) {
  ord <- order(spec$relative_pct, decreasing = TRUE, na.last = NA)
  chosen <- integer()
  for (i in ord) {
    if (length(chosen) >= n) break
    if (!length(chosen) || all(abs(spec$mz[i] - spec$mz[chosen]) >= min_gap)) {
      chosen <- c(chosen, i)
    }
  }
  spec[sort(chosen), , drop = FALSE]
}

plot_spectrum_run <- function(run, spec, out_path, dpi, kind = c("fine", "unit")) {
  kind <- match.arg(kind)
  step <- if (kind == "unit") "1 u" else "0.1 u"
  min_gap <- if (kind == "unit") 7 else 2.5
  labels1 <- top_spectrum_labels(spec, n = 15L, min_gap = min_gap)
  labels1$label <- if (kind == "unit") {
    sprintf("%d", round(labels1$mz))
  } else {
    format(signif(labels1$mz, 4), trim = TRUE)
  }

  if (kind == "unit") {
    second <- spec[spec$relative_pct > 0.02, , drop = FALSE]
    labels2 <- top_spectrum_labels(second[second$relative_pct < 20, , drop = FALSE],
                                   n = 14L, min_gap = 10)
    labels2$label <- sprintf("%d", round(labels2$mz))
    p2 <- ggplot(second, aes(mz, relative_pct)) +
      geom_segment(aes(xend = mz, y = 0.02, yend = relative_pct),
                   color = "#1f4e79", linewidth = 0.25) +
      geom_text(data = labels2, aes(label = label), color = "#b03030",
                size = 2.0, vjust = -0.35, check_overlap = TRUE) +
      labs(title = "full range, log scale - bins under 0.02% omitted",
           x = "m/z", y = "relative intensity (%, log10)") +
      scale_y_log10(limits = c(0.02, 400),
                    expand = expansion(mult = c(0, 0.05))) +
      coord_cartesian(xlim = c(min(spec$mz) - 5, max(spec$mz) + 5)) +
      base_theme()
  } else {
    centre <- spec$mz[which.max(spec$relative_pct)]
    lo <- centre - 12
    hi <- centre + 12
    second <- spec[spec$mz >= lo & spec$mz <= hi, , drop = FALSE]
    labels2 <- top_spectrum_labels(second, n = 14L, min_gap = 0.45)
    labels2$label <- format(signif(labels2$mz, 4), trim = TRUE)
    p2 <- ggplot(second, aes(mz, relative_pct)) +
      geom_segment(aes(xend = mz, y = 0, yend = relative_pct),
                   color = "#1f4e79", linewidth = 0.35) +
      geom_text(data = labels2, aes(label = label), color = "#b03030",
                size = 2.0, vjust = -0.35, check_overlap = TRUE) +
      labs(title = sprintf("zoom on the base peak, %.1f-%.1f u", lo, hi),
           x = "m/z", y = "relative intensity (%)") +
      scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.05))) +
      coord_cartesian(xlim = c(lo, hi)) +
      base_theme()
  }

  p1 <- ggplot(spec, aes(mz, relative_pct)) +
    geom_segment(aes(xend = mz, y = 0, yend = relative_pct),
                 color = "#1f4e79", linewidth = 0.25) +
    geom_text(data = labels1, aes(label = label), color = "#b03030",
              size = 2.1, vjust = -0.35, check_overlap = TRUE) +
    labs(title = sprintf("Summed spectrum, all scans, %s bins - %s", step, run),
         subtitle = sprintf("%d bins", nrow(spec)),
         x = NULL, y = "relative intensity (% of base peak)") +
    scale_y_continuous(limits = c(0, 105), expand = expansion(mult = c(0, 0.05))) +
    coord_cartesian(xlim = c(min(spec$mz) - 5, max(spec$mz) + 5)) +
    base_theme()

  save_stacked_plots(list(p1, p2), out_path, 11, 6.6, dpi)
}

plot_spectrum_grid <- function(runs, specs, out_path, dpi, logy = FALSE) {
  pdata <- do.call(rbind, specs)
  pdata$run <- factor(pdata$run, levels = runs)
  if (logy) {
    pdata <- pdata[pdata$relative_pct > 0, , drop = FALSE]
    pdata$plot_value <- pmax(pdata$relative_pct, 1e-3)
    ylab <- "relative intensity (%, log10)"
  } else {
    pdata$plot_value <- pdata$relative_pct
    ylab <- "relative intensity (%)"
  }
  labels <- do.call(rbind, lapply(specs, function(d) {
    labs <- top_spectrum_labels(d, n = 4L, min_gap = 12)
    labs$plot_value <- if (logy) pmax(labs$relative_pct, 1e-3) else labs$relative_pct
    labs$label <- sprintf("%d", round(labs$mz))
    labs
  }))
  labels$run <- factor(labels$run, levels = runs)
  ncol <- if (length(runs) > 4L) 3L else max(1L, length(runs))
  nrow <- ceiling(length(runs) / ncol)
  p <- ggplot(pdata, aes(mz, plot_value)) +
    geom_segment(aes(xend = mz, y = if (logy) 1e-3 else 0, yend = plot_value),
                 color = "#1f4e79", linewidth = 0.2) +
    geom_text(data = labels, aes(label = label), color = "#b03030",
              size = 2, vjust = -0.3, check_overlap = TRUE) +
    facet_wrap(~ run, ncol = ncol, axes = "all", axis.labels = "all") +
    labs(title = "Summed spectrum, 1 u bins - all runs",
         x = "m/z", y = ylab) +
    base_theme(8) +
    theme(legend.position = "none")
  if (logy) {
    p <- p + scale_y_log10(limits = c(1e-3, 200),
                           expand = expansion(mult = c(0, 0.05)))
  } else {
    p <- p + scale_y_continuous(limits = c(0, 110),
                                expand = expansion(mult = c(0, 0.05)))
  }
  save_plot(p, out_path, 4.6 * ncol, 2.5 * nrow, dpi)
}

eic_columns <- function(d) {
  grep("^eic_", names(d), value = TRUE)
}

long_eic <- function(d, run) {
  cols <- eic_columns(d)
  if (!length(cols)) return(data.frame())
  out <- lapply(cols, function(col) {
    data.frame(run = run, rt_min = num(d$rt_min), ion = sub("^eic_", "m/z ", col),
               mz = sub("^eic_", "", col), value = num(d[[col]]))
  })
  do.call(rbind, out)
}

plot_eic_run <- function(run, eic, peaks, out_path, dpi) {
  cols <- eic_columns(eic)
  if (!length(cols)) return(FALSE)
  plots <- list()
  tic_panel <- data.frame(rt_min = num(eic$rt_min), value = num(eic$tic))
  plots[[1L]] <- ggplot(tic_panel, aes(rt_min, value)) +
    geom_line(color = "#444444", linewidth = 0.3) +
    labs(title = sprintf("EIC - %s", run),
         subtitle = "TIC (reference)",
         x = NULL, y = "TIC") +
    scale_y_continuous(labels = sci_label,
                       limits = c(0, max(tic_panel$value) * 1.05),
                       expand = expansion(mult = c(0, 0.05))) +
    base_theme(8)

  peak_rows <- if (!is.null(peaks) && nrow(peaks)) {
    peaks[peaks$run == run, , drop = FALSE]
  } else {
    NULL
  }
  for (col in cols) {
    mz <- sub("^eic_", "", col)
    d <- data.frame(rt_min = num(eic$rt_min), value = num(eic[[col]]))
    labels <- data.frame(rt_min = numeric(), value = numeric(), label = character())
    if (!is.null(peak_rows) && nrow(peak_rows)) {
      pp <- peak_rows[peak_rows$mz == mz, , drop = FALSE]
      if (nrow(pp)) {
        pp$height <- num(pp$height)
        pp$rt_min <- num(pp$rt_min)
        pp <- spaced_rows(pp, "rt_min", "height", max_labels = 6L, min_gap = 0.4)
        labels <- data.frame(rt_min = pp$rt_min, value = pp$height,
                             label = sprintf("%.2f", pp$rt_min))
      }
    }
    share <- if (sum(eic$tic) > 0) 100 * sum(d$value) / sum(eic$tic) else 0
    plots[[length(plots) + 1L]] <- ggplot(d, aes(rt_min, value)) +
      geom_area(fill = "#1f4e79", alpha = 0.18) +
      geom_line(color = "#1f4e79", linewidth = 0.3) +
      geom_text(data = labels, aes(label = label), color = "#b03030",
                size = 2, vjust = -0.35, check_overlap = TRUE) +
      labs(title = sprintf("m/z %s - %.2f%% of the total ion current", mz, share),
           x = if (col == tail(cols, 1)) "retention time (min)" else NULL,
           y = sprintf("m/z %s", mz)) +
      scale_y_continuous(labels = sci_label,
                         limits = c(0, max(d$value) * 1.08),
                         expand = expansion(mult = c(0, 0.05))) +
      base_theme(8)
  }
  height <- max(4, 1.55 * length(plots) + 1.1)
  save_stacked_plots(plots, out_path, 11, height, dpi)
  TRUE
}

plot_eic_ion <- function(mz, eics_long, runs, out_path, dpi) {
  pdata <- eics_long[eics_long$mz == mz, , drop = FALSE]
  if (!nrow(pdata)) return(FALSE)
  pdata$run <- factor(pdata$run, levels = runs)
  ymax <- max(pdata$value) * 1.05
  plots <- lapply(runs, function(run) {
    d <- pdata[pdata$run == run, , drop = FALSE]
    ggplot(d, aes(rt_min, value)) +
      geom_area(fill = "#1f4e79", alpha = 0.18) +
      geom_line(color = "#1f4e79", linewidth = 0.3) +
      labs(title = if (run == runs[[1L]]) sprintf("EIC m/z %s - all runs", mz) else run,
           subtitle = if (run == runs[[1L]]) "Common intensity scale" else NULL,
           x = if (run == tail(runs, 1)) "retention time (min)" else NULL,
           y = run) +
      scale_y_continuous(labels = sci_label, limits = c(0, ymax),
                         expand = expansion(mult = c(0, 0.05))) +
      base_theme(8)
  })
  save_stacked_plots(plots, out_path, 11, max(4, 1.25 * length(plots) + 1.2), dpi)
  TRUE
}

main <- function() {
  ensure_ggplot2()
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  ingested <- normalizePath(args$ingested, mustWork = FALSE)
  if (is.null(args$out)) args$out <- path_join(ingested, "plots_r")
  data_dir <- path_join(ingested, "data")
  if (!dir.exists(data_dir)) {
    stop(sprintf("no ingested data at %s - run python3 ingest/run_all.py first",
                 data_dir), call. = FALSE)
  }

  runs <- ordered_runs(data_dir, args$only)
  if (!length(runs)) stop("no matching runs", call. = FALSE)

  out <- normalizePath(args$out, mustWork = FALSE)
  mkdir(out)
  mkdir(path_join(out, "tic"))
  mkdir(path_join(out, "avg_spectrum"))
  mkdir(path_join(out, "avg_spectrum_unit"))

  tics <- list()
  unit_specs <- list()
  written <- 0L
  for (run in runs) {
    tic <- read_tic(ingested, run)
    tics[[run]] <- tic
    plot_tic_run(run, tic, path_join(out, "tic", paste0(run, ".png")), args$dpi)

    fine <- read_spectrum(ingested, run, unit = FALSE)
    plot_spectrum_run(run, fine,
                      path_join(out, "avg_spectrum", paste0(run, ".png")),
                      args$dpi, kind = "fine")

    unit <- read_spectrum(ingested, run, unit = TRUE)
    unit_specs[[run]] <- unit
    plot_spectrum_run(run, unit,
                      path_join(out, "avg_spectrum_unit", paste0(run, ".png")),
                      args$dpi, kind = "unit")
    written <- written + 3L

    top <- order(unit$relative_pct, decreasing = TRUE)[seq_len(min(5L, nrow(unit)))]
    message(sprintf("  %-5s TIC max %.3g @ %.2f min | top m/z %s",
                    run, max(tic$tic), tic$rt_min[which.max(tic$tic)],
                    paste(sprintf("%d", round(unit$mz[top])), collapse = ", ")))
  }

  plot_tic_grid(runs, tics, path_join(out, "tic_grid.png"), args$dpi)
  plot_tic_overlay(runs, tics, path_join(out, "tic_overlay.png"), args$dpi)
  plot_spectrum_grid(runs, unit_specs,
                     path_join(out, "avg_spectrum_unit_grid.png"),
                     args$dpi, logy = args$logy)
  written <- written + 3L

  eic_dir <- path_join(ingested, "eic")
  if (!args$skip_eic && dir.exists(eic_dir)) {
    mkdir(path_join(out, "eic"))
    peaks_path <- path_join(eic_dir, "peaks.csv")
    peaks <- if (file.exists(peaks_path)) read_csv(peaks_path) else NULL
    eics_long <- list()
    for (run in runs) {
      eic_path <- path_join(eic_dir, paste0(run, ".csv"))
      if (!file.exists(eic_path)) next
      eic <- read_csv(eic_path)
      if (plot_eic_run(run, eic, peaks,
                       path_join(out, "eic", paste0(run, ".png")),
                       args$dpi)) {
        written <- written + 1L
      }
      eics_long[[run]] <- long_eic(eic, run)
    }
    eics_long <- do.call(rbind, eics_long)
    if (!is.null(eics_long) && nrow(eics_long)) {
      for (mz in unique(eics_long$mz)) {
        if (plot_eic_ion(mz, eics_long, runs,
                         path_join(out, "eic", paste0("mz", mz, ".png")),
                         args$dpi)) {
          written <- written + 1L
        }
      }
    }
  }

  message(sprintf("\n%d PNG files written to %s", written, out))
}

tryCatch(main(), error = function(e) {
  message("ERROR: ", conditionMessage(e))
  quit(status = 1)
})
