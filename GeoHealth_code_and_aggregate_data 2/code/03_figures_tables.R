# GeoHealth manuscript and Supporting Information, Brazil, 2000–2019.
# Usage: Rscript code/03_figures_tables.R [output_directory]
# Independently renders Figures 1–3, S1–S5 (PDF/PNG) and Tables 1–2, S1–S4 (CSV).
# Requires data.table, ggplot2, gtable, patchwork, ragg, RColorBrewer, scales, terra
# and Ghostscript (gs on PATH). Input CSV files are in ../data relative to this script.
suppressPackageStartupMessages(library(data.table))
figure_file <- function(figure, format) file.path(out_dir, paste0(figure, ".", format))
save_plot <- ggplot2::ggsave
pdf_page_size <- function(pdf) {
  gs <- Sys.which("gs")
  input <- normalizePath(pdf, mustWork = TRUE)
  escaped <- gsub("(", "\\(", gsub(")", "\\)", input, fixed = TRUE), fixed = TRUE)
  script <- tempfile(fileext = ".ps")
  writeLines(c("/document <<>> .PDFInit def", sprintf("(%s) document .PDFFile", escaped), "document 0 .PDFPageInfo /MediaBox get {=only ( ) print} forall", "document .PDFClose"), script)
  result <- system2(gs, c("-q", "-dBATCH", "-dNOPAUSE", "-dNODISPLAY", paste0("--permit-file-read=", shQuote(input)), shQuote(script)), stdout = TRUE, stderr = TRUE)
  box <- scan(text = paste(result, collapse = " "), quiet = TRUE)
  stopifnot(length(box) == 4L)
  box[3:4] - box[1:2]
}

render_pdf_png <- function(pdf, png, dpi) {
  gs <- Sys.which("gs")
  pixels <- ceiling(pdf_page_size(pdf) * dpi / 72)
  status <- system2(gs, c("-q", "-dBATCH", "-dNOPAUSE", "-dFIXEDMEDIA", "-sDEVICE=png16m", "-dTextAlphaBits=4", "-dGraphicsAlphaBits=4", paste0("-r", dpi), sprintf("-g%dx%d", pixels[1], pixels[2]), paste0("-sOutputFile=", shQuote(png)), shQuote(pdf)))
  stopifnot(status == 0L, file.exists(png))
  invisible(png)
}

# Simplified IBGE display polygons in SIRGAS 2000. Rings and features are ordered
# within each layer. Municipal codes index the values for Figures 1 and S3.
read_boundaries <- function(layer_name) {
  d <- map_vertices[layer==layer_name]
  id <- unique(d[,.(geom,key)])
  stopifnot(identical(id$geom,seq_len(nrow(id))))
  geom <- as.matrix(d[,.(geom,part,x,y,hole)])
  attributes <- data.frame(value=id$key)
  names(attributes) <- switch(layer_name, municipalities="Code", regions="region", "layer")
  terra::vect(geom,type="polygons",atts=attributes,crs="EPSG:4674")
}

draw_spatial_panels <- function(data, filename) {
  full_width_in <- 190 / 25.4
  font_family <- "sans"
  suppressPackageStartupMessages(library(terra))
  map_data <- copy(data)
  map_data[, Code := as.character(Code)]

  shape <- read_boundaries("municipalities")
  water_shape <- read_boundaries("water")
  region_shape <- read_boundaries("regions")
  country_shape <- read_boundaries("country")
  shape$MeanTemperatureC <- map_data$MeanTemperatureC[
    match(shape$Code, map_data$Code)
  ]
  shape$K550CrudeRatePer100kPY <- map_data$K550CrudeRatePer100kPY[
    match(shape$Code, map_data$Code)
  ]

  temperature_breaks <- seq(14, 30, by = 2)
  temperature_labels <- paste(
    head(temperature_breaks, -1L), tail(temperature_breaks, -1L),
    sep = "-"
  )
  temperature_colours <- rev(RColorBrewer::brewer.pal(
    length(temperature_labels), "Spectral"
  ))
  rate_labels <- c("<0.1", "0.1-0.5", "0.5-1", "1-2", "2-4", ">=4")
  rate_legend_labels <- rate_labels
  rate_colours <- rev(RColorBrewer::brewer.pal(
    length(rate_labels), "Spectral"
  ))
  shape$RateCategory <- cut(
    shape$K550CrudeRatePer100kPY,
    breaks = c(-Inf, .1, .5, 1, 2, 4, Inf),
    labels = rate_labels, right = FALSE
  )

  map_xlim <- c(-75.0, -35.3)
  map_ylim <- c(-34.8, 4.7)
  longitude_breaks <- c(-70, -60, -50, -40)
  latitude_breaks <- c(-30, -20, -10, 0)
  longitude_labels <- paste0(abs(longitude_breaks), "\u00b0W")
  latitude_labels <- ifelse(
    latitude_breaks < 0,
    paste0(abs(latitude_breaks), "\u00b0S"),
    paste0(latitude_breaks, "\u00b0")
  )
  region_labels <- data.frame(
    region = c("North", "Northeast", "Central-West", "Southeast", "South"),
    x = c(-58.7, -41.2, -53.5, -45.5, -52.1),
    y = c(-4.2, -9.7, -15.7, -19.8, -28.1)
  )

  draw_halo_text <- function(x, y, labels, cex = 7.2 / 12) {
    offsets <- expand.grid(dx = c(-.055, 0, .055), dy = c(-.055, 0, .055))
    offsets <- offsets[offsets$dx != 0 | offsets$dy != 0, ]
    for (i in seq_len(nrow(offsets))) {
      text(
        x + offsets$dx[i], y + offsets$dy[i], labels,
        cex = cex, col = "white", font = 2
      )
    }
    text(x, y, labels, cex = cex, col = "grey15", font = 2)
  }

  draw_north_arrow <- function(x, y, height = 2.4, width = .60) {
    polygon(
      c(x, x - width, x), c(y + height, y, y + .62),
      col = "black", border = "black", lwd = .45
    )
    polygon(
      c(x, x + width, x), c(y + height, y, y + .62),
      col = "white", border = "black", lwd = .45
    )
    text(x, y + height + .45, "N", cex = 7 / 12, font = 2)
  }

  add_cartographic_elements <- function(panel_tag) {
    axis(
      1,
      at = longitude_breaks, labels = longitude_labels,
      cex.axis = 7.2 / 12, lwd = .45, lwd.ticks = .35,
      tck = -.018, mgp = c(1.2, .32, 0), padj = -.2
    )
    axis(
      2,
      at = latitude_breaks, labels = latitude_labels,
      cex.axis = 7.2 / 12, las = 1, lwd = .45, lwd.ticks = .35,
      tck = -.018, mgp = c(1.2, .42, 0), hadj = .65
    )
    box(col = "black", lwd = .55)
    draw_halo_text(
      region_labels$x, region_labels$y, region_labels$region
    )
    draw_north_arrow(
      x = -71.5,
      y = 3.2
    )
    terra::sbar(
      d = 1000, xy = c(-72.5, -33.0), type = "bar", divs = 4,
      below = "km", lonlat = TRUE, cex = 6.4 / 12,
      adj = c(.5, -1.0), halo = TRUE, border = "black", xpd = FALSE
    )
    mtext(panel_tag, side = 3, line = .22, adj = 0, cex = .85, font = 2)
  }

  draw_map <- function(filename, device = c("png", "pdf")) {
    device <- match.arg(device)
    if (device == "png") {
      ragg::agg_png(
        filename,
        width = full_width_in, height = 3.28, units = "in",
        res = 600, background = "white"
      )
    } else {
      grDevices::pdf(
        filename,
        width = full_width_in, height = 3.28,
        family = "Helvetica", useDingbats = FALSE, bg = "white"
      )
    }
    old <- par(
      oma = c(.05, 0, .05, 0), family = "sans", xaxs = "i", yaxs = "i"
    )
    # Keep each map close to square while reserving compact legend columns.
    layout(matrix(1:4, nrow = 1L), widths = c(.90, .24, .90, .30))
    on.exit({
      par(old)
      dev.off()
    })

    par(mar = c(1.45, 1.55, 1.15, .08))
    plot(
      shape, "MeanTemperatureC",
      type = "interval",
      breaks = temperature_breaks, col = temperature_colours,
      border = NA, axes = FALSE, legend = FALSE, main = "",
      xlim = map_xlim, ylim = map_ylim
    )
    abline(
      v = longitude_breaks, h = latitude_breaks,
      col = adjustcolor("grey55", alpha.f = .30), lty = 3, lwd = .35
    )
    plot(water_shape, add = TRUE, col = "grey92", border = "grey70", lwd = .2)
    plot(country_shape, add = TRUE, col = NA, border = "grey25", lwd = .35)
    plot(region_shape, add = TRUE, col = NA, border = "black", lwd = .75)
    add_cartographic_elements("a   Annual mean temperature")

    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend(
      "center",
      legend = temperature_labels, fill = temperature_colours,
      title = expression(paste("Temperature (", degree, "C)")),
      bty = "n", cex = 10 / 12, border = NA,
      x.intersp = .55, y.intersp = .88
    )

    par(mar = c(1.45, 1.55, 1.15, .08))
    plot(
      shape,
      col = rate_colours[as.integer(shape$RateCategory)],
      border = NA, axes = FALSE, main = "",
      xlim = map_xlim, ylim = map_ylim
    )
    abline(
      v = longitude_breaks, h = latitude_breaks,
      col = adjustcolor("grey55", alpha.f = .30), lty = 3, lwd = .35
    )
    plot(water_shape, add = TRUE, col = "grey92", border = "grey70", lwd = .2)
    plot(country_shape, add = TRUE, col = NA, border = "grey25", lwd = .35)
    plot(region_shape, add = TRUE, col = NA, border = "black", lwd = .75)
    add_cartographic_elements("b   Crude mortality rate")  # Fresh full-panel composition supplies its own label.

    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend(
      "center",
      legend = rate_legend_labels, fill = rate_colours,
      title = "Average annual\ndeaths per 100,000\npersons", bty = "n", cex = 10 / 12,
      border = NA, x.intersp = .55, y.intersp = .88
    )
  }

  draw_map(filename, "pdf")
}

# Draw the monthly panels and return their replacement geometry.
draw_temporal_panels <- function(data, filename) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(grid)
    library(gtable)
    library(ragg)
  })

  t <- copy(data)
  t[,date:=as.Date(date)]

  blocks <- rbindlist(lapply(2000:2019, function(y) {
    data.table(start = as.Date(sprintf("%d-%s", y, c(
      "01-01",
      "03-01", "06-01", "09-01", "12-01"
    ))), end = as.Date(c(sprintf("%d-%s", y, c(
      "03-01", "06-01", "09-01",
      "12-01"
    )), sprintf("%d-01-01", y + 1))), season = factor(c(
      "Summer", "Autumn", "Winter", "Spring",
      "Summer"
    ), levels = c("Summer", "Autumn", "Winter", "Spring")))
  }))

  season_cols <- c(Summer = "#F6D8CC", Autumn = "#EFE7BD", Winter = "#D6E7F4", Spring = "#DCEEDB")

  base <- theme_classic(base_size = 9, base_family = "sans") + theme(
    axis.line = element_line(
      colour = "#7A8593",
      linewidth = 0.35
    ), axis.ticks = element_line(colour = "#7A8593", linewidth = 0.3), axis.text = element_text(
      size = 8,
      colour = "#20262C"
    ), axis.title = element_text(size = 9), panel.grid.major.y = element_line(
      colour = "#EDF0F3",
      linewidth = 0.25
    ), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(
      size = 9.5,
      face = "bold", hjust = 0, margin = margin(b = 7)
    ), plot.title.position = "panel", plot.margin = margin(
      9,
      8, 5, 6
    ), legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 8, lineheight = 0.95),
    legend.key.size = unit(11, "pt"), legend.key.spacing.y = unit(4, "pt"), legend.spacing.y = unit(3, "pt"), legend.margin = margin(
      0, 0,
      0, 10
    ), legend.box.spacing = unit(8, "pt")
  )

  make <- function(column, title, ytitle, col, showx, legend) {
    ggplot(t, aes(date, .data[[column]])) +
      geom_rect(data = blocks, aes(
        xmin = start, xmax = end, ymin = -Inf,
        ymax = Inf, fill = season
      ), inherit.aes = FALSE, alpha = 0.5, show.legend = legend) +
      geom_line(
        colour = col,
        linewidth = 0.5
      ) +
      scale_fill_manual(values = season_cols,
        labels = c(Summer = "Summer\n(Dec-Feb)", Autumn = "Autumn\n(Mar-May)",
                   Winter = "Winter\n(Jun-Aug)", Spring = "Spring\n(Sep-Nov)"),
        guide = guide_legend(ncol = 1, override.aes = list(alpha = 0.5))) +
      scale_x_date(limits = as.Date(c("2000-01-01", "2020-01-01")), breaks = as.Date(sprintf(
        "%d-01-01",
        seq(2000, 2018, 2)
      )), date_labels = "%Y", expand = c(0, 0)) +
      scale_y_continuous(expand = expansion(mult = c(
        0.04,
        0.07
      ))) +
      labs(title = title, x = if (showx) {
        "Calendar year"
      } else {
        NULL
      }, y = ytitle) +
      base +
      theme(axis.title.y = element_text(size = 8), axis.text.x = if (showx) {
        element_text()
      } else {
        element_blank()
      }, axis.ticks.x = if (showx) {
        element_line()
      } else {
        element_blank()
      })
  }

  C <- make(
    "mean_temperature_c", "c  Monthly mean temperature", "Temperature (°C)", "#C4473A", FALSE,
    TRUE
  )

  D <- make(
    "annualized_monthly_rate_per_100k", "d  Monthly mortality rate", "Annualised mortality rate\n(per 100,000 persons)",
    "#2B6CB0", TRUE, FALSE
  )

  shown <- range(t$displayed_rate_per_100k, na.rm = TRUE)
  new_yrange <- shown + c(-.04, .07) * diff(shown)
  mortality_breaks <- seq(1.1, 1.6, by = .1)

  D_revised <- make(
    "displayed_rate_per_100k", "d  Monthly mortality rate", "Annualised mortality rate\n(per 100,000 persons)",
    "#2B6CB0", TRUE, FALSE
  ) + scale_y_continuous(limits = new_yrange, breaks = mortality_breaks, expand = c(
    0,
    0
  ))

  pdf(NULL, width = 7.48, height = 3.2, family = "Helvetica", useDingbats = FALSE)

  gc <- ggplotGrob(C)
  gd <- ggplotGrob(D)
  newgd <- ggplotGrob(D_revised)

  for (component in c("panel", "axis-l")) {
    gd$grobs[[which(gd$layout$name == component)]] <- newgd$grobs[[which(newgd$layout$name ==
      component)]]
  }

  shared <- unit.pmax(gc$widths, gd$widths)

  gc$widths <- shared

  gd$widths <- shared

  pc <- gc$layout[gc$layout$name == "panel", ]
  pd <- gd$layout[gd$layout$name == "panel", ]

  fixedC <- convertHeight(sum(gc$heights), "bigpts", valueOnly = TRUE)
  fixedD <- convertHeight(sum(gd$heights), "bigpts", valueOnly = TRUE)
  panelH <- (3.2 * 72 - fixedC - fixedD) / 2

  gc$heights[pc$t] <- unit(panelH, "bigpts")

  gd$heights[pd$t] <- unit(panelH, "bigpts")

  heightC <- convertHeight(sum(gc$heights), "bigpts", valueOnly = TRUE)
  heightD <- convertHeight(sum(gd$heights), "bigpts", valueOnly = TRUE)
  dev.off()

  pdf(filename, width = 7.48, height = 3.2, family = "Helvetica", useDingbats = FALSE)

  grid.newpage()

  pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(heightC, heightD), "bigpts"))))

  pushViewport(viewport(layout.pos.row = 1))

  grid.draw(gc)

  upViewport()

  pushViewport(viewport(layout.pos.row = 2))

  grid.draw(gd)

  upViewport(2)

  dev.off()
}

compose_figure1 <- function(spatial_pdf, temporal_pdf, target) {
  gs <- Sys.which("gs")
  a<-pdf_page_size(spatial_pdf); b<-pdf_page_size(temporal_pdf)
  width<-max(a[1],b[1]);height<-a[2]+b[2]
  ps_string<-function(path){path<-normalizePath(path,mustWork=TRUE);path<-gsub("\\","\\\\",path,fixed=TRUE);path<-gsub("(","\\(",path,fixed=TRUE);gsub(")","\\)",path,fixed=TRUE)}
  script<-c(sprintf("<< /PageSize [%.12f %.12f] /PageHasTransparency true >> setpagedevice",width,height),
    "/spatial <<>> .PDFInit def",sprintf("(%s) spatial .PDFFile",ps_string(spatial_pdf)),
    "/temporal <<>> .PDFInit def",sprintf("(%s) temporal .PDFFile",ps_string(temporal_pdf)),
    sprintf("gsave 0 %.12f translate spatial 0 .PDFDrawPage grestore",b[2]),
    "temporal 0 .PDFDrawPage showpage spatial .PDFClose temporal .PDFClose")
  tmp<-tempfile(fileext=".ps");writeLines(script,tmp)
  permitted<-paste(normalizePath(c(spatial_pdf,temporal_pdf)),collapse=":")
  result<-system2(gs,c("-q","-dBATCH","-dNOPAUSE",paste0("--permit-file-read=",shQuote(permitted)),"-sDEVICE=pdfwrite","-dAutoRotatePages=/None",paste0("-sOutputFile=",shQuote(target)),shQuote(tmp)),stdout=TRUE,stderr=TRUE)
  status<-attr(result,"status");if(!is.null(status)&&status!=0L)stop(paste(result,collapse="\n"))
  stopifnot(file.exists(target))
  render_pdf_png(target,figure_file("Figure1","png"),600)
}

draw_figure3 <- function(data, selected_module, output_name) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(ragg)
  })
  d <- copy(data)

  groups <- data.table(
    subgroup = c(
      "Sex", "Age (years)", "Education (years)", "Race", "Per-capita income",
      "Hospital-bed density", "Air-conditioning prevalence", "Automobile prevalence"
    ),
    title = c(
      "Sex", "Age", "Education", "Race", "Per-capita income", "Hospital-bed density",
      "Household air-conditioning\nprevalence", "Household automobile\nownership"
    ),
    first = c("Male", "<65", "<8 years", "White", rep("Low", 4)),
    second = c("Female", ">=65", ">=8 years", "non-White", rep("High", 4)),
    module = c(rep("Individual", 4), rep("Municipality", 4)),
    header_y = c(54, 47.2, 40.4, 33.6, 22.2, 15.4, 8.6, 1.8),
    halfheight = c(rep(.55, 6), .85, .85)
  )
  groups <- groups[module == selected_module]
  groups[, header_y := c(25, 18.2, 11.4, 4.6)]
  d <- d[subgroup %in% groups$subgroup]
  d <- merge(d, groups, by = "subgroup", all.x = TRUE, sort = FALSE)
  d[, level_order := fifelse(level == first, 1L, 2L)]
  stopifnot(all(d$level %in% c(groups$first, groups$second)))
  d[, y := header_y - fifelse(level_order == 1L, 2, 4.2)]
  d[, label := fcase(level == "<65", "<65 years", level == ">=65", "65+ years",
    level == ">=8 years", "8+ years", level == "non-White", "Non-White",
    default = level
  )]
  d[, OR_text := sprintf("%.3f (%.3f, %.3f)", or, lower95, upper95)]
  d[, module_order := fifelse(module == "Individual", 1L, 2L)]
  d[, group_order := match(subgroup, groups$subgroup)]
  setorder(d, module_order, group_order, level_order, exposure)

  # One difference test per variable and exposure, placed on the VARIABLE heading.
  p_data <- d[,
    {
      p <- if (module[1] == "Individual") comparison_p[is.finite(comparison_p)] else unique(contextual_comparison_p)
      stopifnot(length(p) == 1L, is.finite(p))
      list(p_value = p, header_y = header_y[1], module = module[1])
    },
    by = .(subgroup, exposure)
  ]
  p_data[, p_text := sprintf("P for difference = %.3f", p_value)]
  rows <- unique(d[, .(subgroup, level, label, cases, y, module, group_order, level_order)])
  rows[, percent := 100 * cases / 52565]
  rows[, count_label := sprintf("N=%s (%.1f%%)", format(cases, big.mark = ",", trim = TRUE), percent)]
  modules <- data.table(
    y = 28.3,
    label = if (selected_module == "Individual") "Individual characteristics" else "Municipality-level\ncharacteristics"
  )
  yrange <- c(-1.5, 33.65)
  xrange <- c(.70, 2.15)
  breaks <- c(.75, 1, 1.5, 2)
  colors <- c(Cold = "#2B6CB0", Heat = "#C4473A")
  shapes <- c(Cold = 16, Heat = 17)
  base <- theme_bw(base_size = 8.5, base_family = "sans") +
    theme(
      panel.border = element_blank(), panel.grid = element_blank(),
      axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 8.5),
      axis.line.x = element_line(color = "#7A8593", linewidth = .35),
      axis.ticks.x = element_line(color = "#7A8593", linewidth = .35),
      plot.margin = margin(4, 4, 4, 4)
    )
  blank <- theme(
    axis.text.x = element_blank(), axis.title.x = element_blank(),
    axis.ticks.x = element_blank(), axis.line.x = element_blank()
  )
  ys <- scale_y_continuous(limits = yrange, expand = expansion(mult = 0))
  divider <- geom_hline(data = groups, aes(yintercept = header_y), inherit.aes = FALSE, color = "#D3DCE5", linewidth = .28)
  top_rule <- geom_hline(yintercept = 30, color = "#4A5568", linewidth = .45)
  labels <- ggplot(rows, aes(y = y)) +
    geom_segment(
      data = groups, aes(y = header_y - halfheight, yend = header_y + halfheight),
      x = 0, xend = 0, inherit.aes = FALSE, color = "#7C91AA", linewidth = .7
    ) +
    geom_text(aes(x = .04, label = label), hjust = 0, size = 2.85, family = "sans") +
    geom_text(aes(x = .04, y = y - .82, label = count_label),
      hjust = 0, size = 2.55, family = "sans", color = "#647181"
    ) +
    geom_text(
      data = groups, aes(x = .04, y = header_y, label = title), inherit.aes = FALSE,
      hjust = 0, size = 2.8, fontface = "bold", family = "sans", color = "#34465C", lineheight = .95
    ) +
    geom_text(
      data = modules, aes(x = .04, y = y, label = label), inherit.aes = FALSE,
      hjust = 0, size = 2.9, fontface = "bold", family = "sans", color = "#263B53", lineheight = .95
    ) +
    top_rule +
    annotate("text", x = .04, y = 32.15, label = "Subgroup", hjust = 0, size = 3.0, fontface = "bold", family = "sans") +
    annotate("text", x = .04, y = 30.95, label = "Deaths, N (%)", hjust = 0, size = 2.45, family = "sans", color = "#647181") +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
    ys +
    base +
    blank

  forest <- function(which_exposure, panel) {
    individual <- selected_module == "Individual"
    header_label <- if (individual) c(Cold = "Low temperature", Heat = "High temperature")[[which_exposure]] else which_exposure
    header_size <- if (individual) 2.50 else 3.15
    header_x <- if (individual) 1.32 else sqrt(prod(xrange))
    if (individual) panel <- tolower(panel)
    dd <- d[exposure == which_exposure]
    ggplot(dd, aes(x = or, y = y)) +
      geom_segment(
        data = data.table(x = breaks), aes(x = x, xend = x), y = yrange[1], yend = 30,
        inherit.aes = FALSE, color = "#EDF0F3", linewidth = .25
      ) +
      divider +
      annotate("segment", x = 1, xend = 1, y = yrange[1], yend = 30, linetype = 2, color = "#8793A2", linewidth = .4) +
      top_rule +
      geom_errorbar(aes(xmin = lower95, xmax = upper95),
        orientation = "y", width = .25,
        linewidth = .55, color = colors[[which_exposure]]
      ) +
      geom_point(shape = shapes[[which_exposure]], size = 1.95, color = colors[[which_exposure]]) +
      annotate("text", x = .725, y = 31.55, label = panel, hjust = 0, size = 3.15, fontface = "bold", family = "sans") +
      annotate("text",
        x = header_x, y = 31.55, label = header_label, size = header_size, fontface = "bold",
        family = "sans", color = colors[[which_exposure]]
      ) +
      scale_x_log10(limits = xrange, breaks = breaks, labels = sprintf("%.2f", breaks), expand = expansion(mult = 0)) +
      ys +
      labs(x = "Odds ratio") +
      base +
      theme(axis.text.x = element_text(size = 7.5))
  }
  values <- function(which_exposure) {
    dd <- d[exposure == which_exposure]
    ggplot(dd, aes(y = y)) +
      top_rule +
      geom_text(aes(x = .02, label = OR_text), hjust = 0, size = 2.7, family = "sans", color = colors[[which_exposure]]) +
      annotate("text", x = .02, y = 31.55, label = "OR (95% CI)", hjust = 0, size = 2.9, fontface = "bold", family = "sans") +
      scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
      ys +
      base +
      blank
  }
  pcolumn <- function(which_exposure) {
    pp <- p_data[exposure == which_exposure]
    ggplot(pp) +
      top_rule +
      geom_text(aes(x = .5, y = header_y, label = sprintf("%.3f", p_value)),
        size = 2.55, family = "sans", color = "#526276"
      ) +
      annotate("text",
        x = .5, y = 31.55, label = "P for\ndifference", size = 2.5,
        fontface = "bold", family = "sans", color = "#526276", lineheight = .95
      ) +
      scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
      ys +
      base +
      blank
  }
  cold <- forest("Cold", "A")
  heat <- forest("Heat", "B")
  cold_values <- values("Cold")
  heat_values <- values("Heat")
  cold_p <- pcolumn("Cold")
  heat_p <- pcolumn("Heat")
  p <- (labels + cold + cold_values + cold_p + heat + heat_values + heat_p) +
    plot_layout(widths = c(2.5, 1.5, 1.7, .8, 1.5, 1.7, .8)) +
    plot_annotation(theme = theme(plot.margin = margin(4, 4, 4, 4)))
  figure_height <- 4.8
  figure_height_pdf <- 4.8
  prefix <- file.path(out_dir, output_name)
  ggsave(paste0(prefix, ".png"), p, width = 7.48, height = figure_height, dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(paste0(prefix, ".pdf"), p,
    width = 7.48, height = figure_height_pdf, device = grDevices::pdf,
    family = "Helvetica", useDingbats = FALSE, bg = "white"
  )
  stopifnot(all(d$lower95 >= xrange[1]), all(d$upper95 <= xrange[2]))

  invisible(p_data)
}

draw_figure_s1 <- function(data) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(ragg)
    library(scales)
  })
  plot_data <- copy(data)
  display_age_levels <- c("0-19", paste0(seq(20, 85, 5), "-", seq(24, 89, 5)), "90+")
  plot_data[, `:=`(sex = factor(sex, levels = c("Female", "Male")), age_group = factor(age_group, levels = display_age_levels))]
  left_max <- max(plot_data$deaths_thousands, na.rm = TRUE)
  right_max <- max(plot_data$mortality_rate_per_100k, na.rm = TRUE)
  rate_to_count_scale <- (left_max * 1.08) / right_max

  plot_data[, x_position := as.integer(age_group) + fifelse(
    sex == "Female", -0.18, 0.18
  )]

  sex_colours <- c(Female = "#2166AC", Male = "#B2182B")
  p <- ggplot(plot_data, aes(x = x_position)) +
    geom_col(
      aes(y = deaths_thousands, fill = sex),
      width = 0.33, alpha = 0.30, colour = NA
    ) +
    geom_line(
      aes(
        y = mortality_rate_per_100k * rate_to_count_scale,
        colour = sex, group = sex
      ),
      linewidth = 0.72
    ) +
    geom_point(
      aes(
        y = mortality_rate_per_100k * rate_to_count_scale,
        colour = sex
      ),
      shape = 21, fill = "white", stroke = 0.65, size = 1.65
    ) +
    scale_x_continuous(
      breaks = seq_along(display_age_levels),
      labels = display_age_levels,
      expand = expansion(mult = c(0.015, 0.025))
    ) +
    scale_y_continuous(
      name = "Deaths (1000)",
      labels = label_number(accuracy = 1),
      expand = expansion(mult = c(0, 0.10)),
      sec.axis = sec_axis(
        ~ . / rate_to_count_scale,
        name = "Average annual mortality rate (per 100 000 persons)",
        labels = label_number(accuracy = 0.1)
      )
    ) +
    scale_fill_manual(values = sex_colours, name = "Sex") +
    scale_colour_manual(values = sex_colours, name = "Sex") +
    guides(
      fill = guide_legend(order = 1, override.aes = list(alpha = 0.30)),
      colour = "none"
    ) +
    labs(
      x = "Age group (years)",
      caption = NULL
    ) +
    theme_bw(base_size = 10, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(
        colour = "grey82", linewidth = 0.35, linetype = "dotted"
      ),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.55),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8.2),
      axis.text = element_text(colour = "black", size = 9.5),
      axis.title = element_text(colour = "black", size = 10),
      axis.title.y.left = element_text(margin = margin(r = 8)),
      axis.title.y.right = element_text(margin = margin(l = 8)),
      legend.position = c(0.075, 0.90),
      legend.justification = c(0, 1),
      legend.background = element_rect(fill = alpha("white", 0.88), colour = "grey60"),
      legend.key.height = grid::unit(0.38, "cm"),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9.5),
      plot.margin = margin(8, 12, 8, 8)
    )

  png_file <- figure_file("FigureS1", "png")
  pdf_file <- figure_file("FigureS1", "pdf")
  save_plot(
    png_file, p,
    width = 7.48, height = 4.80, dpi = 600,
    bg = "white", device = ragg::agg_png
  )
  save_plot(
    pdf_file, p,
    width = 7.48, height = 4.80,
    bg = "white", device = pdf, family = "Helvetica", useDingbats = FALSE
  )
}

draw_figure_s2 <- function(data) {
  suppressPackageStartupMessages(library(ragg))
  surface <- copy(data)
  required <- c("temperature_c", "lag", "log_or", "or", "reference_c")
  if (!all(required %chin% names(surface))) stop("Unexpected surface data structure")
  if (anyNA(surface[, ..required])) stop("Surface data contain missing values")

  reference_temperature <- unique(surface$reference_c)
  if (length(reference_temperature) != 1L) stop("Expected one reference temperature")

  temperatures <- sort(unique(surface$temperature_c))
  lags <- sort(unique(surface$lag))
  or_wide <- dcast(surface, temperature_c ~ lag, value.var = "or")
  log_or_wide <- dcast(surface, temperature_c ~ lag, value.var = "log_or")
  setorder(or_wide, temperature_c)
  setorder(log_or_wide, temperature_c)

  # The prediction grid is much denser than is useful in a printed 3D figure.
  # Retain enough temperature slices to preserve the fitted surface while avoiding
  # a dark, moire-like mesh in high-resolution output.
  temperature_index <- unique(sort(c(
    round(seq(1L, nrow(or_wide), length.out = 76L)),
    which.min(abs(or_wide$temperature_c - reference_temperature))
  )))
  or_wide <- or_wide[temperature_index]
  log_or_wide <- log_or_wide[temperature_index]
  temperatures <- or_wide$temperature_c
  or_matrix <- as.matrix(or_wide[, -"temperature_c"])
  log_or_matrix <- as.matrix(log_or_wide[, -"temperature_c"])
  if (!isTRUE(all.equal(as.numeric(names(or_wide)[-1L]), as.numeric(lags)))) {
    stop("Lag grid is incomplete")
  }

  # Colour each surface facet according to the mean log odds ratio at its corners.
  facet_log_or <- (
    log_or_matrix[-nrow(log_or_matrix), -ncol(log_or_matrix)] +
      log_or_matrix[-1L, -ncol(log_or_matrix)] +
      log_or_matrix[-nrow(log_or_matrix), -1L] +
      log_or_matrix[-1L, -1L]
  ) / 4

  colour_limit <- max(abs(surface$log_or))
  colour_values <- seq(-colour_limit, colour_limit, length.out = 256L)
  surface_palette <- grDevices::colorRampPalette(
    c("#2B6CB0", "#F7F7F7", "#C33A3A"),
    space = "Lab"
  )(length(colour_values))
  colour_index <- round(
    (facet_log_or + colour_limit) / (2 * colour_limit) *
      (length(surface_palette) - 1L)
  ) + 1L
  colour_index <- pmax(1L, pmin(length(surface_palette), colour_index))
  facet_colours <- surface_palette[colour_index]
  or_range <- range(or_matrix)
  z_limits <- c(
    floor((or_range[1L] - 0.01) * 20) / 20,
    ceiling((or_range[2L] + 0.01) * 20) / 20
  )
  legend_breaks <- c(0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15, 1.20)
  legend_breaks <- legend_breaks[
    log(legend_breaks) >= -colour_limit & log(legend_breaks) <= colour_limit
  ]

  draw_colour_bar <- function() {
    par(mar = c(4.7, 0.0, 2.5, 2.8), family = "sans", xaxs = "i", yaxs = "i")
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(-colour_limit, colour_limit))
    lower <- head(colour_values, -1L)
    upper <- tail(colour_values, -1L)
    rect(
      xleft = 0.18, ybottom = lower, xright = 0.56, ytop = upper,
      col = surface_palette[-length(surface_palette)], border = NA
    )
    rect(
      xleft = 0.18, ybottom = -colour_limit,
      xright = 0.56, ytop = colour_limit,
      border = "black", lwd = 0.8
    )
    axis(
      side = 4, at = log(legend_breaks), labels = sprintf("%.2f", legend_breaks),
      las = 1, line = -0.48, tcl = -0.20, cex.axis = 0.76
    )
    mtext("Odds ratio", side = 3, line = 0.48, at = 0.37, cex = 0.80)
  }

  draw_surface <- function() {
    layout(matrix(c(1L, 2L), nrow = 1L), widths = c(8.9, 1.1))

    par(
      mar = c(3.5, 3.2, 0.35, 0.05), family = "sans",
      mgp = c(1.85, 0.55, 0), tcl = -0.22,
      cex.axis = 0.78, cex.lab = 0.84
    )
    projection <- persp(
      x = temperatures,
      y = lags,
      z = or_matrix,
      xlim = range(temperatures),
      ylim = range(lags),
      zlim = z_limits,
      theta = -47,
      phi = 28,
      expand = 0.62,
      col = as.vector(facet_colours),
      border = grDevices::adjustcolor("grey30", alpha.f = 0.22),
      lwd = 0.22,
      shade = 0.05,
      ticktype = "detailed",
      nticks = 5,
      xlab = "Daily mean temperature (\u00b0C)",
      ylab = "Lag (days)",
      zlab = "Odds ratio"
    )

    reference_line <- trans3d(
      x = rep(reference_temperature, length(lags)),
      y = lags,
      z = rep(1, length(lags)),
      pmat = projection
    )
    lines(reference_line,
      col = grDevices::adjustcolor("black", alpha.f = 0.82),
      lty = 3, lwd = 0.90
    )

    draw_colour_bar()
    layout(1)
  }

  png_file <- figure_file("FigureS2", "png")
  pdf_file <- figure_file("FigureS2", "pdf")
  ragg::agg_png(
    png_file,
    width = 7.48, height = 5.35, units = "in", res = 600,
    background = "white"
  )
  draw_surface()
  dev.off()

  grDevices::pdf(
    pdf_file,
    width = 7.48, height = 5.35,
    family = "Helvetica", useDingbats = FALSE
  )
  draw_surface()
  dev.off()
}

draw_figure_s3 <- function(data) {
  suppressPackageStartupMessages({
    library(terra)
    library(ragg)
  })
  cutoffs <- data$cutoffs[, .(modifier, cutoff)]
  cutoff_for <- function(id) {
    ans <- cutoffs[modifier == id, unique(cutoff)]
    if (length(ans) != 1L || !is.finite(ans)) {
      stop("Missing or non-unique subgroup cutoff for ", id)
    }
    ans
  }
  cutoff_label_for <- function(id) {
    formatC(
      cutoff_for(id),
      format = "f", digits = 2L
    )
  }

  municipal_values <- copy(data$values)
  specs <- list(
    income_per_capita = list(
      tag = "A",
      modifier = "economic",
      title = "Per-capita income",
      unit = "BRL (2010)",
      transform = "log1p", tick_digits = 0L,
      cutoff_label = cutoff_label_for("economic")
    ),
    hospital_beds = list(
      tag = "B",
      modifier = "hospital_beds",
      title = "Hospital-bed density",
      unit = "per 1000 residents (2002)",
      transform = "log1p", tick_digits = 1L,
      cutoff_label = cutoff_label_for("hospital_beds")
    ),
    air_conditioning = list(
      tag = "C",
      modifier = "air_conditioning",
      title = "Household air conditioning",
      unit = "% of households (2000)",
      transform = "log1p", tick_digits = 1L,
      cutoff_label = cutoff_label_for("air_conditioning")
    ),
    automobile = list(
      tag = "D",
      modifier = "transport",
      title = "Household automobile ownership",
      unit = "% of households (2000)",
      transform = "identity", tick_digits = 0L,
      cutoff_label = cutoff_label_for("transport")
    )
  )

  # Use the same municipality geometry and cartographic framing as Figure 1.
  shape <- read_boundaries("municipalities")
  water_shape <- read_boundaries("water")
  region_shape <- read_boundaries("regions")
  country_shape <- read_boundaries("country")
  low_colours <- grDevices::colorRampPalette(c("#08519C", "#FFFFFF"))(256L)
  high_colours <- grDevices::colorRampPalette(c("#FFFFFF", "#A63603"))(256L)
  missing_colour <- "grey90"
  transform_values <- function(x, method) {
    if (method == "log1p") log1p(x) else x
  }

  inverse_transform <- function(x, method) {
    if (method == "log1p") expm1(x) else x
  }

  continuous_map_colours <- function(values, cutoff, method) {
    ans <- rep(missing_colour, length(values))
    valid <- is.finite(values)
    if (!any(valid)) {
      return(ans)
    }

    low <- valid & values <= cutoff
    high <- valid & values > cutoff
    low_min <- min(values[low], na.rm = TRUE)
    high_max <- max(values[high], na.rm = TRUE)
    transformed_cutoff <- transform_values(cutoff, method)

    if (any(low)) {
      z <- transform_values(values[low], method)
      z0 <- transform_values(low_min, method)
      score <- (z - z0) / (transformed_cutoff - z0)
      idx <- pmax(1L, pmin(256L, 1L + floor(score * 255)))
      ans[low] <- low_colours[idx]
    }
    if (any(high)) {
      z <- transform_values(values[high], method)
      z1 <- transform_values(high_max, method)
      score <- (z - transformed_cutoff) / (z1 - transformed_cutoff)
      idx <- pmax(1L, pmin(256L, 1L + floor(score * 255)))
      ans[high] <- high_colours[idx]
    }
    ans
  }

  # Shift the mapped geography about 1.8% of the panel width to the right while
  # retaining the same map scale. Header elements and legends use panel-relative
  # coordinates and therefore remain fixed.
  map_xlim <- c(-75.7, -36.0)
  # Reserve a compact cartographic strip above the northernmost municipalities
  # for the tag, panel label, and north arrow. A near-square map panel prevents
  # the plotting device from adding a large automatic upper margin.
  map_ylim <- c(-34.8, 7.3)
  longitude_breaks <- c(-70, -60, -50, -40)
  latitude_breaks <- c(-30, -20, -10, 0)
  longitude_labels <- paste0(abs(longitude_breaks), "°W")
  latitude_labels <- ifelse(
    latitude_breaks < 0,
    paste0(abs(latitude_breaks), "°S"),
    paste0(latitude_breaks, "°")
  )

  draw_north_arrow <- function(x, y, height = 2.25, width = .55) {
    polygon(
      c(x, x - width, x), c(y + height, y, y + .58),
      col = "black", border = "black", lwd = .4
    )
    polygon(
      c(x, x + width, x), c(y + height, y, y + .58),
      col = "white", border = "black", lwd = .4
    )
    text(x, y + height + .42, "N", cex = 7 / 12, font = 2)
  }

  draw_map_panel <- function(id) {
    sp <- specs[[id]]
    values <- municipal_values[[id]][match(shape$Code, municipal_values$Code)]
    polygon_colours <- continuous_map_colours(
      values, cutoff_for(sp$modifier), sp$transform
    )

    par(mar = c(1.35, 1.55, .18, .05))
    plot(
      shape,
      col = polygon_colours, border = NA, axes = FALSE, main = "",
      xlim = map_xlim, ylim = map_ylim
    )
    abline(
      v = longitude_breaks, h = latitude_breaks,
      col = adjustcolor("grey55", alpha.f = .30), lty = 3, lwd = .35
    )
    plot(water_shape, add = TRUE, col = "grey94", border = "grey72", lwd = .2)
    plot(country_shape, add = TRUE, col = NA, border = "grey25", lwd = .4)
    plot(region_shape, add = TRUE, col = NA, border = "black", lwd = .72)
    axis(
      1,
      at = longitude_breaks, labels = longitude_labels,
      cex.axis = 7 / 12, lwd = .45, lwd.ticks = .35,
      tck = -.018, mgp = c(1.2, .32, 0), padj = -.2
    )
    axis(
      2,
      at = latitude_breaks, labels = latitude_labels,
      cex.axis = 7 / 12, las = 1, lwd = .45, lwd.ticks = .35,
      tck = -.018, mgp = c(1.2, .42, 0), hadj = .65
    )
    box(col = "black", lwd = .55)

    # Match the relative placement in main-text Figure 1: tag and analytical
    # panel label at the upper-left, north arrow below them, and map centred in
    # the same longitude/latitude frame.
    usr <- par("usr")
    dx <- diff(usr[1:2])
    dy <- diff(usr[3:4])
    tag_x <- usr[1] + .045 * dx
    label_x <- usr[1] + .105 * dx
    # terra preserves the map aspect ratio inside the framed plotting region,
    # so place the header slightly above the returned map coordinates to keep
    # it close to the panel's upper border, as in main-text Figure 1.
    top_y <- usr[4] + .065 * dy
    text(
      tag_x, top_y, sp$tag,
      adj = c(0, 1), cex = 10 / 12, font = 2, xpd = NA
    )
    panel_label <- paste0(sp$title, " (", sp$unit, ")")
    text(
      label_x, top_y, panel_label,
      adj = c(0, 1), cex = 8.5 / 12, font = 2, xpd = NA
    )
    draw_north_arrow(
      x = usr[1] + .095 * dx,
      y = usr[4] - .130 * dy,
      height = .080 * dy,
      width = .017 * dx
    )
    terra::sbar(
      d = 1000, xy = c(-72.5, -33.0), type = "bar", divs = 4,
      below = "km", lonlat = TRUE, cex = 6.2 / 12,
      adj = c(.5, -1.0), halo = TRUE, border = "black", xpd = FALSE
    )
  }

  draw_legend_panel <- function(id) {
    sp <- specs[[id]]
    cutoff <- cutoff_for(sp$modifier)
    values <- municipal_values[[id]][match(shape$Code, municipal_values$Code)]
    low_min <- min(values[is.finite(values) & values <= cutoff], na.rm = TRUE)
    high_max <- max(values[is.finite(values) & values > cutoff], na.rm = TRUE)
    format_tick <- function(x) {
      if (abs(x - cutoff) < 1e-8) {
        return(sp$cutoff_label)
      }
      formatC(x, format = "f", digits = sp$tick_digits, big.mark = ",")
    }

    draw_vertical_gradient <- function(colours, x0, x1, y0, y1) {
      y_edges <- seq(y0, y1, length.out = length(colours) + 1L)
      rect(
        x0, head(y_edges, -1L), x1, tail(y_edges, -1L),
        col = colours, border = NA
      )
      rect(x0, y0, x1, y1, col = NA, border = "grey40", lwd = .35)
    }

    par(mar = c(.25, .15, .25, .10))
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(0, 1))

    # One vertical continuous bar: blue values below the population-weighted
    # median, white at the cutoff, and orange values above it.
    draw_vertical_gradient(
      c(low_colours, high_colours), .12, .36, .18, .82
    )
    segments(.12, .50, .36, .50, col = "grey25", lwd = .5)

    # Seven ticks give the continuous scale enough quantitative detail while
    # retaining the population-weighted median as the central reference tick.
    low_z <- seq(
      transform_values(low_min, sp$transform),
      transform_values(cutoff, sp$transform),
      length.out = 4L
    )
    high_z <- seq(
      transform_values(cutoff, sp$transform),
      transform_values(high_max, sp$transform),
      length.out = 4L
    )
    tick_values <- c(
      inverse_transform(low_z, sp$transform),
      inverse_transform(high_z[-1L], sp$transform)
    )
    tick_y <- c(seq(.18, .50, length.out = 4L), seq(.50, .82, length.out = 4L)[-1L])
    for (i in seq_along(tick_y)) {
      segments(.36, tick_y[i], .43, tick_y[i], lwd = .35)
      text(
        .48, tick_y[i], format_tick(tick_values[i]),
        adj = c(0, .5), cex = 7.6 / 12
      )
    }
  }

  draw_figure <- function(filename, device = c("png", "pdf")) {
    device <- match.arg(device)
    width <- 190 / 25.4
    height <- 158 / 25.4
    if (device == "png") {
      ragg::agg_png(
        filename,
        width = width, height = height, units = "in", res = 600,
        background = "white"
      )
    } else {
      grDevices::pdf(
        filename,
        width = width, height = height,
        family = "Helvetica", useDingbats = FALSE, bg = "white"
      )
    }

    old <- par(
      oma = c(.10, .05, .10, .05), family = "sans", xaxs = "i", yaxs = "i"
    )
    layout(
      matrix(1:8, nrow = 2L, byrow = TRUE),
      widths = c(1, .20, 1, .20), heights = c(1, 1)
    )
    on.exit({
      par(old)
      dev.off()
    })

    for (id in names(specs)) {
      draw_map_panel(id)
      draw_legend_panel(id)
    }
  }

  stem <- sub("[.]pdf$", "", figure_file("FigureS3", "pdf"))
  draw_figure(paste0(stem, ".png"), "png")
  draw_figure(paste0(stem, ".pdf"), "pdf")
}

draw_figure_s5 <- function(data) {
  suppressPackageStartupMessages(library(grid))
  categories <- c("Total non-optimal", "Cold", "Moderate cold", "Extreme cold", "Heat", "Moderate heat", "Extreme heat")
  sections <- list(
    list(group = "Sex", title = "Sex", levels = c("Male", "Female"), labels = c("Male", "Female")),
    list(group = "Age (years)", title = "Age", levels = c("<65", ">=65"), labels = c("<65 years", "65+ years")),
    list(group = "Education (years)", title = "Education", levels = c(">=8 years", "<8 years"), labels = c("8+ years", "<8 years")),
    list(group = "Race/skin colour", title = "Race / skin colour", levels = c("White", "non-White"), labels = c("White", "Non-White"))
  )
  dark <- "#202C38"
  muted <- "#617080"
  text_mm <- function(x, y, label, size = 8.1, bold = FALSE, colour = dark, align = "left") {
    grid.text(label, x = unit(x, "mm"), y = unit(y, "mm"), just = c(align, "bottom"), gp = gpar(fontfamily = "Arial", fontsize = size, fontface = if (bold) "bold" else "plain", col = colour))
  }
  rule_mm <- function(y, x0 = 5, x1 = 186, colour = "#AAB5BF", width = .4) {
    grid.lines(x = unit(c(x0, x1), "mm"), y = unit(c(y, y), "mm"), gp = gpar(col = colour, lwd = width / .75))
  }
  integer_text <- function(x) format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
  ci_text <- function(lower, upper, integer = FALSE) {
    values <- if (integer) integer_text(c(lower, upper)) else sprintf("%.2f", c(lower, upper))
    gsub("-", "−", paste0("(", values[1], ", ", values[2], ")"), fixed = TRUE)
  }

  target <- figure_file("FigureS5", "pdf")
  # Embed the system Arial font without requiring X11 on macOS.
  if (capabilities("aqua")) {
    grDevices::quartz(type = "pdf", file = target, width = 190 / 25.4, height = 247 / 25.4, family = "Arial", bg = "white")
  } else {
    grDevices::cairo_pdf(target, width = 190 / 25.4, height = 247 / 25.4, family = "Arial", bg = "white")
  }
  on.exit(dev.off(), add = TRUE)
  grid.newpage()
  for (section_index in seq_along(sections)) {
    section <- sections[[section_index]]
    top <- 240 - 60 * (section_index - 1L)
    first_row <- top - 23.7
    text_mm(5, top, LETTERS[section_index], size = 12.8, bold = TRUE)
    text_mm(13, top, section$title, size = 10.5, bold = TRUE)
    text_mm(5, top - 16.5, "Temperature category", size = 8.2, bold = TRUE)
    for (level_index in seq_along(section$levels)) {
      level_name <- section$levels[level_index]
      sample <- data[subgroup == section$group & level == level_name & scenario == "Total non-optimal"]
      center <- c(83, 152)[level_index]
      text_mm(center, top - 6, section$labels[level_index], size = 9.4, bold = TRUE, align = "center")
      text_mm(center, top - 10, paste0("N = ", integer_text(sample$deaths)), size = 7.8, colour = muted, align = "center")
      rule_mm(top - 12.7, x0 = c(50, 120)[level_index], x1 = c(116, 186)[level_index])
      text_mm(c(65, 134)[level_index], top - 16.5, "EF, % (95% CI)", size = 7.8, bold = TRUE, align = "center")
      text_mm(c(99, 168)[level_index], top - 16.5, "EN, n (95% CI)", size = 7.8, bold = TRUE, align = "center")
    }
    rule_mm(top - 19.2, width = .45, colour = dark)
    for (category_index in seq_along(categories)) {
      category <- categories[category_index]
      y <- first_row - (category_index - 1L) * 4.5
      parent <- category %in% c("Total non-optimal", "Cold", "Heat")
      label <- switch(category,
        "Total non-optimal" = "Non-optimal temperature",
        "Cold" = "Cold, total",
        "Heat" = "Heat, total",
        category
      )
      if (category == "Total non-optimal") grid.rect(x = unit(5, "mm"), y = unit(y - 1.4, "mm"), width = unit(181, "mm"), height = unit(4.4, "mm"), just = c("left", "bottom"), gp = gpar(fill = "#F4F6F8", col = NA))
      colour <- switch(category,
        "Cold" = "#2B6CB0",
        "Heat" = "#B83D32",
        dark
      )
      text_mm(if (parent) 5 else 8, y, label, size = 8.6, bold = parent, colour = colour)
      for (level_index in seq_along(section$levels)) {
        level_name <- section$levels[level_index]
        record <- data[subgroup == section$group & level == level_name & scenario == category]
        stopifnot(nrow(record) == 1L)
        ef <- sprintf("%.2f", record$attributable_fraction_percent)
        ef_ci <- ci_text(record$attributable_fraction_lower95, record$attributable_fraction_upper95)
        en <- integer_text(record$attributable_deaths)
        en_ci <- ci_text(record$attributable_deaths_lower95, record$attributable_deaths_upper95, integer = TRUE)
        text_mm(c(60, 129)[level_index], y, ef, size = 8.7, bold = TRUE, align = "right")
        text_mm(c(61.8, 130.8)[level_index], y, ef_ci, size = 8, colour = muted)
        text_mm(c(92.7, 161.7)[level_index], y, en, size = 8.7, bold = TRUE, align = "right")
        text_mm(c(94.5, 163.5)[level_index], y, en_ci, size = 8, colour = muted)
      }
    }
    rule_mm(top - 53.3, width = .45)
  }
  dev.off()
  on.exit(NULL)
  render_pdf_png(target, figure_file("FigureS5", "png"), 220)
}

# Figure 2. Cumulative and lag-specific associations, with the case-day temperature distribution.
draw_figure2 <- function() {
  curve <- fread(file.path(data_dir, "Figure2_curve.csv"))[order(TemperatureC)]
  bins <- fread(file.path(data_dir, "Figure2_histogram.csv"))
  lags <- fread(file.path(data_dir, "Figure2_lags.csv"))
  ref <- unique(curve$ReferenceC)
  stopifnot(length(ref)==1, identical(unique(curve$CumulativeLag), "0-5"),
            curve$TemperatureC[which.min(curve$OR)]==ref)
  colors <- list(low=c("#295870","#d7e4e9"), high=c("#934833","#eedbd3"))
  stroke <- function(col="#646766", lwd=.65, lty=1) grid::gpar(col=col,lwd=lwd,lty=lty,lineend="round")
  text_gp <- function(size, bold=FALSE) grid::gpar(col="#282a29",fontsize=size,fontfamily="sans",fontface=if(bold) "bold" else "plain")
  u <- grid::unit
  panel <- function(x,y,w,h,xrange,yrange) {
    grid::pushViewport(grid::viewport(x=x,y=y,width=w,height=h,just=c("left","bottom"),
                                      xscale=xrange,yscale=yrange,clip="off"))
  }
  heading <- function(tag,title) {
    y <- u(1,"npc")+u(10,"pt")
    grid::grid.text(tag,x=u(-6,"pt"),y=y,just=c("right","bottom"),gp=text_gp(11.5,TRUE))
    grid::grid.text(title,x=u(2,"pt"),y=y,just=c("left","bottom"),gp=text_gp(10))
  }
  axis_y <- function(values,labels,title,title_size=8.7) {
    grid::grid.segments(x0=0,x1=0,y0=0,y1=1,gp=stroke())
    y <- u(values,"native")
    grid::grid.segments(x0=u(-3,"pt"),x1=0,y0=y,y1=y,gp=stroke())
    grid::grid.text(labels,x=u(-5.5,"pt"),y=y,just="right",gp=text_gp(8))
    grid::grid.text(title,x=u(-27,"pt"),y=.5,rot=90,gp=text_gp(title_size))
  }
  axis_x <- function(values,label=NULL,tick_labels=TRUE) {
    x <- u(values,"native")
    grid::grid.segments(x0=0,x1=1,y0=0,y1=0,gp=stroke())
    grid::grid.segments(x0=x,x1=x,y0=0,y1=u(-3,"pt"),gp=stroke())
    if(tick_labels) grid::grid.text(values,x=x,y=u(-5,"pt"),just="top",gp=text_gp(8))
    if(!is.null(label)) grid::grid.text(label,x=.5,y=u(-19,"pt"),just="top",gp=text_gp(8.7))
  }
  reference_line <- function() grid::grid.segments(x0=0,x1=1,y0=u(1,"native"),y1=u(1,"native"),gp=stroke("#686b69",.7,"33"))
  draw_curve <- function(x,y,lo,hi,palette) {
    grid::grid.polygon(x=u(c(x,rev(x)),"native"),y=u(c(lo,rev(hi)),"native"),
                       gp=grid::gpar(col=NA,fill=palette[2]))
    grid::grid.lines(u(x,"native"),u(y,"native"),gp=stroke(palette[1],1.7))
  }
  draw_page <- function() {
    grid::grid.newpage()
    panel(.081,.345,.493,.52,c(10,30),c(.975,1.385))
    for (side in names(colors)) {
      d <- if(side=="low") curve[TemperatureC<=ref] else curve[TemperatureC>=ref]
      draw_curve(d$TemperatureC,d$OR,d$Lower95,d$Upper95,colors[[side]])
    }
    reference_line()
    axis_y(seq(1,1.3,.1),sprintf("%.1f",seq(1,1.3,.1)),"Cumulative odds ratio")
    heading("a","Cumulative association, lags 0–5")
    ends <- curve[c(1L,.N)]
    for(i in 1:2) {
      d <- ends[i]
      x <- u(d$TemperatureC,"native"); side <- c("left","right")[i]
      grid::grid.text(sprintf("%s   %.2f °C",c("P1","P99")[i],d$TemperatureC),x,u(1.372,"native"),just=c(side,"bottom"),gp=text_gp(8.3))
      grid::grid.text(sprintf("%.3f (%.3f–%.3f)",d$OR,d$Lower95,d$Upper95),x,u(1.345,"native"),just=c(side,"bottom"),gp=text_gp(7.5))
    }
    grid::grid.segments(x0=u(ref,"native"),x1=u(ref,"native"),y0=u(1.006,"native"),y1=u(1.089,"native"),gp=stroke())
    grid::grid.points(u(ref,"native"),u(1,"native"),pch=16,size=u(3,"pt"),gp=stroke("#282a29"))
    grid::grid.text(sprintf("MMT  %.2f °C",ref),u(ref,"native"),u(1.095,"native"),just="bottom",gp=text_gp(8))
    grid::popViewport()

    panel(.081,.145,.493,.155,c(10,30),c(0,6.3))
    grid::grid.rect(x=u(bins$xmin,"native"),y=u(0,"native"),
                    width=u(bins$xmax-bins$xmin,"native"),height=u(100*bins$proportion,"native"),
                    just=c("left","bottom"),gp=grid::gpar(fill="#adb1af",col="white",lwd=.34))
    axis_y(c(0,3,6),c("0","3","6"),"Case days (%)",8)
    axis_x(seq(10,30,5),"Daily mean temperature (°C)")
    grid::popViewport()
    for (i in 1:2) {
      contrast_name <- c("P1 cold","P99 heat")[i]
      d <- lags[contrast==contrast_name][order(lag)]
      stopifnot(identical(d$lag,0:21),all(d$reference_c==ref))
      panel(.685,c(.585,.145)[i],.296,.28,c(0,21),c(.78,1.36))
      draw_curve(d$lag,d$or,d$lower95,d$upper95,colors[[i]])
      reference_line()
      axis_y(c(.8,1,1.2),c("0.8","1.0","1.2"),"Lag-specific OR")
      axis_x(seq(0,21,3),if(i==2) "Lag (days)" else NULL,i==2)
      heading(c("b","c")[i],sprintf("%s, %.2f °C",c("Low temperature","High temperature")[i],unique(d$temperature_c)))
      grid::popViewport()
    }
  }
  for (format in c("pdf","png")) {
    if(format=="pdf") grDevices::pdf(figure_file("Figure2","pdf"),width=8.4,height=4.25,
                   family="Helvetica",encoding="WinAnsi.enc",useDingbats=FALSE)
    else ragg::agg_png(figure_file("Figure2","png"),width=8.4,height=4.25,units="in",res=420,background="white")
    tryCatch(draw_page(), finally=grDevices::dev.off())
  }
}

# Tables 1–2 and S1–S4. Apply the manuscript rounding at display time.
render_tables <- function() {
  load_table <- function(id) fread(file.path(data_dir,paste0("Table",id,".csv")),na.strings=NULL)
  save_table <- function(id,x,headers) {
    setnames(x,headers)
    fwrite(x,file.path(out_dir,paste0("Table",id,".csv")),na="")
  }
  counts <- function(x) formatC(x, format="f", digits=0, big.mark=",")
  count_pct <- function(n,p) sprintf("%s (%.1f)",counts(n),p)
  blank_repeated <- function(x) replace(x,c(FALSE,tail(x,-1)==head(x,-1)),"")
  ci <- function(x,lo,hi,digits) {
    f <- function(v) formatC(v,format="f",digits=digits,big.mark=if(digits==0L) "," else "")
    sprintf("%s (%s, %s)",f(x),f(lo),f(hi))
  }
  a <- load_table("1")
  save_table("1", a[,.(blank_repeated(Characteristic), Category, count_pct(Deaths,Percent))],
             c("Characteristic","Category","Deaths, n (%)"))
  b <- load_table("2")
  save_table("2",b[,.(temperature_category,
             ci(attributable_fraction_percent,attributable_fraction_lower95,attributable_fraction_upper95,2L),
             ci(attributable_deaths,attributable_deaths_lower95,attributable_deaths_upper95,0L))],
             c("Temperature category","Excess fraction, % (95% CI)","Excess number of deaths, n (95% CI)"))
  s <- load_table("S1")
  save_table("S1",s[,.(characteristic,definition,source,reference_year,unit,weighting_population,sprintf("%.2f",cutoff))],
             c("Characteristic","Definition","Source","Reference year","Unit","Population used for weighting","Population-weighted median"))
  s <- load_table("S2")
  panels <- lapply(unique(s$panel),function(label) {
    rbind(data.table(condition=label,code="",deaths=""),
          s[panel==label,.(condition,code=icd10_definition,deaths=count_pct(deaths,percent))])
  })
  save_table("S2",rbindlist(panels),c("Clinical role and condition","ICD-10 code(s)","Deaths, n (%)"))
  s <- load_table("S3")
  s[,Meteorology:=blank_repeated(Meteorology)]
  cols <- c("Mean","SD","Min","P1","Median","P99","Max")
  s[, (cols):=lapply(.SD,function(v)sprintf("%.2f",v)),.SDcols=cols]
  save_table("S3",s,c("Meteorological variable","Macroregion","Mean","SD","Minimum","1st percentile","Median","99th percentile","Maximum"))
  s <- load_table("S4")
  save_table("S4",s[,.(model_specification,ci(low_or,low_lower95,low_upper95,3L),ci(high_or,high_lower95,high_upper95,3L))],
             c("Model specification","Low temperature: P1 vs MMT, OR (95% CI)","High temperature: P99 vs MMT, OR (95% CI)"))
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)>1L) stop("Usage: Rscript code/03_figures_tables.R [output_directory]")
  script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
  root <- dirname(dirname(normalizePath(script)))
  data_dir <- file.path(root,"data")
  out_dir <- if(length(args)) args[1] else paste0(root,"_outputs")
  protected <- c(file.path(root,"code"), data_dir)
  candidate <- normalizePath(out_dir, mustWork=FALSE)
  if(candidate==root || any(candidate==protected) || any(startsWith(candidate,paste0(protected,"/"))))
    stop("Choose an output directory separate from code and data.")
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  out_dir <- normalizePath(out_dir,mustWork=TRUE)
  if(!nzchar(Sys.which("gs"))) stop("Ghostscript (gs) must be available on PATH.")
  render_tables()
  map_vertices <- fread(file.path(data_dir,"Figures1_S3_map_polygons.csv"),colClasses=c(key="character"))
  spatial <- tempfile(fileext=".pdf")
  temporal <- tempfile(fileext=".pdf")
  draw_spatial_panels(fread(file.path(data_dir,"Figure1_municipal.csv"),colClasses=c(Code="character")),spatial)
  draw_temporal_panels(fread(file.path(data_dir,"Figure1_monthly.csv")),temporal)
  compose_figure1(spatial,temporal,figure_file("Figure1","pdf"))
  draw_figure2()
  forest <- fread(file.path(data_dir,"Figures3_S4_forest.csv"))
  draw_figure3(forest,"Individual","Figure3")
  draw_figure3(forest,"Municipality","FigureS4")
  draw_figure_s1(fread(file.path(data_dir,"FigureS1.csv")))
  draw_figure_s2(fread(file.path(data_dir,"FigureS2.csv")))
  draw_figure_s3(list(values=fread(file.path(data_dir,"FigureS3.csv"),colClasses=c(Code="character")),
                      cutoffs=fread(file.path(data_dir,"TableS1.csv"))))
  draw_figure_s5(fread(file.path(data_dir,"FigureS5.csv")))
  message("Eight figures and six tables written to ",out_dir)
}
