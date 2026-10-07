# Shared spatial plotting geometry for scripts 02, 04, 05, and 06.
#
# Biological spot diameters and arrow shafts are expressed in low-resolution
# image pixels.  Physical units are reserved for purely graphical styling such
# as outline width and arrow-head size.

spatial_plot_core_api_version <- "1.0.0"

sp_read_context <- function(spaceranger_dir, load_image = TRUE) {
  spatial_dir <- file.path(spaceranger_dir, "spatial")
  scale_path <- file.path(spatial_dir, "scalefactors_json.json")
  if (!file.exists(scale_path)) {
    stop("Missing Space Ranger scalefactors JSON: ", scale_path, call. = FALSE)
  }

  scale_factors <- jsonlite::fromJSON(scale_path)
  required <- c("spot_diameter_fullres", "tissue_lowres_scalef")
  missing_values <- required[vapply(required, function(x) {
    is.null(scale_factors[[x]]) || !is.finite(as.numeric(scale_factors[[x]]))
  }, logical(1))]
  if (length(missing_values) > 0) {
    stop(
      "Space Ranger scalefactors JSON is missing finite value(s): ",
      paste(missing_values, collapse = ", "),
      call. = FALSE
    )
  }

  lowres_scale <- as.numeric(scale_factors$tissue_lowres_scalef)
  spot_diameter_lowres <- as.numeric(scale_factors$spot_diameter_fullres) * lowres_scale
  if (!is.finite(spot_diameter_lowres) || spot_diameter_lowres <= 0) {
    stop("Calculated low-resolution spot diameter must be positive.", call. = FALSE)
  }

  image_path <- file.path(spatial_dir, "tissue_lowres_image.png")
  has_image <- isTRUE(load_image) && file.exists(image_path)
  image <- NULL
  image_grob <- NULL
  image_width <- NA_integer_
  image_height <- NA_integer_
  if (has_image) {
    image <- png::readPNG(image_path)
    image_height <- dim(image)[[1]]
    image_width <- dim(image)[[2]]
    image_grob <- grid::rasterGrob(
      image,
      interpolate = FALSE,
      width = grid::unit(1, "npc"),
      height = grid::unit(1, "npc")
    )
  }

  list(
    spatial_dir = spatial_dir,
    scale_path = scale_path,
    image_path = image_path,
    has_image = has_image,
    image = image,
    image_grob = image_grob,
    image_width = image_width,
    image_height = image_height,
    tissue_lowres_scalef = lowres_scale,
    spot_diameter_lowres = spot_diameter_lowres,
    spot_radius_lowres = spot_diameter_lowres / 2
  )
}

sp_make_spot_polygons <- function(
    data,
    context,
    x_col = "X",
    y_col = "Y",
    segments = 24L,
    flip_y = TRUE) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (nrow(data) == 0) {
    data$.spot_group <- integer()
    data$.spot_x <- numeric()
    data$.spot_y <- numeric()
    return(data)
  }
  missing_cols <- setdiff(c(x_col, y_col), colnames(data))
  if (length(missing_cols) > 0) {
    stop("Spot data is missing coordinate column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  segments <- as.integer(segments)
  if (is.na(segments) || segments < 12L) {
    stop("Spot circle segments must be an integer >= 12.", call. = FALSE)
  }

  x <- as.numeric(data[[x_col]])
  y <- as.numeric(data[[y_col]])
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    stop("Spot coordinates must all be finite.", call. = FALSE)
  }
  if (isTRUE(flip_y)) y <- -y

  angles <- seq(0, 2 * pi, length.out = segments + 1L)[-1L]
  row_index <- rep(seq_len(nrow(data)), each = segments)
  angle_index <- rep(seq_len(segments), times = nrow(data))
  out <- data[row_index, , drop = FALSE]
  rownames(out) <- NULL
  out$.spot_group <- row_index
  out$.spot_x <- x[row_index] + context$spot_radius_lowres * cos(angles[angle_index])
  out$.spot_y <- y[row_index] + context$spot_radius_lowres * sin(angles[angle_index])
  out
}

sp_spot_layer <- function(
    spot_polygons,
    fill_col = NULL,
    colour_col = NULL,
    fill = NA,
    colour = "grey20",
    alpha = 0.82,
    linewidth = 0.1,
    ...) {
  mapping <- ggplot2::aes(x = .spot_x, y = .spot_y, group = .spot_group)
  if (!is.null(fill_col)) mapping$fill <- as.name(fill_col)
  if (!is.null(colour_col)) mapping$colour <- as.name(colour_col)

  args <- list(
    data = spot_polygons,
    mapping = mapping,
    inherit.aes = FALSE,
    alpha = alpha,
    linewidth = linewidth,
    na.rm = TRUE
  )
  if (is.null(fill_col)) args$fill <- fill
  if (is.null(colour_col)) args$colour <- colour
  args <- c(args, list(...))
  do.call(ggplot2::geom_polygon, args)
}

sp_spatial_canvas <- function(context, show_image = TRUE, title = NULL) {
  p <- ggplot2::ggplot()
  if (isTRUE(show_image) && isTRUE(context$has_image)) {
    p <- p + ggplot2::annotation_custom(
      grob = context$image_grob,
      xmin = 0,
      xmax = context$image_width,
      ymin = -context$image_height,
      ymax = 0
    )
  }

  if (is.finite(context$image_width) && is.finite(context$image_height)) {
    p <- p + ggplot2::coord_fixed(
      ratio = 1,
      xlim = c(0, context$image_width),
      ylim = c(-context$image_height, 0),
      expand = FALSE,
      clip = "on"
    )
  } else {
    p <- p + ggplot2::coord_fixed(ratio = 1, expand = FALSE, clip = "on")
  }

  p <- p +
    ggplot2::theme_void(base_size = 11, base_family = "Arial") +
    ggplot2::theme(
      legend.position = "right",
      legend.box = "vertical",
      legend.justification = "left",
      plot.title = ggplot2::element_text(hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.margin = ggplot2::margin(4, 4, 4, 4)
    )
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p
}

sp_boundary_base <- function(
    spot_polygons,
    context,
    boundary_cols,
    title = NULL,
    show_image = TRUE,
    legend_name = "Location",
    alpha = 0.82,
    border_colour = "grey20",
    border_linewidth = 0.1) {
  sp_spatial_canvas(context, show_image = show_image, title = title) +
    sp_spot_layer(
      spot_polygons,
      fill_col = "Location",
      colour = border_colour,
      alpha = alpha,
      linewidth = border_linewidth,
      key_glyph = "point"
    ) +
    ggplot2::scale_fill_manual(values = boundary_cols, drop = FALSE, name = legend_name) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        override.aes = list(
          shape = 21,
          size = 2,
          colour = border_colour,
          alpha = 1
        )
      )
    )
}

sp_prepare_arrows <- function(
    data,
    arrow_length_scale,
    x_col = "X",
    y_col = "Y",
    vx_col = "vx.u",
    vy_col = "vy.u",
    flip_y = TRUE) {
  out <- as.data.frame(data, stringsAsFactors = FALSE)
  if (nrow(out) == 0) {
    out$X_end <- numeric()
    out$Y_end <- numeric()
    return(out)
  }
  required <- c(x_col, y_col, vx_col, vy_col)
  missing_cols <- setdiff(required, colnames(out))
  if (length(missing_cols) > 0) {
    stop("Arrow data is missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  arrow_length_scale <- as.numeric(arrow_length_scale)
  if (length(arrow_length_scale) != 1L || !is.finite(arrow_length_scale) || arrow_length_scale <= 0) {
    stop("arrow_length_scale must be one finite positive number.", call. = FALSE)
  }

  start_x <- as.numeric(out[[x_col]])
  start_y <- as.numeric(out[[y_col]])
  end_x <- start_x + as.numeric(out[[vx_col]]) * arrow_length_scale
  end_y <- start_y + as.numeric(out[[vy_col]]) * arrow_length_scale
  out$X <- start_x
  out$Y <- if (isTRUE(flip_y)) -start_y else start_y
  out$X_end <- end_x
  out$Y_end <- if (isTRUE(flip_y)) -end_y else end_y
  out
}

sp_fixed_layout <- function(plot, map_width = 6.45, legend_width = 1.55) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for fixed spatial plot layout.", call. = FALSE)
  }
  (plot + patchwork::guide_area()) +
    patchwork::plot_layout(widths = c(map_width, legend_width), guides = "collect") &
    ggplot2::theme(legend.position = "right")
}

sp_save_plot <- function(
    plot,
    filename,
    width = 8,
    height = 7,
    fixed_layout = TRUE,
    map_width = 6.45,
    legend_width = 1.55,
    ...) {
  output_plot <- if (isTRUE(fixed_layout)) {
    sp_fixed_layout(plot, map_width = map_width, legend_width = legend_width)
  } else plot
  extension <- tolower(tools::file_ext(filename))
  if (identical(extension, "pdf") && isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(filename, width = width, height = height, family = "Arial", onefile = TRUE)
    device_open <- TRUE
    tryCatch(
      print(output_plot),
      finally = {
        if (isTRUE(device_open)) grDevices::dev.off()
      }
    )
  } else if (identical(extension, "svg")) {
    if (!requireNamespace("svglite", quietly = TRUE)) {
      stop("Package 'svglite' is required for SVG export.", call. = FALSE)
    }
    svglite::svglite(filename, width = width, height = height)
    device_open <- TRUE
    tryCatch(
      print(output_plot),
      finally = {
        if (isTRUE(device_open)) grDevices::dev.off()
      }
    )
  } else {
    ggplot2::ggsave(filename, output_plot, width = width, height = height, ...)
  }
  invisible(filename)
}
