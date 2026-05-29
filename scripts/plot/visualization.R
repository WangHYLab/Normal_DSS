library(ggplot2)
library(dplyr)
library(forcats)
library(ggridges)
library(viridis)

my_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.text = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold",size = 14),
      legend.key = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )
}

plot_density_ridges_group <- function(
    df,
    group_var,
    value_var,
    title = "Density distribution across groups",
    font_size = 13,
    ref_group = "All drugs",
    ref_fill = "grey80",
    other_fill = "orange",
    show_legend = FALSE
) {
  group_var <- rlang::ensym(group_var)
  value_var <- rlang::ensym(value_var)

  df_plot <- df %>%
    dplyr::filter(!is.na(!!value_var)) %>%
    dplyr::group_by(!!group_var) %>%
    dplyr::mutate(median_value = median(!!value_var, na.rm = TRUE)) %>%
    dplyr::ungroup()

  df_plot[[rlang::as_string(group_var)]] <- forcats::fct_reorder(
    df_plot[[rlang::as_string(group_var)]],
    df_plot$median_value
  )

  df_plot$fill_group <- ifelse(
    as.character(df_plot[[rlang::as_string(group_var)]]) == ref_group,
    ref_group,
    "other"
  )

  fill_values <- c(ref_fill, other_fill)
  names(fill_values) <- c(ref_group, "other")

  p <- ggplot(
    df_plot,
    aes(x = !!value_var, y = !!group_var, fill = fill_group)
  ) +
    ggridges::geom_density_ridges(
      scale = 1.2,
      rel_min_height = 0.01,
      color = "white",
      size = 0.3,
      alpha = 0.95
    ) +
    scale_fill_manual(
      values = fill_values,
      guide = if (show_legend) "legend" else "none"
    ) +
    ggridges::theme_ridges(font_size = font_size, grid = TRUE) +
    theme(
      axis.title.y = element_blank()
    ) +
    labs(
      title = title,
      x = rlang::as_string(value_var),
      y = NULL
    )

  return(p + my_theme())
}

get_upset <- function(listInput,
                     height = unit(3, "cm"),
                     width  = unit(5, "cm"),
                     comb_col = NULL) {

  require(ComplexHeatmap)
  require(grid)

  comb = make_comb_mat(listInput)
  comb_sets = lapply(comb_name(comb), function(nm) extract_comb(comb, nm))

  if (is.null(comb_col)) {
    c_col = viridis::cividis(max(comb_degree(comb)))
    comb_col = c_col[comb_degree(comb)]
  }
  comb_ordering = order(comb_size(comb), decreasing = TRUE)

  UpSet(
    comb,
    comb_order = comb_ordering,
    pt_size = unit(4, "mm"),
    lwd = 3,
    height = height,
    width  = width,
    comb_col = comb_col,
    top_annotation = upset_top_annotation(comb, add_numbers = TRUE),
    right_annotation = upset_right_annotation(
      comb,
      add_numbers = TRUE,
      width = unit(1.5, "cm")
    )
  )
}

plot_stacked_bar <- function(df,
                             x_col,
                             y_col,
                             fill_col = NULL,
                             agg_func = dplyr::n_distinct,
                             plot_title = "",
                             xlab = NULL,
                             ylab = NULL,
                             fill_legend = NULL,
                             palette = "Set2") {
  library(dplyr)
  library(ggplot2)
  library(rlang)
  library(forcats)
  library(RColorBrewer)

  x_sym <- sym(x_col)
  y_sym <- sym(y_col)

  if (!is.null(fill_col)) {
    fill_sym <- sym(fill_col)

    summary_df <- df %>%
      group_by(!!x_sym, !!fill_sym) %>%
      summarise(value = agg_func(!!y_sym), .groups = "drop")

    total_height <- summary_df %>%
      group_by(!!x_sym) %>%
      summarise(total = sum(value), .groups = "drop")

    summary_df <- summary_df %>%
      left_join(total_height, by = rlang::as_string(x_sym)) %>%
      mutate(!!x_sym := fct_reorder(!!x_sym, total, .desc = TRUE))

    n_fill <- n_distinct(summary_df[[fill_col]])

    p <- ggplot(summary_df, aes(x = !!x_sym, y = value, fill = !!fill_sym)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(
        title = plot_title,
        x = xlab %||% x_col,
        y = ylab %||% paste0("Aggregated ", y_col),
        fill = fill_legend %||% fill_col
      )

    if (n_fill <= 8) {
      p <- p + scale_fill_brewer(palette = palette)
    } else {
      p <- p + scale_fill_manual(values = scales::hue_pal()(n_fill))
    }

    p <- p + geom_text(
      data = total_height,
      aes_string(x = x_col, y = "total", label = "total"),
      inherit.aes = FALSE,
      vjust = -0.3,
      size = 2.5
    )

  } else {
    summary_df <- df %>%
      group_by(!!x_sym) %>%
      summarise(value = agg_func(!!y_sym), .groups = "drop") %>%
      mutate(!!x_sym := fct_reorder(!!x_sym, value, .desc = TRUE))

    p <- ggplot(summary_df, aes(x = !!x_sym, y = value)) +
      geom_bar(stat = "identity", fill = "#69b3a2") +
      geom_text(aes(label = value), vjust = -0.3) +
      theme_minimal(base_size = 14) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
      labs(
        title = plot_title,
        x = xlab %||% x_col,
        y = ylab %||% paste0("Aggregated ", y_col)
      )
  }

  return(p+my_theme())
}

plot_heatmap_with_annotation <- function(
  mat,
  title_base="DSE",
  anno_df = NULL,
  anno_cols = NULL,
  row_anno_df = NULL,
  row_anno_cols = NULL,
  zscore = FALSE,
  zscore_by = c("column", "row"),
  mark_cols = NULL,
  mark_labels = NULL,
  mark_side = "bottom",
  col_fun = NULL,
  ...
) {
  library(ComplexHeatmap)
  library(circlize)

  zscore_by <- match.arg(zscore_by)

  mat <- as.matrix(mat)

  if (zscore) {
    if (zscore_by == "column") {
      mat_plot <- scale(mat)
    } else {
      mat_plot <- t(scale(t(mat)))
    }

    if (is.null(col_fun)) {
      col_fun <- colorRamp2(
        c(-2, 0, 2),
        c("lightblue", "orange", "red")
      )
    }

    legend_title <- paste0(title_base, "\n(Z-score)")

  } else {
    mat_plot <- mat

    if (is.null(col_fun)) {
      vmin <- min(mat_plot, na.rm = TRUE)
      vmax <- max(mat_plot, na.rm = TRUE)
      vmid <- (vmin + vmax) / 2

      col_fun <- colorRamp2(
        c(vmin, vmid, vmax),
        c("lightblue", "orange", "red")
      )
    }

    legend_title <- title_base
  }

  top_anno <- NULL
  if (!is.null(anno_df)) {
    anno_df <- as.data.frame(anno_df)

    if (is.null(anno_cols)) {
      anno_cols <- colnames(anno_df)
    }

    anno_df <- anno_df[, anno_cols, drop = FALSE]
    anno_df <- anno_df[intersect(rownames(anno_df), colnames(mat_plot)), , drop = FALSE]
    anno_df <- anno_df[colnames(mat_plot), , drop = FALSE]

    anno_col_list <- lapply(anno_df, function(x) {
      if (is.numeric(x)) {
        colorRamp2(
          c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)),
          c("white", "darkblue")
        )
      } else {
        vals <- unique(na.omit(x))
        setNames(
          colorRampPalette(
            RColorBrewer::brewer.pal(12, "Paired")
          )(length(vals)),
          vals
        )
      }
    })

    top_anno <- HeatmapAnnotation(
      df = anno_df,
      col = anno_col_list,
      annotation_name_side = "left"
    )
  }

  right_anno <- NULL

  if (!is.null(row_anno_df)) {
    row_anno_df <- as.data.frame(row_anno_df)

    if (is.null(rownames(row_anno_df))) {
      stop("row_anno_df must have row names matching mat row names")
    }

    if (is.null(row_anno_cols)) {
      row_anno_cols <- colnames(row_anno_df)
    }

    row_anno_df <- row_anno_df[, row_anno_cols, drop = FALSE]

    common_rows <- intersect(rownames(mat_plot), rownames(row_anno_df))
    row_anno_df <- row_anno_df[common_rows, , drop = FALSE]
    mat_plot <- mat_plot[common_rows, , drop = FALSE]

    anno_col_list <- lapply(row_anno_df, function(x) {
      if (is.numeric(x)) {
        colorRamp2(
          c(min(x, na.rm = TRUE), max(x, na.rm = TRUE)),
          c("white", "darkblue")
        )
      } else {
        vals <- unique(na.omit(x))
        setNames(
          colorRampPalette(
            brewer.pal(min(12, length(vals)), "Set3")
          )(length(vals)),
          vals
        )
      }
    })

    right_anno <- rowAnnotation(
      df = row_anno_df,
      col = anno_col_list,
      annotation_name_side = "top"
    )
  }

  mark_anno <- NULL
  if (!is.null(mark_cols)) {

    mark_cols <- intersect(mark_cols, colnames(mat_plot))

    if (length(mark_cols) > 0) {
      at_idx <- match(mark_cols, colnames(mat_plot))

      if (is.null(mark_labels)) {
        mark_labels <- mark_cols
      }

      mark_anno <- HeatmapAnnotation(
        mark = anno_mark(
          at = at_idx,
          labels = mark_labels,
          side = mark_side
        ),
        which = "column"
      )
    }
  }

  Heatmap(
    mat_plot,
    name = legend_title,
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    na_col = "grey90",
    top_annotation = top_anno,
    bottom_annotation = mark_anno,
    right_annotation = right_anno,
    ...
  )

}

draw_nulldist_density_plots <- function(NullDist_list, title_text = "Non-immune null distribution") {
  NullDist <- NullDist_list$NullDist
  threshold_r <- NullDist_list$threshold_r
  threshold_s <- NullDist_list$threshold_s

  df <- data.frame(NullDist = NullDist)

  max_dens <- max(density(df$NullDist)$y)

  ggplot(df, aes(x = NullDist)) +
    geom_density(fill = "lightblue", alpha = 0.5, color = "blue", size = 0.8) +

    geom_vline(xintercept = threshold_r, linetype = "dashed", color = "red", size = 0.8) +
    geom_vline(xintercept = threshold_s, linetype = "dashed", color = "red", size = 0.8) +

    annotate("text",
             x = threshold_s - 0.2,
             y = max_dens * 0.7,
             label = paste0("Threshold_s = ", round(threshold_s, 3)),
             hjust = 1, size = 5) +
    annotate("text",
             x = threshold_r + 0.2,
             y = max_dens * 0.7,
             label = paste0("Threshold_r = ", round(threshold_r, 3)),
             hjust = 0, size = 5) +

    labs(title = title_text,
         x = "NES value",
         y = "Density") +

    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14, color = "black"),
      axis.text = element_text(size = 12, color = "black")
    )
}

plot_drug_rank_consistency <- function(df,
                                       drug_col = "drug_id_name",
                                       sample_col = "sample",
                                       rank_col = "DSE_rank",
                                       drugs_to_plot = NULL,
                                       title = "Rank consistency across samples",
                                       smooth = FALSE,
                                       base_size = 14) {
  library(dplyr)
  library(ggplot2)
  library(scales)

  df <- df %>%
    dplyr::select(all_of(c(drug_col, sample_col, rank_col,'drug_name'))) %>%
    dplyr::filter(!is.na(.data[[rank_col]]))

  if (is.null(drugs_to_plot)) {
    drugs_to_plot <- unique(df[[drug_col]])
  }

  df <- df %>% filter(.data[[drug_col]] %in% drugs_to_plot)

  max_rank_all <- max(df[[rank_col]], na.rm = TRUE)
  df <- df %>% mutate(rank_percentile = .data[[rank_col]] / max_rank_all)

  df_cdf <- df %>%
    group_by(.data[[drug_col]]) %>%
    arrange(rank_percentile, .by_group = TRUE) %>%
    mutate(cum_prop = seq_along(rank_percentile) / n()) %>%
    ungroup()

  p <- ggplot(df_cdf, aes(x = rank_percentile, y = cum_prop, color = .data[['drug_name']])) +
    {
      if (smooth) geom_smooth(se = FALSE, size = 1.3, span = 0.5) else geom_line(size = 1.2)
    } +
    geom_point(size = 1.2, alpha = 0.6) +
    scale_color_brewer(palette = "Set1") +
    theme_bw(base_size = base_size) +
    labs(
      title = ifelse(length(drugs_to_plot) == 1,
                     paste0("Rank consistency of ", drugs_to_plot),
                     title),
      x = "Rank percentile",
      y = "Cumulative fraction of samples"
    )

  return(p+my_theme()+
           theme(legend.title = element_blank()))
}
