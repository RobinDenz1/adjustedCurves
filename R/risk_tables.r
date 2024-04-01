
## calculate unadjusted or weighted number at risk or number of
## censored observations for the entire supplied dataset
get_risk_table.all <- function(times, data, ev_time, event=NULL,
                               type="n_at_risk", weights=NULL) {

  # set weights to all 1 if not specified
  if (is.null(weights)) {
    weights <- rep(1, nrow(data))
  }

  # number at risk
  if (type=="n_at_risk") {
    est <- vapply(times, function(x){sum(weights[data[, ev_time]>=x])},
                  FUN.VALUE=numeric(1))
  # cumulative number of right-censored observations
  } else if (type=="n_cens") {
    est <- vapply(times, FUN=function(x){sum(weights[data[, ev_time] < x &
                                             data[, event]==0])},
                       FUN.VALUE=numeric(1))
  # cumulative number of events
  } else if (type=="n_events") {
    est <- vapply(times, FUN=function(x){sum(weights[data[, ev_time] < x &
                                                     data[, event]==1])},
                  FUN.VALUE=numeric(1))
  }
  out <- data.frame(time=times, est=est)

  return(out)
}

## calculate number at risk or number of censored observations
## stratified by levels of variable
get_risk_table.groups <- function(times, data, variable, ev_time, event=NULL,
                                  type="n_at_risk", weights=NULL) {

  levs <- levels(data[, variable])
  out <- vector(mode="list", length=length(levs))
  for (i in seq_len(length(levs))) {
    dat_temp <- data[data[, variable]==levs[i],]

    if (!is.null(weights)) {
      weights_i <- weights[data[, variable]==levs[i]]
    } else {
      weights_i <- NULL
    }

    out_i <- get_risk_table.all(times=times, data=dat_temp, ev_time=ev_time,
                                event=event, type=type, weights=weights_i)

    out_i$group <- levs[i]
    out[[i]] <- out_i
  }
  out <- dplyr::bind_rows(out)
  out$group <- factor(out$group, levels=levs)

  return(out)
}

## calculate (weighted) risk table from multiply imputed data
get_risk_table.mi <- function(times, data, ev_time, variable, event,
                              type, weights) {

  time <- est <- group <- NULL

  # get fully imputed datasets
  data_all <- mice::complete(data, action="long", include=FALSE)

  # calculate risk table for each one
  out <- vector(mode="list", length=max(data_all$.imp))
  for (i in seq_len(max(data_all$.imp))) {
    data_i <- data_all[data_all$.imp==i, ]

    if (!is.null(weights)) {
      weights_i <- weights[[i]]
    } else {
      weights_i <- NULL
    }

    out[[i]] <- get_risk_table(times=times, data=data_i, ev_time=ev_time,
                               variable=variable, event=event, type=type,
                               weights=weights_i)
  }
  out <- dplyr::bind_rows(out)

  # pool results
  if (is.null(variable)) {
    out <- out %>%
      dplyr::group_by(time) %>%
      dplyr::summarise(est=mean(est),
                       .groups="drop_last")
    out <- out[order(out$time), ]
  } else {
    out <- out %>%
      dplyr::group_by(time, group) %>%
      dplyr::summarise(est=mean(est),
                       .groups="drop_last")
    out <- out[order(out$group, out$time), ]
  }

  return(as.data.frame(out))
}

## estimate general risk table
get_risk_table <- function(times, data, ev_time, variable=NULL, event=NULL,
                           type="n_at_risk", weights=NULL, digits=1) {

  ## with multiple imputation
  if (inherits(data, "mids")) {
    out <- get_risk_table.mi(times=times, data=data, variable=variable,
                             ev_time=ev_time, event=event, type=type,
                             weights=weights)
  ## regular analysis
  } else {
    if (!is.null(variable)) {
      out <- get_risk_table.groups(times=times, data=data, variable=variable,
                                   ev_time=ev_time, event=event, type=type,
                                   weights=weights)
    } else {
      out <- get_risk_table.all(times=times, data=data, ev_time=ev_time,
                                event=event, type=type, weights=weights)
    }
  }

  # remove NA values
  out <- stats::na.omit(out)

  # round values (only relevant if weighted or MI)
  out$est <- round(out$est, digits=digits)

  return(out)
}

## plot unstratified risk table
plot_risk_table.all <- function(plotdata, breaks,
                                gg_theme=ggplot2::theme_classic(),
                                text_size=4.2, text_alpha=1, text_color="black",
                                text_family="sans", text_fontface="plain",
                                color_groups=FALSE, vjust=5,
                                reverse_order=FALSE, custom_colors=NULL) {

  p <- ggplot2::ggplot(data=plotdata,
                       ggplot2::aes(x=.data$time, y=1, label=.data$est)) +
    ggplot2::geom_text(size=text_size, alpha=text_alpha, color=text_color,
                       family=text_family, fontface=text_fontface) +
    gg_theme +
    ggplot2::scale_x_continuous(breaks=breaks) +
    ggplot2::scale_y_continuous(breaks=1) +
    ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                   axis.title.y=ggplot2::element_text(vjust=vjust))

  return(p)
}

## plot stratified risk table
plot_risk_table.groups <- function(plotdata, breaks,
                                   gg_theme=ggplot2::theme_classic(),
                                   text_size=4.2, text_alpha=1,
                                   text_color="black", text_family="sans",
                                   text_fontface="plain", color_groups=TRUE,
                                   reverse_order=TRUE, custom_colors=NULL,
                                   vjust=5) {

  mapping <- ggplot2::aes(x=.data$time,
                          y=.data$group,
                          color=.data$group,
                          label=.data$est)

  if (!color_groups) {
    mapping$colour <- NULL
  }

  if (color_groups) {
    main_geom <- ggplot2::geom_text(size=text_size, alpha=text_alpha,
                                    family=text_family, fontface=text_fontface)
  } else {
    main_geom <- ggplot2::geom_text(size=text_size, alpha=text_alpha,
                                    family=text_family, fontface=text_fontface,
                                    color=text_color)
  }

  p <- ggplot2::ggplot(data=plotdata, mapping) +
    main_geom +
    gg_theme +
    ggplot2::scale_x_continuous(breaks=breaks) +
    ggplot2::theme(legend.position="none")

  if (!is.null(custom_colors)) {
    p <- p + ggplot2::scale_colour_manual(values=custom_colors)
  }

  if (reverse_order) {
    p <- p + ggplot2::scale_y_discrete(limits=rev)
  }

  return(p)
}

## general function to plot all supported types of risk tables
plot_risk_table <- function(data, ev_time, event=NULL, variable=NULL,
                            type, times, xlab="Time", ylab="default",
                            title="default", title_position="middle",
                            title_size=14, weights=NULL, digits=1,
                            text_format=TRUE, additional_layers=list(), ...) {

  plotdata <- get_risk_table(data=data, times=times, ev_time=ev_time,
                             event=event, variable=variable, type=type,
                             weights=weights, digits=digits)
  if (text_format) {
    plotdata$est <- format(plotdata$est, trim=TRUE)
  }

  if (!is.null(variable)) {
    p <- plot_risk_table.groups(plotdata=plotdata, breaks=times, ...)
  } else {
    p <- plot_risk_table.all(plotdata=plotdata, breaks=times, ...)
  }

  defaults <- get_default_labs_rt(type=type, ylab=ylab, title=title,
                                  weights=weights, variable=variable)
  title <- defaults$title
  ylab <- defaults$ylab

  # add axis labels
  p <- p + ggplot2::labs(y=ylab, x=xlab, title=title)

  ## changing title position
  if (title_position=="left") {
    hjust <- 0
  } else if (title_position=="right") {
    hjust <- 1
  } else if (title_position=="middle") {
    hjust <- 0.5
  }
  p <- p + ggplot2::theme(plot.title=ggplot2::element_text(hjust=hjust,
                                                           size=title_size))

  # potentially add more stuff to the plot
  if (length(additional_layers) > 0) {
    for (i in seq_len(length(additional_layers))) {
      p <- p + additional_layers[[i]]
    }
  }

  return(p)
}

## get default values for title and ylab from given values
get_default_labs_rt <- function(type, ylab, title, weights, variable) {

  # get default title / ylab depending on type and whether weights
  # were used
  if (!is.null(weights)) {
    if (type=="n_at_risk") {
      title_default <- "Weighted Number at Risk"
      ylab_default <- "Weighted Number\nat Risk"
    } else if (type=="n_cens") {
      title_default <- "Weighted Cumulative Number Censored"
      ylab_default <- "Weighted Cum.\nNumber Censored"
    } else if (type=="n_events") {
      title_default <- "Weighted Cumulative Number of Events"
      ylab_default <- "Weighted Cum.\nNumber of Events"
    }
  } else {
    if (type=="n_at_risk") {
      title_default <- "Number at Risk"
      ylab_default <- "Number at\nRisk"
    } else if (type=="n_cens") {
      title_default <- "Cumulative Number Censored"
      ylab_default <- "Cum. Number\nCensored"
    } else if (type=="n_events") {
      title_default <- "Cumulative Number of Events"
      ylab_default <- "Cum. Number\nof Events"
    }
  }

  # set default title if not specified
  if (!is.null(variable) && !is.null(title) && title=="default") {
    title <- title_default
  } else if (!is.null(title) && title=="default") {
    title <- NULL
  }

  # set default ylab if not specified
  if (!is.null(variable) && !is.null(ylab) && ylab=="default") {
    ylab <- "Group"
  } else if (!is.null(ylab) && ylab=="default") {
    ylab <- ylab_default
  }

  return(list(title=title, ylab=ylab))
}

## add a risk table to an existing plot of survival curves
add_risk_table <- function(p_surv, ..., height=0.25) {

  requireNamespace("cowplot")

  # x-axis breaks used in survival curve plot
  breaks_p_surv <-
    ggplot2::ggplot_build(p_surv)$layout$panel_params[[1]]$x$breaks

  # create risk table plot
  p_risk_table <- plot_risk_table(times=breaks_p_surv, ...)

  # extract limits from both plots
  limits_rt <-
    ggplot2::ggplot_build(p_risk_table)$layout$panel_params[[1]]$x$limits
  limits_surv <- ggplot2::ggplot_build(p_surv)$layout$panel_params[[1]]$x$limits

  # re-scale plots to have the same limits
  if (limits_rt[2] > limits_surv[2]) {
    p_surv$scales$scales[[1]]$limits <- limits_rt
    p_surv$scales$scales[[1]]$breaks <- breaks_p_surv
  } else if (limits_rt[2] < limits_surv[2]) {
    p_risk_table$scales$scales[[1]]$limits <- limits_surv
    p_risk_table$scales$scales[[1]]$breaks <- breaks_p_surv
  }

  # put plots together
  p <- cowplot::plot_grid(p_surv, p_risk_table, ncol=1, align="v",
                          rel_heights=c(1-height, height), axis="lr")
  return(p)
}
