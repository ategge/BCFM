#' @title A heatmap of group assignments, Z using ggplot2
#' @description A heatmap of group assignments, Z using ggplot2. It first sorts by the largest group with the most assigned.
#'
#' @param Gibbs Gibbs sample from \code{BCFM} function
#' @param true.val Table of true group assignments, if applicable
#' @param burnin Number of burn-in period. If not specified, it uses the first tenths sample as burn-in.
#' @param main.title Title of the plot. Default is "Cluster Assignment Heatmap"
#' @param x.label X-axis label. Default is "Cluster"
#' @param y.label Y-axis label. Default is "Subject-Time"
#'
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot aes geom_tile geom_line scale_fill_gradient theme ggtitle labs theme element_blank
#'
#' @return A ggplot object
#' @export ggplot_Zit.heatmap
#' @export

ggplot_Zit.heatmap <- function(Gibbs,
                                true.val = NA,
                                burnin = NA,
                                main.title = "Cluster Assignment Heatmap",
                                x.label = "Cluster",
                                y.label = "Subject-Time"){

  # true.val should be the original group, ordered by the first subject's observations by all time
  if(is.null(Gibbs$Zit)){Gibbs$Zit <- Gibbs$Z}
  n.iter <- dim(Gibbs$Zit)[1]
  n <- dim(Gibbs$Zit)[2]
  times <- dim(Gibbs$Zit)[3]
  G <- dim(Gibbs$Omega)[2]
  if(is.na(burnin)){burnin <- round(n.iter / 10)}
  turns <- seq(burnin + 1, n.iter)

  Z.long <- matrix(NA, n*times, G)
  Z.subject <- paste("S", rep(1:n, each = times), ".", rep(1:times, n), sep = "")
  Z.group <- paste("G", 1:G, sep = "")
  colnames(Z.long) <- Z.group
  rownames(Z.long) <- Z.subject

  if(is.na(true.val[1])){
    for(i in 1:n){
      for(tt in 1:times){
        current <- times*(i-1) + tt
        z.factor <- factor(Gibbs$Zit[turns, i, tt], levels = 1:G)
        Z.long[current,] <- table(z.factor)/length(turns)
      }
    }
    Z.long <- as.data.frame(Z.long)
    Z.long.data <- rbind(rep(NA, G))
    Z.long.subject <- c()
    Z.long$subject <- Z.subject

    Z.order <- apply(Z.long[,1:G], 1, which.max)
    for(k1 in 1:G){
      Z.mat.current <- as.matrix(Z.long[Z.order == k1, 1:G])
      Z.subject.current <- Z.subject[Z.order == k1]
      for(k2 in G:1){
        if(k2 == k1){next}
        Z.mat.current <- Z.mat.current[order(Z.mat.current[,k2], decreasing = TRUE),]
        Z.subject.current <- Z.subject.current[order(Z.mat.current[,k2], decreasing = TRUE)]
      }
      Z.mat.current <- Z.mat.current[order(Z.mat.current[,k1], decreasing = TRUE),]
      Z.subject.current <- Z.subject.current[order(Z.mat.current[,k1], decreasing = TRUE)]
      Z.long.data <- rbind(Z.long.data, Z.mat.current)
      Z.long.subject <- c(Z.long.subject, Z.subject.current)
    }
    Z.long.data <- Z.long.data[-1,]
    Z.long[,1:G] <- Z.long.data
    Z.long$subject <- Z.long.subject

    Z.long <- tidyr::gather(Z.long, group, assigned, -c(subject))
    Z.long$subject <- factor(Z.long$subject, levels = (Z.long$subject)[(n*times):1])

    p <- ggplot(Z.long, aes(x = group, y = subject, fill = assigned)) +
      scale_fill_gradient(low = "gray80", high = "black", name = "Probability") +
      geom_tile() +
      theme_minimal() +
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid = element_blank()) +
      labs(title = main.title,
           x = x.label,
           y = y.label)
  }

  if(!is.na(true.val)[1]){
    for(i in 1:n){
      for(tt in 1:times){
        current <- times*(i-1) + tt
        z.factor <- factor(Gibbs$Zit[turns, i, tt], levels = 1:G)
        Z.long[current,] <- table(z.factor)/length(turns)
      }
    }
    Z.long <- Z.long[order(true.val),]
    Z.long <- as.data.frame(Z.long)
    Z.long$subject <- Z.subject
    group.original.long <- c(true.val)
    start.point <- c(1, cumsum(table(group.original.long)) + 1)[-(G+1)]
    end.point <- cumsum(table(group.original.long))
    for(j in 1:G){
      correct.probs <- Z.long[seq(start.point[j], end.point[j]), j]
      Z.long[seq(start.point[j], end.point[j]),] <- Z.long[seq(start.point[j], end.point[j])[order(correct.probs, decreasing = TRUE)],]
    }
    Z.long <- tidyr::gather(Z.long, group, assigned, -c(subject))
    Z.long$subject <- factor(Z.long$subject, levels = (Z.long$subject)[(n*times):1])
    Z.line <- data.frame(x = c(G, 0) + 0.5, y = rep(n*times - end.point, each = 2) + 0.5)

    p <- ggplot(Z.long) +
      scale_fill_gradient(low = "gray80", high = "black", name = "Probability") +
      geom_tile(aes(x = group, y = subject, fill = assigned)) +
      theme_minimal() +
      theme(axis.ticks.y = element_blank(),
            axis.text.y = element_blank(),
            panel.grid = element_blank()) +
      labs(title = main.title,
           x = x.label,
           y = y.label) +
      geom_line(data = Z.line, aes(x = x, y = y, group = y), color = "blue", linewidth = 1)
  }

  return(p)
}
