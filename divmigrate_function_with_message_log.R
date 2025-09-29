#modified Divmigrate function to show progress update:

divMigrate <- function (infile = NULL, outfile = NULL, boots = 0, stat = "all", 
          filter_threshold = 0, plot_network = FALSE, plot_col = "darkblue", 
          para = FALSE) 
{
  dat <- rgp(infile)
  npops <- length(dat$genos)
  nloci <- length(dat$af)
  dat$af <- lapply(dat$af, function(x) {
    cs <- colSums(x)
    x[, cs == 0] <- NA
    return(x)
  })
  
  message("Computing pairwise Ht and Hs values...")

  
  if (!is.null(outfile)) {
    dir.create(path = paste(getwd(), "/", outfile, "-[divMigrate]", 
                            "/", sep = ""))
    of <- paste(getwd(), "/", outfile, "-[divMigrate]", 
                "/", sep = "")
  }
  pw <- combn(npops, 2)
  hths <- lapply(dat$af, pwHt, pw = pw - 1)
  ht <- lapply(hths, "[[", "ht")
  hs <- lapply(hths, "[[", "hs")
  if (stat == "d" || stat == "all" || stat == "Nm") {
    message("Calculating Nei’s D and relative migration...")
    d <- function(ht, hs) {
      return(((ht - hs)/(1 - hs)) * 2)
    }
  }
  if (stat == "gst" || stat == "all" || stat == "Nm") {
    message("Calculating Gst and relative migration...")
    g <- function(ht, hs) {
      ot <- (ht - hs)/ht
      diag(ot) <- 0
      return(ot)
    }
  }
  if (stat == "Nm" || stat == "all") {
    Nm <- function(g, d, n) {
      t1 <- (1 - g)/g
      t2 <- ((n - 1)/n)^2
      t3 <- ((1 - d)/(1 - ((n - 1)/n) * d))
      return(0.25 * t1 * t2 * t3)
    }
  }
  if (stat == "d" || stat == "all" || stat == "Nm") {
    dloc <- mapply(d, ht = ht, hs = hs, SIMPLIFY = "array")
    dloc[is.nan(dloc)] <- 1
    hrmD <- apply(dloc, c(1, 2), function(x) {
      mn <- mean(x, na.rm = TRUE)
      vr <- var(x, na.rm = TRUE)
      return(1/((1/mn) + vr * (1/mn)^3))
    })
    dMig <- (1 - hrmD)/hrmD
    dMig[is.infinite(dMig)] <- NA
    dRel <- dMig/max(dMig, na.rm = TRUE)
    dRel[is.nan(dRel)] <- NA
  }
  if (stat == "gst" || stat == "all" || stat == "Nm") {
    g <- function(ht, hs) {
      ot <- (ht - hs)/ht
      diag(ot) <- 0
      return(ot)
    }
    hsAr <- array(unlist(hs), dim = c(npops, npops, nloci))
    mnHs <- apply(hsAr, c(1, 2), mean, na.rm = TRUE)
    htAr <- array(unlist(ht), dim = c(npops, npops, nloci))
    mnHt <- apply(htAr, c(1, 2), mean, na.rm = TRUE)
    hrmGst <- g(mnHt, mnHs)
    gMig <- ((1/hrmGst) - 1)/4
    gMig[is.infinite(gMig)] <- NA
    gRel <- gMig/max(gMig, na.rm = TRUE)
  }
  if (stat == "all" || stat == "Nm") {
    nm <- Nm(hrmGst, hrmD, 2)
    diag(nm) <- NA
    nmRel <- nm/max(nm, na.rm = TRUE)
  }
  if (plot_network || boots != 0L) {
    message("Preparing data for plotting or bootstrapping...")
    if (stat == "d" || stat == "all" || stat == "Nm") {
      dRelPlt <- dRel
      dRelPlt[dRelPlt < filter_threshold] <- 0
    }
    if (stat == "gst" || stat == "all" || stat == "Nm") {
      gRelPlt <- gRel
      gRelPlt[gRelPlt < filter_threshold] <- 0
    }
    if (stat == "all" || stat == "Nm") {
      nmRelPlt <- nmRel
      nmRelPlt[nmRelPlt < filter_threshold] <- 0
    }
  }
  if (plot_network) {
    message("Plotting relative migration network...")
    if (stat == "d") {
      qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                 "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; D)", sep = ""))
      if (!is.null(outfile)) {
        pdf(paste(of, "Relative_migration.pdf", sep = ""), 
            paper = "a4r")
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; D)", sep = ""))
      }
    }
    if (stat == "gst") {
      qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                 "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Gst)", sep = ""))
      if (!is.null(outfile)) {
        pdf(paste(of, "Relative_migration.pdf", sep = ""), 
            paper = "a4r")
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Gst)", sep = ""))
      }
    }
    if (stat == "Nm") {
      qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                  "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Nm)", sep = ""))
      if (!is.null(outfile)) {
        pdf(paste(of, "Relative_migration.pdf", sep = ""), 
            paper = "a4r")
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Nm)", sep = ""))
      }
    }
    if (stat == "all") {
      qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                 "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; D)", sep = ""))
      qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                 "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Gst)", sep = ""))
      qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                  "[", 1), legend = TRUE, posCol = plot_col, edge.labels = TRUE, 
                     mar = c(2, 2, 5, 5), curve = 2.5)
      title(paste("\n Relative migration network \n (Filter threshold = ", 
                  filter_threshold, "; Nm)", sep = ""))
      if (!is.null(outfile)) {
        pdf(paste(of, "Relative_migration.pdf", sep = ""), 
            paper = "a4r")
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; D)", sep = ""))
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Gst)", sep = ""))
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       edge.labels = TRUE, mar = c(2, 2, 5, 5), curve = 2.5)
        title(paste("\n Relative migration network \n (Filter threshold = ", 
                    filter_threshold, "; Nm)", sep = ""))
      }
    }
    if (boots == 0L && !is.null(outfile)) {
      dev.off()
    }
  }
  if (boots != 0L) {
    message(paste("Running bootstraps (", boots, " iterations)...", sep = ""))
    ps <- sapply(dat$indnms, length)
    idx <- lapply(1:boots, function(i) {
      lapply(ps, function(x) {
        return(sample(x, size = x, replace = TRUE))
      })
    })
    if (para) {
      message("Running bootstraps in parallel...")
      cl <- parallel::makeCluster(detectCores())
      parallel::clusterExport(cl, c("bsFun", "dat", "pw", 
                                    "stat"), envir = environment())
      bsStat <- parallel::parLapply(cl, idx, function(x) {
        return(bsFun(genos = dat$genos, idx = x, af = dat$af, 
                     pw = pw, stat = stat))
      })
      parallel::stopCluster(cl)
    }
    else {
      message("Running bootstraps sequentially...")
      bsStat <- vector("list", length = boots)
      for (i in seq_len(boots)) {
        message(paste("Running bootstrap", i, "of", boots))
        bsStat[[i]] <- bsFun(genos = dat$genos, idx = idx[[i]], af = dat$af, pw = pw, stat = stat)
      }
    }
    message("Bootstrapping complete. Computing significance...")
    if (stat == "d" || stat == "all") {
      bsD <- sapply(bsStat, "[[", "dRel", simplify = "array")
    }
    if (stat == "gst" || stat == "all") {
      bsG <- sapply(bsStat, "[[", "gRel", simplify = "array")
    }
    if (stat == "Nm" || stat == "all") {
      bsNm <- sapply(bsStat, "[[", "nmRel", simplify = "array")
    }
    sigDiff <- function(x, y) {
      if (x[1] < y[1] && x[2] < y[1]) {
        return(TRUE)
      }
      else {
        return(FALSE)
      }
    }
    if (stat == "d" || stat == "all") {
      sigMatD <- matrix(NA, nrow = ncol(dRel), ncol(dRel))
      for (i in 1:ncol(pw)) {
        p1 <- quantile(bsD[pw[1, i], pw[2, i], ], prob = c(0.025, 
                                                           0.975))
        p2 <- quantile(bsD[pw[2, i], pw[1, i], ], prob = c(0.025, 
                                                           0.975))
        sigMatD[pw[2, i], pw[1, i]] <- sigDiff(p1, p2)
        sigMatD[pw[1, i], pw[2, i]] <- sigDiff(p2, p1)
      }
      dRelPlt[!sigMatD] <- 0
    }
    if (stat == "gst" || stat == "all") {
      sigMatG <- matrix(NA, nrow = ncol(gRel), ncol(gRel))
      for (i in 1:ncol(pw)) {
        p1 <- quantile(bsG[pw[1, i], pw[2, i], ], prob = c(0.025, 
                                                           0.975))
        p2 <- quantile(bsG[pw[2, i], pw[1, i], ], prob = c(0.025, 
                                                           0.975))
        sigMatG[pw[2, i], pw[1, i]] <- sigDiff(p1, p2)
        sigMatG[pw[1, i], pw[2, i]] <- sigDiff(p2, p1)
      }
      gRelPlt[!sigMatG] <- 0
    }
    if (stat == "Nm" || stat == "all") {
      sigMatNm <- matrix(NA, nrow = ncol(nmRel), ncol(nmRel))
      for (i in 1:ncol(pw)) {
        p1 <- quantile(bsNm[pw[1, i], pw[2, i], ], prob = c(0.025, 
                                                            0.975))
        p2 <- quantile(bsNm[pw[2, i], pw[1, i], ], prob = c(0.025, 
                                                            0.975))
        sigMatNm[pw[2, i], pw[1, i]] <- sigDiff(p1, 
                                                p2)
        sigMatNm[pw[1, i], pw[2, i]] <- sigDiff(p2, 
                                                p1)
      }
      nmRelPlt[!sigMatNm] <- 0
    }
    if (plot_network) {
      if (stat == "d" && !is.null(outfile)) {
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; D method)", sep = ""))
        dev.off()
      }
      if (stat == "gst" && !is.null(outfile)) {
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Gst method)", sep = ""))
        dev.off()
      }
      if (stat == "Nm" && !is.null(outfile)) {
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Nm method)", sep = ""))
        dev.off()
      }
      if (stat == "all" && !is.null(outfile)) {
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; D method)", sep = ""))
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps)", sep = ""))
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Nm method)", sep = ""))
        dev.off()
      }
      if (stat == "d") {
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; D method)", sep = ""))
      }
      if (stat == "gst") {
        qgraph::qgraph(gRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Gst method)", sep = ""))
      }
      if (stat == "Nm") {
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Nm method)", sep = ""))
      }
      if (stat == "all") {
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; D method)", sep = ""))
        qgraph::qgraph(dRelPlt, nodeNames = sapply(dat$indnms, 
                                                   "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Gst Method)", sep = ""))
        qgraph::qgraph(nmRelPlt, nodeNames = sapply(dat$indnms, 
                                                    "[", 1), legend = TRUE, posCol = plot_col, 
                       label.color = plot_col, edge.labels = TRUE, 
                       curve = 2.5, mar = c(2, 2, 5, 5))
        title(paste("Significant relative migration network \n (", 
                    boots, " bootstraps; Nm method)", sep = ""))
      }
    }
  }
  if (boots != 0L) {
    if (stat == "d") {
      list(dRelMig = dRel, dRelMigSig = dRelPlt)
    }
    else if (stat == "gst") {
      list(gRelMig = gRel, gRelMigSig = gRelPlt)
    }
    else if (stat == "Nm") {
      list(nmRelMig = nmRel, nmRelMigSig = nmRelPlt)
    }
    else if (stat == "all") {
      list(dRelMig = dRel, dRelMigSig = dRelPlt, gRelMig = gRel, 
           gRelMigSig = gRelPlt, nmRelMig = nmRel, nmRelMigSig = nmRelPlt)
    }
  }
  else {
    if (stat == "d") {
      list(dRelMig = dRel)
    }
    else if (stat == "gst") {
      list(gRelMig = gRel)
    }
    else if (stat == "Nm") {
      list(nmRelMig = nmRel)
    }
    else if (stat == "all") {
      list(dRelMig = dRel, gRelMig = gRel, nmRelMig = nmRel)
    }
  }
}
