plot_signal <-
function(signal_map) {
  ch_1 <- signal_map$nucl
  ch_2 <- signal_map$cyto
  ch_ori <- signal_map$input
  my_x <- ncol(ch_1)
  my_y <- nrow(ch_1)
  ch_0 <- matrix(0, ncol = my_x, nrow = my_y)
  ch_1[ch_1 == 0] <- NA
  ch_2[ch_2 == 0] <- NA
  #
  #faux palettes
  greys <- colorpanel(n = 15, low = "black", high = "grey99")
  cyans <- colorpanel(n = 15, low = "black", high = "cyan")
  greens <- colorpanel(n = 15, low = "black", high = "limegreen")
  #
  curr_par <- par(no.readonly = T)
  par(mfrow = c(2,2))
  image(ch_ori, col = greys, add = F, axes = F, main = "input image", useRaster = T)
  #
  image(ch_0, col = "black", add = F, axes = F, main = "faux colors", useRaster = T)
  image(ch_1, col = greens, add = T, axes = F, useRaster = T)
  image(ch_2, col = cyans, add = T, axes = F, useRaster = T)
  #
  image(ch_0, col = "black", add = F, axes = F, main = "nucl. signal", useRaster = T)
  image(ch_1, col = greys, add = T, axes = F, useRaster = T)
  #
  image(ch_0, col = "black", add = F, axes = F, main = "cyto. signal", useRaster = T)
  image(ch_2, col = greys, add = T, axes = F, useRaster = T)
  #
  par(curr_par)
}
