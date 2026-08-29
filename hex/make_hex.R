# ------------------------------------------------------------------
# orchaRd hex sticker
#
# Concept: an orchard plot drawn as a fruit tree. The vertical zero line
# of the orchard plot is the trunk, each moderator level is a branch
# (prediction interval = thin branch, confidence interval = thick branch,
# point estimate = the big fruit at the tip), and the individual effect
# sizes are the fruit clustered around each branch, sized by precision.
#
# Run from the package root:  Rscript hex/make_hex.R
# Outputs: hex/orchaRd_hex.png (hi-res), hex/orchaRd_hex.svg,
#          man/figures/logo.png (240 px, for README / pkgdown)
# ------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(showtext)
})

set.seed(2019)  # year the orchard plot was born

# ---- fonts ---------------------------------------------------------
sysfonts::font_add_google("Nunito", "nunito", regular.wt = 800, bold.wt = 900)
showtext_auto()
showtext_opts(dpi = 600)

# ---- palette -------------------------------------------------------
pal <- list(
  bg_top    = "#2E8B62",   # leafy green (top of hex)
  bg_bottom = "#0B3A2A",   # deep forest (bottom of hex)
  border    = "#F2D680",   # warm gold rim
  cream     = "#FBF6E9",   # branches / text
  cherry    = "#E8433F",
  apricot   = "#F7A038",
  lemon     = "#F9D34C",
  plum      = "#B76BD1"
)

# ---- hexagon geometry (pointy top, circumradius 1) -----------------
hex_pts <- function(r = 1) {
  a <- seq(90, 450, by = 60)[1:6] * pi / 180
  data.frame(x = r * cos(a), y = r * sin(a))
}
outer_hex <- hex_pts(1.00)
inner_hex <- hex_pts(0.935)

# ---- the "tree": branches (rows of the orchard plot) ---------------
# Each row: vertical position, mean estimate, CI / PI half-widths,
# number of fruit and its colour.
branches <- data.frame(
  y      = c( 0.44,  0.16, -0.12),
  mean   = c( 0.24, -0.20,  0.10),
  ci     = c( 0.08,  0.07,  0.06),
  pi     = c( 0.36,  0.30,  0.33),
  n      = c( 26,    22,    30),
  colour = c(pal$cherry, pal$apricot, pal$lemon),
  stringsAsFactors = FALSE
)

# A tiny deterministic beeswarm: bin on x, stack alternately on y
swarm <- function(x, bin = 0.045, dy = 0.036) {
  b   <- floor(x / bin)
  off <- ave(seq_along(x), b, FUN = function(i) {
    k <- seq_along(i) - 1
    ((k + 1) %/% 2) * ifelse(k %% 2 == 0, 1, -1)
  })
  off * dy
}

fruit <- do.call(rbind, lapply(seq_len(nrow(branches)), function(i) {
  b  <- branches[i, ]
  x  <- rnorm(b$n, b$mean, 0.19)
  x  <- pmax(pmin(x, b$mean + 0.42), b$mean - 0.42)
  # precision: points close to the mean tend to be bigger (as in real data)
  prec <- pmax(0.25, 1 - abs(x - b$mean) / 0.45) * runif(b$n, 0.55, 1.15)
  o  <- order(x)
  x  <- x[o]; prec <- prec[o]
  data.frame(x = x, y = b$y + swarm(x), size = prec, colour = b$colour)
}))

# ---- leaves at the branch tips -------------------------------------
leaf <- function(cx, cy, len = 0.09, wid = 0.035, angle = 0, n = 40) {
  t <- seq(0, 2 * pi, length.out = n)
  # almond / leaf outline
  x <- (len / 2) * cos(t)
  y <- (wid / 2) * sin(t) * (1 - 0.35 * cos(t))
  a <- angle * pi / 180
  data.frame(x = cx + x * cos(a) - y * sin(a),
             y = cy + x * sin(a) + y * cos(a))
}
leaf_spec <- data.frame(cx = c(-0.062, 0.066, 0.010),
                        cy = c( 0.690, 0.705, 0.760),
                        angle = c(52, -48, 95))
leaves <- do.call(rbind, lapply(seq_len(nrow(leaf_spec)), function(i)
  cbind(leaf(leaf_spec$cx[i], leaf_spec$cy[i], len = 0.105, wid = 0.042,
             angle = leaf_spec$angle[i]), id = i)))
veins <- transform(leaf_spec,
  x    = cx - 0.045 * cos(angle * pi / 180), y    = cy - 0.045 * sin(angle * pi / 180),
  xend = cx + 0.040 * cos(angle * pi / 180), yend = cy + 0.040 * sin(angle * pi / 180))

# ---- build the plot -------------------------------------------------
bg_grad <- grid::linearGradient(
  colours = c(pal$bg_top, pal$bg_bottom),
  x1 = 0.5, y1 = 1, x2 = 0.5, y2 = 0
)
glow <- grid::radialGradient(
  colours = c("#FFFFFF33", "#FFFFFF00"),
  cx1 = 0.5, cy1 = 0.62, r1 = 0, cx2 = 0.5, cy2 = 0.62, r2 = 0.55
)

p <- ggplot() +
  # hex body + rim
  geom_polygon(data = outer_hex, aes(x, y), fill = pal$border) +
  geom_polygon(data = inner_hex, aes(x, y), fill = bg_grad) +
  geom_polygon(data = inner_hex, aes(x, y), fill = glow) +
  # trunk (the zero line of the orchard plot)
  annotate("segment", x = 0, xend = 0, y = -0.30, yend = 0.66,
           colour = pal$cream, linewidth = 2.6, lineend = "round",
           alpha = 0.95) +
  annotate("segment", x = 0, xend = 0, y = -0.30, yend = 0.66,
           colour = pal$bg_bottom, linewidth = 0.55, lineend = "round",
           linetype = "42", alpha = 0.28) +
  # leaves on top of the trunk
  geom_polygon(data = leaves, aes(x, y, group = id),
               fill = "#7FD59A", colour = pal$cream, linewidth = 0.35) +
  geom_segment(data = veins, aes(x = x, y = y, xend = xend, yend = yend),
               colour = pal$cream, linewidth = 0.3, alpha = 0.8, lineend = "round") +
  # prediction intervals (thin branches)
  geom_segment(data = branches,
               aes(x = mean - pi, xend = mean + pi, y = y, yend = y),
               colour = pal$cream, linewidth = 1.3, lineend = "round",
               alpha = 0.9) +
  # confidence intervals (thick branches)
  geom_segment(data = branches,
               aes(x = mean - ci, xend = mean + ci, y = y, yend = y),
               colour = pal$cream, linewidth = 3.2, lineend = "round") +
  # the fruit (individual effect sizes)
  geom_point(data = fruit, aes(x, y, size = size, fill = colour),
             shape = 21, colour = alpha("#0B3A2A", 0.55), stroke = 0.35,
             alpha = 0.92) +
  # point estimates (the big fruit on each branch)
  geom_point(data = branches, aes(x = mean, y = y, fill = colour),
             shape = 21, size = 8.2, colour = pal$cream, stroke = 1.4) +
  geom_point(data = branches, aes(x = mean - 0.012, y = y + 0.012),
             shape = 16, size = 2.2, colour = alpha("white", 0.65)) +
  scale_fill_identity() +
  scale_size_continuous(range = c(1.6, 4.6), guide = "none") +
  # package name
  ggtext::geom_richtext(
    data = data.frame(x = 0, y = -0.525),
    aes(x, y),
    label = paste0("orcha<span style='color:", pal$border, "'>R</span>d"),
    family = "nunito", size = 15, colour = pal$cream,
    fill = NA, label.colour = NA, label.padding = grid::unit(0, "pt")
  ) +
  coord_fixed(xlim = c(-sqrt(3) / 2, sqrt(3) / 2), ylim = c(-1, 1),
              expand = FALSE, clip = "off") +
  theme_void() +
  theme(plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = margin(0, 0, 0, 0))

# ---- write out ------------------------------------------------------
w_in <- 2 * sqrt(3) / 2 * 2   # 3.46 in wide
h_in <- 2 * 2                 # 4.00 in tall (hexb.in ratio 1.73 : 2)

ragg::agg_png("hex/orchaRd_hex.png", width = w_in, height = h_in,
              units = "in", res = 600, background = "transparent")
print(p); invisible(dev.off())

svglite::svglite("hex/orchaRd_hex.svg", width = w_in, height = h_in, bg = "transparent")
print(p); invisible(dev.off())

# README / pkgdown logo (usethis::use_logo convention: 240 px wide)
img <- magick::image_read("hex/orchaRd_hex.png")
magick::image_write(magick::image_resize(img, "240x"), "man/figures/logo.png")

cat("Wrote hex/orchaRd_hex.png, hex/orchaRd_hex.svg, man/figures/logo.png\n")
