library(metafor)

library(devtools)

load_all()


# Base model for the leave one out
res <- metafor::rma.mv(lnrr, lnrr_vi,
                       random = list(~ 1 | es_ID,
                                     ~ 1 | paper_ID),
                       data = fish)

# ----------------------------------------
# This part is for testing: how would I run the code manually

papers_ids <- unique(fish$paper_ID)
res_all <- data.frame(left_out = papers_ids,
                      b        = NA,
                      ci_ub    = NA,
                      ci_lb    = NA)


for (i in seq_along(papers_ids)) {
  paper_id <- papers_ids[i]

  fish_sub <- fish[fish$paper_ID != paper_id, ]

  print(dim(fish_sub)) 

  res_sub <- metafor::rma.mv(lnrr, lnrr_vi,
                             random = list(~ 1 | es_ID,
                                           ~ 1 | paper_ID),
                             data = fish_sub)

  res_all$b[i] <- res_sub$b[1]
  res_all$ci_ub[i] <- res_sub$ci.ub[1]
  res_all$ci_lb[i] <- res_sub$ci.lb[1]
}

# This is how it shoul look like:
res_all

# ----------------------------------------


outputs <- run_leave1out(res, group = "paper_ID")

outputs

bla <- leave_one_out(res, group = "paper_ID")

names(bla)

names(res_all)


library(patchwork)
library(ggplot2)

p1 <- ggplot(res_all, aes(x = left_out, y = b)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lb, ymax = ci_ub), width = 0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip() +
  labs(title = "Manually")


p2 <- ggplot(bla, aes(x = left_out, y = b)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lb, ymax = ci_ub), width = 0.2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  coord_flip() + 
  labs(title = "Function")

p1 + p2
