df <- iris
df$spp_id <- as.numeric(iris$Species)

test_that(".column_exists works", {
  expect_true(.column_exists(data = df, column_name = "spp_id"))
  expect_false(.column_exists(data = df, column_name = "foo"))
})

test_that(".column_is_categorical works", {
  expect_true(.column_is_categorical(col = df$Species))
  expect_true(.column_is_categorical(col = df$spp_id))
  expect_false(.column_is_categorical(col = df$Sepal.Width))
})

test_that(".is_group_valid works", {
  expect_no_error(.is_group_valid(data = df, column_name = "Species"))
  expect_no_error(.is_group_valid(data = df, column_name = "spp_id"))
  expect_error(.is_group_valid(data = df, column_name = "Sepal.Width"))
  expect_error(.is_group_valid(df$Sepal.Width))
})

test_that(".is_group_valid rejects categorical columns containing NAs", {
  df_na <- df
  df_na$Species[c(1, 5, 10)] <- NA
  expect_error(.is_group_valid(data = df_na, column_name = "Species"),
               "can't have NAs", fixed = TRUE)
})

# --- .is_metafor_object / .is_model_valid coverage ---

data(lim)
lim$vi <- (1/sqrt(lim$N - 3))^2
lim_MR <- metafor::rma.mv(yi = yi, V = vi, mods = ~ Phylum - 1,
                          random = list(~1 | Article, ~1 | Datapoint), data = lim)

test_that(".is_metafor_object recognises metafor classes and rejects others", {
  expect_true(.is_metafor_object(lim_MR))
  expect_false(.is_metafor_object(stats::lm(yi ~ Phylum, data = lim)))
  expect_false(.is_metafor_object(list()))
})

test_that(".is_model_valid errors on missing, NULL, or non-metafor objects", {
  expect_true(.is_model_valid(lim_MR))
  expect_error(.is_model_valid(),     "Incorrect argument 'model'", fixed = TRUE)
  expect_error(.is_model_valid(NULL), "Incorrect argument 'model'", fixed = TRUE)
  expect_error(.is_model_valid(stats::lm(yi ~ Phylum, data = lim)),
               "Incorrect argument 'model'", fixed = TRUE)
})

