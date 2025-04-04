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

