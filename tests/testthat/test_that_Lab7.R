context("ridgereg")

data("iris")

test_that("lenreg rejects errounous input", {
  expect_error(ridgereg_mod <- ridgereg$new(formula = Petal.Length~Sepdsal.Width+Sepal.Length, data=iris))
  expect_error(ridgereg_mod <- ridgereg$new(formula = Petal.Length~Sepdsal.Width+Sepal.Length, data=irfsfdis))
})

test_that("coef result is true", {
  library (MASS)
  lmR <- lm.ridge(Petal.Length ~ Species, iris)
  ridgereg <-  ridgereg$new(Petal.Length ~ Species, data=iris)
  expect_equal(round(as.vector(lmR$coef),2),round(as.vector(ridgereg$Coef),2))    
})