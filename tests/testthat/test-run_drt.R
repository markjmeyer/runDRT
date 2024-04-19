test_that("test G is vector", {
  expect_error(run_drt(t(fda::CanadianWeather$dailyAv[,,'Temperature.C']) ~ 1))
})

test_that("test default", {
  expect_s3_class(run_drt(t(fda::CanadianWeather$dailyAv[,,'Temperature.C']) ~ fda::CanadianWeather$region), "htest")
})

test_that("test X is matrix or array", {
  expect_error(run_drt(1 ~ fda::CanadianWeather$region))
})

test_that("test method typo stuff.rank", {
  expect_error(run_drt(t(fda::CanadianWeather$dailyAv[,,'Temperature.C']) ~ fda::CanadianWeather$region, method = 'stuff.rank'))
})

test_that("test method typo suf.rank", {
  expect_error(run_drt(t(fda::CanadianWeather$dailyAv[,,'Temperature.C']) ~ fda::CanadianWeather$region, method = 'suf.rank'))
})

test_that("test method typo ag.rank", {
  expect_error(run_drt(t(fda::CanadianWeather$dailyAv[,,'Temperature.C']) ~ fda::CanadianWeather$region, method = 'ag.rank'))
})

test_that("test method = suff.rank", {
  expect_s3_class(run_drt(t(fda::CanadianWeather$dailyAv[,,'Temperature.C']) ~ fda::CanadianWeather$region, method = 'suff.rank'), "htest")
})

test_that("test method = avg.rank", {
  expect_s3_class(run_drt(t(fda::CanadianWeather$dailyAv[,,'Temperature.C']) ~ fda::CanadianWeather$region, method = 'avg.rank'), "htest")
})
