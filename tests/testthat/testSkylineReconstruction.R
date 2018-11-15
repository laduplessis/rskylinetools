context("Skyline Reconstruction")

# Load simulated data and example outputs
load("../testDataBackwardSkylines.RData")


test_that("Skyline reconstruction works on a grid (small sample)", {
    grid_t    <- seq(0, 10, length.out=5)
    expect_equal(reconstructBackwardSkylineSlow(linearSample$N, linearSample$t, grid_t), linearSampleExpectedGridded)
})


test_that("Skyline reconstruction works at specific times and boundary conditions (small sample)", {
    check_t <- c(1.5, 3, 4.2, 5.4, 7.2)
    expect_equal(reconstructBackwardSkylineSlow(linearSample$N, linearSample$t, check_t), linearSampleExpectedChecked)
})


test_that("Skyline reconstruction works on a bigger grid with more samples (periodic function) - Careful! Example output produced by function!", {
    grid_t    <- seq(0, 10, length.out=100)
    expect_equal(reconstructBackwardSkylineSlow(periodicSample$N, periodicSample$t, grid_t), periodicSampleExpectedGridded)
})


test_that("Skyline reconstruction works on an uneven grid with more samples (exponential function) - Careful! Example output produced by function!", {
  expect_equal(reconstructBackwardSkylineSlow(exponentialSample$N,exponentialSample$t, randomTimes), exponentialSampleExpectedRandom)
})

