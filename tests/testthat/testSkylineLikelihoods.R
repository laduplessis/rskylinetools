context("Generalized Skyline Plot likelihood")

test_that("skyline_lk() gives the same results as ape::skyline (homochronous trees only)", {

  # Classic skyline, small tree (groupSizes = 1)
  tree.test <- ape::read.tree(text = "((D4Philip56:30.0,(D4Philip64:23.0,D4Philip84:23.0):7.0):10.0,(D4SLanka78:25.0,(D4Thai78:11.0,D4Thai84:11.0):14.0):15.0);") # load tree
  sk1 <- ape::skyline(tree.test)
  expect_equal(skyline_lk(tree.test, sk1$population.size), sk1$logL)

  # Classic skyline, big tree (groupSizes = 1)
  data("hivtree.newick", package="ape") # example tree in NH format
  tree.hiv <- ape::read.tree(text = hivtree.newick) # load tree
  sk1      <- ape::skyline(tree.hiv)
  expect_equal(skyline_lk(tree.hiv, sk1$population.size), sk1$logL, tolerance=1e-6)

  # Generalized skyline, big tree
  cl2      <- ape::collapsed.intervals(ape::coalescent.intervals(tree.hiv),0.0119)
  sk2      <- ape::skyline(cl2)
  groupSizes <- c()
  for (i in 1:max(cl2$collapsed.interval)) {
    groupSizes <- c(groupSizes, sum(cl2$collapsed.interval == i))
  }
  expect_equal(skyline_lk(tree.hiv, sk2$population.size, groupSizes), sk2$logL, tolerance=1e-6)
})


test_that("skyline_lk() gives the same results as the BSP", {
  expect_equal(skyline_lk(ape::read.tree("../testTree1.tree"), c(1,2), c(2,3), full=FALSE), -305.57944154167984)
  expect_equal(skyline_lk(ape::read.tree("../testTree2.tree"), c(1,2), c(2,2), full=FALSE), -24.386294361119894)
  expect_equal(skyline_lk(ape::read.tree("../testTree3.tree"), c(1,2), c(2,3), full=FALSE), -95.57944154167984)
  expect_equal(skyline_lk(ape::read.tree("../testTree4.tree"), c(1,2), c(2,3), full=FALSE), -93.57944154167984)
  expect_equal(skyline_lk(ape::read.tree("../testTree5.tree"), c(1,2), c(2,3), full=FALSE), -103.57944154167984)
})


test_that("esp_lk() gives the same results as those calculated by hand (and by BESP)", {
  expect_equal(esp_lk(ape::read.tree("../testTree1.tree"), c(1,2), c(8,3), 1, 6, full=FALSE), -305.57944154167984)
  expect_equal(esp_lk(ape::read.tree("../testTree2.tree"), c(1,2), c(6,3), 2, 5, full=FALSE), -56.2274112777602)
  expect_equal(esp_lk(ape::read.tree("../testTree3.tree"), c(2,3), c(6,5), c(1,2,3), c(3,2,1), full=FALSE), -183.8716748273099)
  expect_equal(esp_lk(ape::read.tree("../testTree4.tree"), c(3,2,1), c(3,4,4), c(2,3), c(4,2), full=FALSE), -205.234349834419)
  expect_equal(esp_lk(ape::read.tree("../testTree5.tree"), c(2,3,4), c(5,2,4), c(2,3,4), c(3,2,1), full=FALSE), -289.280186700424)
})


