comps <- compartmentData(c("MAPK1", "MAPK3"))$enrichment
comps_t <- compartmentData(c("MAPK1", "MAPK3"), aspect ="Endosomal system")$enrichment

test_that("returns plot", {
  expect_equal(class(runSubcellulaRvis(comps, "blue", "red")), c("gg", "ggplot"))
  expect_equal(class(
    runSubcellulaRvis(
      comps_t,
      "blue",
      "red",
      aspect ="Endosomal system"
    )
  ), c("gg", "ggplot"))
  expect_equal(class(runSubcellulaRvis(comps, "blue", "red", text_size = 15)),
               c("gg", "ggplot"))
  expect_equal(class(runSubcellulaRvis(comps, "blue", "red", legend = F)),
               c("gg", "ggplot"))
  expect_equal(class(runSubcellulaRvis(comps, "blue", "red", legend.pos = "left")),
               c("gg", "ggplot"))
  expect_equal(class(runSubcellulaRvis(comps, "blue", "red", labels = T)),
               c("gg", "ggplot"))
})

test_that("wrong data returns helpful error", {
  expect_error(runSubcellulaRvis(data.frame("uh", "oh"), "blue", "red"))
  expect_error(runSubcellulaRvis(comps[, 1:3], "blue", "red"))
  expect_error(runSubcellulaRvis(comps[1:3,], "blue", "red"))
  expect_error(runSubcellulaRvis(comps, "blue", "red", aspect ="Endosomal system"))
  expect_error(runSubcellulaRvis(compartmentData(c("MAPK1", "MAPK3"),
                                                 trafficking = T), 
                                 "blue", "red"))
  
})
