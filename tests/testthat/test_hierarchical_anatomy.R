context("makeMICeDefsHierachical")

# Regression test for issue #335: label propagation failed on R >= 4.2 with
# "'length = 2' in coercion to 'logical(1)'" from data.tree::Aggregate.

make_abi_json <- function(path) {
  abi <- list(msg = list(list(
    id = 1, name = "root", acronym = "root", color_hex_triplet = "FFFFFF",
    children = list(list(
      id = 2, name = "Brain", acronym = "B", color_hex_triplet = "AAAAAA",
      children = list(
        list(id = 3, name = "Cortex", acronym = "Ctx",
             color_hex_triplet = "BBBBBB", children = list()),
        list(id = 4, name = "Cerebellum", acronym = "Cb",
             color_hex_triplet = "CCCCCC", children = list())
      )
    ))
  )))
  writeLines(rjson::toJSON(abi), path)
}

make_defs_csv <- function(path) {
  defs <- data.frame(
    Structure = c("cortex", "cerebellum"),
    right.label = c(2, 3),
    left.label = c(1, 3),
    ABI = c("Cortex", "Cerebellum")
  )
  write.csv(defs, path, row.names = FALSE)
}

test_that("makeMICeDefsHierachical propagates labels up the tree", {
  skip_if_not_installed("data.tree")
  skip_if_not_installed("rjson")

  abijson <- tempfile(fileext = ".json")
  defscsv <- tempfile(fileext = ".csv")
  on.exit(unlink(c(abijson, defscsv)))
  make_abi_json(abijson)
  make_defs_csv(defscsv)

  hdefs <- makeMICeDefsHierachical(defscsv, abijson)

  expect_setequal(hdefs$label, c(1, 2, 3))
  expect_setequal(data.tree::FindNode(hdefs, "Brain")$label, c(1, 2, 3))
  expect_setequal(data.tree::FindNode(hdefs, "Cortex")$label, c(1, 2))
  expect_equal(data.tree::FindNode(hdefs, "left Cortex")$label, 1)
  expect_equal(data.tree::FindNode(hdefs, "right Cortex")$label, 2)
  expect_equal(data.tree::FindNode(hdefs, "Cerebellum")$label, 3)
  expect_equal(data.tree::FindNode(hdefs, "Cortex")$acronym, "Ctx")
  expect_equal(data.tree::FindNode(hdefs, "left Cortex")$acronym, "lCtx")
})

make_slash_files <- function(abijson, defscsv, abi_names) {
  abi <- list(msg = list(list(
    id = 1, name = "root", acronym = "root", color_hex_triplet = "FFFFFF",
    children = list(list(
      id = 2, name = "Brain", acronym = "B", color_hex_triplet = "AAAAAA",
      children = list(
        list(id = 3, name = "Visual area, layer 2/3", acronym = "VIS2/3",
             color_hex_triplet = "BBBBBB", children = list()),
        list(id = 4, name = "Cerebellum", acronym = "Cb",
             color_hex_triplet = "CCCCCC", children = list())
      )
    ))
  )))
  writeLines(rjson::toJSON(abi), abijson)
  defs <- data.frame(
    Structure = c("visual 2/3", "cerebellum"),
    right.label = c(2, 3),
    left.label = c(1, 3),
    ABI = abi_names
  )
  write.csv(defs, defscsv, row.names = FALSE)
}

test_that("makeMICeDefsHierachical accepts names with a slash (#221)", {
  skip_if_not_installed("data.tree")
  skip_if_not_installed("rjson")

  abijson <- tempfile(fileext = ".json")
  defscsv <- tempfile(fileext = ".csv")
  on.exit(unlink(c(abijson, defscsv)))
  make_slash_files(abijson, defscsv, c("Visual area, layer 2/3", "Cerebellum"))

  hdefs <- makeMICeDefsHierachical(defscsv, abijson)

  vis <- data.tree::FindNode(hdefs, "Visual area, layer 2 and 3")
  expect_false(is.null(vis))
  expect_setequal(vis$label, c(1, 2))
  expect_equal(vis$acronym, "VIS2/3")
  # no extra level is made from the part after the slash
  expect_null(data.tree::FindNode(hdefs, "3"))
  expect_setequal(hdefs$label, c(1, 2, 3))
})

test_that("makeMICeDefsHierachical names ABI names that are not in the JSON (#304)", {
  skip_if_not_installed("data.tree")
  skip_if_not_installed("rjson")

  abijson <- tempfile(fileext = ".json")
  defscsv <- tempfile(fileext = ".csv")
  on.exit(unlink(c(abijson, defscsv)))
  make_slash_files(abijson, defscsv, c("Visual area, layer 2/3", "Cerebelum"))

  expect_error(
    makeMICeDefsHierachical(defscsv, abijson),
    "not in the ABI hierarchy.*\"Cerebelum\""
  )
})

test_that("hanatLm and hanatFDR store one vector of results in each node", {
  skip_if_not_installed("data.tree")
  skip_if_not_installed("rjson")

  abijson <- tempfile(fileext = ".json")
  defscsv <- tempfile(fileext = ".csv")
  on.exit(unlink(c(abijson, defscsv)))
  make_abi_json(abijson)
  make_defs_csv(defscsv)
  hdefs <- makeMICeDefsHierachical(defscsv, abijson)

  set.seed(2)
  gfh <- data.frame(IV = factor(rep(c("a", "b"), each = 10)))
  vols <- cbind(
    rnorm(20, 10) + (gfh$IV == "b"),
    rnorm(20, 10),
    rnorm(20, 20) - (gfh$IV == "b")
  )
  attr(vols, "anatIDs") <- c(1, 2, 3)
  hvols <- addVolumesToHierarchy(hdefs, vols)

  verboseRun(hlm <- hanatLm(~IV, gfh, hvols))
  cortex <- data.tree::FindNode(hlm, "Cortex")
  expect_null(dim(cortex$model))
  expect_true("tvalue-IVb" %in% names(cortex$model))
  verboseRun(
    expected <- anatLm(~IV, gfh, matrix(cortex$volumes, ncol = 1))
  )
  expect_equal(
    unname(cortex$model[["tvalue-IVb"]]),
    unname(unclass(expected)[1, "tvalue-IVb"])
  )

  verboseRun(hq <- hanatFDR(hlm))
  expect_null(dim(data.tree::FindNode(hq, "Cortex")$aFDR))
  expect_false(is.null(attr(hq, "thresholds")))
})
