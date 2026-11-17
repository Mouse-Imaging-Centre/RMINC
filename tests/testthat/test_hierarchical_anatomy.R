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
