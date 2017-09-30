if (require("testthat")) {
  test_check("TreeWAS")
} else {
  warning("testthat not available. Skipping unittests!")
}
