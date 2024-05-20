test_that("visualize PTM-SEA result works", {
  expect_true(ggplot2::is.ggplot(visualize_PTMSEA(test_path("testdata", "ccRCC-combined.gct"), "Tumor_vs_NAT_diff")))
  expect_true(ggplot2::is.ggplot(visualize_PTMSEA(test_path("testdata", "ccRCC-combined.gct"),
                   "Tumor_vs_NAT_diff",
                   selected_concepts=c("KINASE-PSP_PKCA/PRKCA", "KINASE-PSP_PKCB/PRKCB", "KINASE-PSP_PKCG/PRKCG", "KINASE-PSP_PKCB_iso2/PRKCB", "KINASE-PSP_PKCD/PRKCD", "KINASE-PSP_PKCI/PRKCI", "KINASE-PSP_PKCT/PRKCQ",
                                       "KINASE-PSP_PKCH/PRKCH", "KINASE-PSP_PKCE/PRKCE", "KINASE-PSP_PKCZ/PRKCZ"))))
})
