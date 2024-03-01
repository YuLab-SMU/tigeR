test_that("multiplication works", {
  expect_no_error(Signature_calculation(SE=MEL_GSE78220))

  expect_no_error({
    train_set <- list(MEL_GSE78220, GBM_PRJNA482620, nonsqNSCLC_GSE93157)
    mymodel <- build_Model(Model='NB', SE=train_set, feature_genes=Stem.Sig,
                           response_NR = TRUE)
  })

  expect_no_error(build_Model(Model='SURV',
                              SE=MEL_GSE78220,
                              feature_genes=Stem.Sig,
                              PT_drop= TRUE))
  expect_error({ ##error: due to the abcent survival information of NSCLC_GSE135222
    mymodel <- build_Model(Model='SURV', SE=MEL_GSE93157,
                           feature_genes=Stem.Sig,
                           PT_drop= TRUE)
    test_Model(mymodel, SE=NSCLC_GSE135222, PT_drop=FALSE)
  })

  expect_no_error(Immunotherapy_Response(SE=MEL_GSE78220,
                                      geneSet ="CD274"))

  expect_no_error(plt_diff(SE=MEL_GSE78220,
                           gene='CD274',
                           type='Treatment'))
  expect_no_error(
    Immunotherapy_Response(SE=MEL_GSE78220, geneSet ="CD274")
  )

  expect_warning(
    Immunotherapy_Response(SE=MEL_GSE93157, geneSet ="CD274")
  )
  expect_no_error(plt_diff(SE=MEL_GSE78220,gene='CD274',type='Treatment'))
  expect_no_error(plt_surv(SE=MEL_GSE78220,gene=c('CD274')))
  expect_no_error(CIBERSORT(sig_matrix=LM22,SE=GBM_PRJNA482620,
                            perm=10, QN=T, style = 'elegant'))
  expect_no_error(CIBERSORT(sig_matrix=LM22,SE=GBM_PRJNA482620,
                            perm=10, QN=T, style = 'raw'))
})


