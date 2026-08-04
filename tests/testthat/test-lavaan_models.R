test_that(".build_mimic generates valid syntax", {
    sel <- structure(
        list(
            variables = c("y1", "y2", "y3"),
            indicator_groups = list(1:3),
            direct_effects = character(0),
            K = 1L,
            fit = list(X = matrix(0, 10, 3,
                                  dimnames = list(NULL, c("(Intercept)", "group", "sex"))))
        ),
        class = "biomimic_selection"
    )

    syntax <- build_lavaan_model(sel, type = "mimic")
    expect_true(grepl("F1 =~ y1 \\+ y2 \\+ y3", syntax))
    expect_true(grepl("F1 ~ group \\+ sex", syntax))
})

test_that(".build_multigroup requires group_var", {
    sel <- structure(list(variables = c("y1", "y2", "y3")),
                     class = "biomimic_selection")
    expect_error(build_lavaan_model(sel, type = "multigroup"),
                 "group_var")
})

test_that(".build_direct adds direct-effect paths", {
    sel <- structure(
        list(
            variables = c("y1", "y2", "y3"),
            indicator_groups = list(1:3),
            direct_effects = c("y2"),
            K = 1L,
            fit = list(X = matrix(0, 10, 2,
                                  dimnames = list(NULL, c("(Intercept)", "group"))))
        ),
        class = "biomimic_selection"
    )

    syntax <- build_lavaan_model(sel, type = "direct")
    expect_true(grepl("y2 ~ group", syntax))
})
