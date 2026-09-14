#!/usr/bin/env Rscript
# User explicitly authorized qs and required dependencies in the project library.
args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", args[grepl("^--file=", args)])
here <- dirname(normalizePath(script))
lib <- file.path(here, "r-library")
sources <- file.path(here, "external", "r-package-sources")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
dir.create(sources, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(timeout = 300)
Sys.setenv(MAKEFLAGS = "-j2")
# Use an explicit compiler standard within this process only.
makevars <- file.path(here, "inputs", "qs_Makevars")
writeLines(c("CXX = g++ -std=gnu++17", "CXX_STD = CXX17"), makevars)
Sys.setenv(R_MAKEVARS_USER = makevars)
packages <- c("RcppParallel", "RApiSerialize", "stringfish", "qs")
versions <- c("6.2.1", "0.1.4", "0.16.0", "0.27.3")
# stringfish 0.19.2 changed the C++ interface used by the archived qs release.
urls <- c(paste0("https://cloud.r-project.org/src/contrib/", packages[1:2], "_", versions[1:2], ".tar.gz"),
          "https://cran.r-project.org/src/contrib/Archive/stringfish/stringfish_0.16.0.tar.gz",
          "https://cran.r-project.org/src/contrib/Archive/qs/qs_0.27.3.tar.gz")
for (i in seq_along(packages)) {
    target <- file.path(sources, basename(urls[i]))
    if (!file.exists(target)) download.file(urls[i], target, mode = "wb")
    description <- file.path(lib, packages[i], "DESCRIPTION")
    if (!file.exists(description) || read.dcf(description, fields = "Version")[1] != versions[i]) {
        install.packages(target, lib = lib, repos = NULL, type = "source", Ncpus = 2)
    }
    stopifnot(requireNamespace(packages[i], quietly = TRUE),
              as.character(packageVersion(packages[i])) == versions[i],
              normalizePath(find.package(packages[i])) == normalizePath(file.path(lib, packages[i])))
}
example <- list(x = c(0, 1, NA_real_), label = c("DA", "other neurons"), f = factor(c("HC", "PD")))
roundtrip <- tempfile(fileext = ".qs")
qs::qsave(example, roundtrip)
stopifnot(identical(example, qs::qread(roundtrip, use_alt_rep = FALSE)))
unlink(roundtrip)
manifest <- data.frame(package = packages, version = versions, source_url = urls,
                       source_md5 = unname(tools::md5sum(file.path(sources, basename(urls)))),
                       installed_path = vapply(packages, find.package, character(1)))
write.csv(manifest, file.path(here, "inputs", "qs_installation_manifest.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(here, "inputs", "qs_R_sessionInfo.txt"))
cat("PASS: isolated qs installation and object round-trip\n")
print(manifest[, c("package", "version", "installed_path")], row.names = FALSE)
