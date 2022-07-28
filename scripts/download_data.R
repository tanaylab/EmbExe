# script for downloading and initializing the single-cell RNA-seq database for the metacell package
options(timeout = 1e9)


download_scrna_db <- function() {
    download.file("https://embexe.s3.eu-west-1.amazonaws.com/scrna_db_embexe.tar.gz", "scrna_db_embexe.tar.gz")

    system("tar -xzvf scrna_db_embexe.tar.gz")

    file.remove("scrna_db_embexe.tar.gz")
}

download_embexe_data <- function() {
    download.file("https://embexe.s3.eu-west-1.amazonaws.com/embexe_data.tar.gz", "embexe_data.tar.gz")

    if (!dir.exists("data")) {
        dir.create("data")
    }

    system("tar -xzvf embexe_data.tar.gz")

    file.remove("embexe_data.tar.gz")
}

download_full_data <- function() {
    download_embexe_data()
    download_scrna_db()
}


if (!dir.exists("figs")) {
    dir.create("figs")
}