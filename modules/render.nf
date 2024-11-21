process BOOK_RENDER {
    tag "Rendering report from analysis..."
    stageInMode 'copy'

    publishDir "$params.outdir/book", mode:'copy', pattern: "_book/*"
    publishDir "$params.outdir/book", mode:'copy', pattern: "report.zip"

    memory '16 GB'
    cpus 1

    input:
    path scripts
    path book_assets
    path quarto_assets



    output:
    path "_book/*"
    path "report.zip"


    script:
    """
    unzip -o '*.zip'
    Rscript filter_yaml.R
    quarto render index.qmd
    zip -r report.zip _book
    """

}


