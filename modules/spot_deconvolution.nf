process SPOT_DECONVOLUTION {
    tag "Performing spatial deconvolution..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_annotated.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"


    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val integrate_datasets
 

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_annotated.RDS", emit:seurat
    path "_freeze", emit: freeze

    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P integrate_datasets:${integrate_datasets} 


    bash chapter_package.sh "${rmd.baseName}"
    """

}

