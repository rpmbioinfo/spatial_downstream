process ISOLATE_GRANULOMAS {
    tag "Reclustering granulomas ..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_granulomas_clustered.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"


    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val clustering_res
    val integrate_datasets
    val outcomes
 

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_granulomas_clustered.RDS", emit:seurat
    path "_freeze", emit: freeze

    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P clustering_res:${clustering_res} \
                    -P integrate_datasets:${integrate_datasets} \
                    -P outcomes:"${outcomes}"


    bash chapter_package.sh "${rmd.baseName}"
    """

}