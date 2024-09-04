process MULTIMODAL_INTEGRATION {
    tag "Performing multimodal integration..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_integrated_wnn.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"
    
    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val umap_ndims
    val integrate_datasets
 

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_integrated_wnn.RDS", emit:seurat


    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P umap_ndims:${umap_ndims} \
                    -P integrate_datasets:${integrate_datasets}


    bash chapter_package.sh "${rmd.baseName}"
    """

}