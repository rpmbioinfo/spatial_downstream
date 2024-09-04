process DIMENSION_REDUCTION {
    tag "Performing dimension reduction..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_dim_reduced.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"

    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val umap_ndims
 
 

    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_dim_reduced.RDS", emit:seurat
    path rmd, emit:script

    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P umap_ndims:${umap_ndims}

    bash chapter_package.sh "${rmd.baseName}"
    """

}