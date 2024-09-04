process CLUSTERING {
    tag "Performing clustering..."
    stageInMode 'copy'

    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_clustered.RDS"
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
    path "seurat_clustered.RDS", emit:seurat
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

process SPATIAL_FEATURES {
    tag "Finding spatially variable features..."
    stageInMode 'copy'

    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"


    input:
    path seurat
    path rmd
    path book_assets
    val pipeline
    val integrate_datasets
    val selected_var_features
 

    output:
    path "*_freeze.zip", emit: quarto
    path "_freeze", emit: freeze

    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P seurat:"${seurat}" \
                    -P pipeline:"${pipeline}" \
                    -P integrate_datasets:${integrate_datasets} \
                    -P selected_var_features:"${selected_var_features}"


    bash chapter_package.sh "${rmd.baseName}"
    """

}


