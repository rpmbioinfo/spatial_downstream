process QUALITY_CONTROL {
    tag "Performing Quality Control..."

    stageInMode 'copy'
    publishDir "$params.outdir/seurat", mode:'copy', pattern: "seurat_after_qc.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "_quarto_full.yml"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "index.qmd"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "run_funcs.R"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "favicon-32x32.png"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "rpm_en_logo_lowres.jpg"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.qmd"
    publishDir "$params.outdir/quarto/", mode:'copy', pattern: "_freeze/${rmd.baseName}/*"


    input:
    path seurat_solo
    path rmd
    path book_assets
    val pipeline
    path metadata
    val genome
    val nCount_Spatial_min
    val nCount_Spatial_max
    val nFeature_Spatial_min
    val nFeature_Spatial_max
    val sample_exclusion



    output:
    path "*_freeze.zip", emit: quarto
    path "seurat_after_qc.RDS", emit:seurat
    path rmd, emit:script
    path book_assets, emit: assets


    script:
    """
    Rscript filter_yaml.R ${rmd}
    quarto render ${rmd} \
                    -P pipeline:"${pipeline}" \
                    -P metadata:"${metadata}" \
                    -P genome:"${genome}" \
                    -P nCount_Spatial_min:"${nCount_Spatial_min}" \
                    -P nCount_Spatial_max:"${nCount_Spatial_max}" \
                    -P nFeature_Spatial_min:"${nFeature_Spatial_min}" \
                    -P nFeature_Spatial_max:"${nFeature_Spatial_max}" \
                    -P sample_exclusion:"${sample_exclusion}" 


    bash chapter_package.sh "${rmd.baseName}"
    """

}
