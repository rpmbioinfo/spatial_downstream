process CREATE_SEURAT_VISIUM {
    tag "Making Seurat objects... (${sample}"

    stageInMode 'copy'
    publishDir "$params.outdir/seurat", mode:'copy', pattern: "*.RDS"
    publishDir "$params.outdir/quarto", mode:'copy', pattern: "*.Rmd"



    memory '16 GB'
    cpus 2

    input:
    tuple val(sample), path(raw_matrix), path(filt_matrix), path(image)
    val pipeline
    path rmd


    output:
    path "*_raw_counts_seurat.RDS", emit: raw
    path "*_filtered_counts_seurat.RDS", emit: filt
    path rmd, emit:script


    script:
    """
    mkdir -p /home/data/spatial
    tar -xvzf ${image} --directory /home/data/spatial  
    cp *bc_matrix.h5 /home/data/
    Rscript -e 'rmarkdown::render("${rmd}", params = list(pipeline = "${pipeline}", sample = "${sample}"))'
    """

}

process PACKAGE_IMAGES {
    tag "Staging images... (${sample}"

    stageInMode 'copy'


    memory '16 GB'
    cpus 2

    input:
    tuple val(sample), path(files)



    output:
    tuple val(sample), path("*.tar.gz")


    script:
    """
    tar -hczvf spatial.tar.gz ${files}
    """

}
