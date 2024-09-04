import groovy.json.JsonOutput

book_assets = channel.fromPath( 'assets/book_assets/*').collect()

include { QUALITY_CONTROL } from './modules/qc.nf'
include { PACKAGE_IMAGES } from './modules/create_seurat.nf'
include { CREATE_SEURAT_VISIUM } from './modules/create_seurat.nf'
include { DIMENSION_REDUCTION } from './modules/dim_reduce.nf'
include { INTEGRATE_DATASETS } from './modules/integrate_datasets.nf'
include { MULTIMODAL_INTEGRATION } from './modules/multimodal_integration.nf'
include { CLUSTERING } from './modules/clustering.nf'
include { SPATIAL_FEATURES } from './modules/clustering.nf'
include { CELL_ANNOTATION } from './modules/cell_annotation.nf'
include { SPOT_DECONVOLUTION } from './modules/spot_deconvolution.nf'
include { ISOLATE_GRANULOMAS } from './modules/isolate_granulomas.nf'
include { BOOK_RENDER } from './modules/render.nf'


scripts_ch = channel.fromPath("bin/*")
             .collect()


workflow {

    Channel.fromPath(  "${params.input}/**/outs/filtered_*bc_matrix.h5", checkIfExists : true )
    | map { it -> tuple( it.parent.parent.name, it)}
    | set { filt_counts_ch }

    Channel.fromPath(  "${params.input}/**/outs/spatial/*", checkIfExists : true )
    | map { it -> tuple( it.parent.parent.parent.name, it)}
    | groupTuple
    | set { image_ch }

    PACKAGE_IMAGES( image_ch)
    | set { pkg_image_ch }


    Channel.fromPath(  "${params.input}/**/outs/raw_*bc_matrix.h5", checkIfExists : true )
    | map { it -> tuple( it.parent.parent.name, it)}
    | join(filt_counts_ch)
    | join(pkg_image_ch)
    | set { counts_ch }

 
    CREATE_SEURAT_VISIUM(counts_ch, 
                         params.pipeline, 
                         params.create_seurat)
    | set { seurat_ch }


    if (params.matrix == "filtered") {
        seurat_ch.filt
        | collect 
        | set { seurat_combined_ch}
        
    } else {
        seurat_ch.raw
        | collect 
        | set { seurat_combined_ch}
    }

   


    QUALITY_CONTROL(seurat_combined_ch,
                    params.qc_script,
                    book_assets,
                    params.pipeline,
                    params.metadata,
                    params.genome,
                    params.nCount_Spatial_min,
                    params.nCount_Spatial_max,
                    params.nFeature_Spatial_min,
                    params.nFeature_Spatial_max,
                    params.sample_exclusion).set {qc_ch}
 
    quarto_ch = qc_ch.quarto

   

    DIMENSION_REDUCTION(qc_ch.seurat,
                       params.dimreduc_script,
                       book_assets,
                       params.pipeline,
                       params.umap_ndims)
    | set {umap_ch }
    quarto_ch = quarto_ch.concat(umap_ch.quarto)
    
    
    if (params.integrate_datasets ) {
        INTEGRATE_DATASETS(umap_ch.seurat,
                        params.int_script,
                        book_assets,
                        params.pipeline,
                        params.integrate_datasets,
                        params.integration_method,
                        params.integrate_by,
                        params.view_batch,
                        params.umap_ndims)
        | set { int_ch }
        quarto_ch = quarto_ch.concat(int_ch.quarto)
        
    } else {
        int_ch = umap_ch
    }


    if (params.pipeline == 'cite') {
        MULTIMODAL_INTEGRATION(int_ch.seurat,
                            params.int_multimod_script,
                            book_assets,
                            params.pipeline,
                            params.umap_ndims,
                            params.integrate_datasets)
        | set { intwnn_ch }
        quarto_ch = quarto_ch.concat(intwnn_ch.quarto)
        
    } else {
        intwnn_ch = int_ch
    }

    CLUSTERING(intwnn_ch.seurat,
                params.clustering_script,
                book_assets,
                params.pipeline,
                params.clustering_res,
                params.integrate_datasets,
                params.outcomes)
    | set { cluster_ch }

    quarto_ch = quarto_ch.concat(cluster_ch.quarto)



    SPOT_DECONVOLUTION(cluster_ch.seurat,
                        params.annotation_script,
                        book_assets,
                        params.pipeline,
                        params.integrate_datasets)
    | set { deconv_ch }

    quarto_ch = quarto_ch.concat(deconv_ch.quarto)

   
    SPATIAL_FEATURES(cluster_ch.seurat,
                params.spat_var_script,
                book_assets,
                params.pipeline,
                params.integrate_datasets,
                params.selected_var_features)
    | set { spat_var_ch }

    quarto_ch = quarto_ch.concat(spat_var_ch.quarto)

    /*
    ISOLATE_GRANULOMAS(deconv_ch.seurat,
                        params.annotation_script,
                        book_assets,
                        params.pipeline,
                        params.integrate_datasets)
    | set { deconv_ch }

    quarto_ch = quarto_ch.concat(deconv_ch.quarto)


    CELL_ANNOTATION(cluster_ch.seurat,
                    params.annotation_script,
                    book_assets,
                    params.pipeline,
                    params.integrate_datasets,
                    params.annotation_resolution,
                    params.run_azimuth,
                    params.azimuth_reference,
                    params.run_singler,
                    params.singler_reference,
                    params.run_scimilarity,
                    params.scimilarity_reference,
                    params.selected_method,
                    params.markers_rna,
                    params.markers_adt,
                    params.outcomes)
    | set { cell_ann_ch }

    quarto_ch = quarto_ch.concat(cell_ann_ch.quarto)

    */


    BOOK_RENDER(scripts_ch, 
                     book_assets,
                     quarto_ch.collect())
    
    
}


workflow.onComplete {
    log.info ( workflow.success ? "\nDone!" : "Oops .. something went wrong" )
}
