process FILTER{
    publishDir = [
        path: { "results/seurat" },
        mode: 'symlink',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]



    input: 
        path rds
    
    output:
        path("*.rds"), emit:rds

    script:
        """
            filter.R ${rds}
        """
}

process SEURAT{
    publishDir = [
        path: { "results/seurat" },
        mode: 'symlink',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]



    input: 
        path rds
    
    output:
        path("*.rds"), emit:rds

    script:
        """
            seurat.R ${rds}
        """
}

process MARKER{
    publishDir = [
        path: { "results/marker" },
        mode: 'symlink',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
    clusterOptions = "-p low"
    input: 
        path txt
    output:
        path("*.rds"), emit:rds

    script:
        """
           marker.R   ${txt} 
        """
}


process FEATURE_PLOT {
    publishDir = [
        path: { "results/feature_plot" },
        mode: 'copy',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]
    clusterOptions = "-p low"
    input: 
        path seuratObjPath
        path markerPath
    output:
        path("*/*.png"), emit:png

    script:
        """
           featurePlot.R   ${seuratObjPath}  ${markerPath}
        """
}


process CELL_ASSIGN{
    publishDir = [
        path: { "results/assign" },
        mode: 'symlink',
        saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
    ]


    input: 
        path rds
        path marker 
    output:
        path("*.rds"), emit:rds

    script:
        """
            cellassign.R  ${rds} ${marker}
        """
}



input = file(params.input) 
marker = file(params.marker) 
workflow{
    ch_input = Channel.fromPath( input )
    // // ch_input =  Channel.fromPath( 'pmbc_1683691139445.rds' )
    FILTER(ch_input)
    // // ch_input.view()
    SEURAT(FILTER.out.rds)
    MARKER(marker) 
    FEATURE_PLOT(SEURAT.out.rds, MARKER.out.rds)
    CELL_ASSIGN(FILTER.out.rds,MARKER.out.rds)  
}


