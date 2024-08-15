#Average bigwig

bigWigMerge S2DOX221117-1.coverage.bw S2DOX221117-2.coverage.bw S2DOX221117-3.coverage.bw S2DOX221117-4.coverage.bw S2DOX221117-5.coverage.bw S2DOX221117-6.coverage.bw S2DOX221117.avg.coverage.bedGraph
sort -k1,1 -k2,2n S2DOX221117.avg.coverage.bedGraph > S2DOX221117.sort.avg.coverage.bedGraph
bedGraphToBigWig S2DOX221117.sort.avg.coverage.bedGraph /scratch1/jodie/wolbachia/Micro-C/reference_genomes/Dmel_ref.genome S2DOX221117.avg.coverage.bw

bigWigMerge S2wMel221117-1.coverage.bw S2wMel221117-2.coverage.bw S2wMel221117-3.coverage.bw S2wMel221117-4.coverage.bw S2wMel221117-5.coverage.bw S2wMel221117-6.coverage.bw S2wMel221117.avg.coverage.bedGraph
sort -k1,1 -k2,2n S2wMel221117.avg.coverage.bedGraph > S2wMel221117.sort.avg.coverage.bedGraph
bedGraphToBigWig S2wMel221117.sort.avg.coverage.bedGraph /scratch1/jodie/wolbachia/Micro-C/reference_genomes/Dmel_ref.genome S2wMel221117.avg.coverage.bw

bigWigMerge JW18wMel221117-1.coverage.bw JW18wMel221117-2.coverage.bw JW18wMel221117-3.coverage.bw JW18wMel221117-4.coverage.bw JW18wMel221117-5.coverage.bw JW18wMel221117-6.coverage.bw JW18wMel221117.avg.coverage.bedGraph
sort -k1,1 -k2,2n JW18wMel221117.avg.coverage.bedGraph > JW18wMel221117.sort.avg.coverage.bedGraph
bedGraphToBigWig JW18wMel221117.sort.avg.coverage.bedGraph /scratch1/jodie/wolbachia/Micro-C/reference_genomes/Dmel_ref.genome JW18wMel221117.avg.coverage.bw

bigWigMerge JW18DOX221117-1.coverage.bw JW18DOX221117-2.coverage.bw JW18DOX221117-3.coverage.bw JW18DOX221117-4.coverage.bw JW18DOX221117-5.coverage.bw JW18DOX221117-6.coverage.bw JW18DOX221117.avg.coverage.bedGraph
sort -k1,1 -k2,2n JW18DOX221117.avg.coverage.bedGraph > JW18DOX221117.sort.avg.coverage.bedGraph
bedGraphToBigWig JW18DOX221117.sort.avg.coverage.bedGraph /scratch1/jodie/wolbachia/Micro-C/reference_genomes/Dmel_ref.genome JW18DOX221117.avg.coverage.bw