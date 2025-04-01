#Average bigwig

bigWigMerge JW18-DOX-1.coverage.bw JW18-DOX-2.coverage.bw JW18-DOX.avg.coverage.bedGraph
sort -k1,1 -k2,2n  JW18-DOX.avg.coverage.bedGraph >  JW18-DOX.sort.avg.coverage.bedGraph
bedGraphToBigWig  JW18-DOX.sort.avg.coverage.bedGraph /scratch1/jodie/wolbachia/Micro-C/reference_genomes/Dmel_ref.genome  JW18-DOX.avg.coverage.bw

bigWigMerge JW18-wMel-1.coverage.bw JW18-wMel-2.coverage.bw JW18-wMel.avg.coverage.bedGraph
sort -k1,1 -k2,2n  JW18-wMel.avg.coverage.bedGraph >  JW18-wMel.sort.avg.coverage.bedGraph
bedGraphToBigWig  JW18-wMel.sort.avg.coverage.bedGraph /scratch1/jodie/wolbachia/Micro-C/reference_genomes/Dmel_ref.genome  JW18-wMel.avg.coverage.bw

bigWigMerge JW18-wRi-1.coverage.bw JW18-wRi-2.coverage.bw JW18-wRi.avg.coverage.bedGraph
sort -k1,1 -k2,2n  JW18-wRi.avg.coverage.bedGraph >  JW18-wRi.sort.avg.coverage.bedGraph
bedGraphToBigWig  JW18-wRi.sort.avg.coverage.bedGraph /scratch1/jodie/wolbachia/Micro-C/reference_genomes/Dmel_ref.genome  JW18-wRi.avg.coverage.bw

bigWigMerge JW18-wWil-1.coverage.bw JW18-wWil-2.coverage.bw JW18-wWil.avg.coverage.bedGraph
sort -k1,1 -k2,2n  JW18-wWil.avg.coverage.bedGraph >  JW18-wWil.sort.avg.coverage.bedGraph
bedGraphToBigWig  JW18-wWil.sort.avg.coverage.bedGraph /scratch1/jodie/wolbachia/Micro-C/reference_genomes/Dmel_ref.genome  JW18-wWil.avg.coverage.bw