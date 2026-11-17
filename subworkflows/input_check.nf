//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../modules/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv

    main:
    SAMPLESHEET_CHECK ( samplesheet )
        .csv
        .splitCsv ( header:true, sep:',' )
        .map { create_fastq_channel(it) }
        .set { reads }

    emit:
    reads                                     // channel: [ val(meta), [ reads ] ]
}

// A relative FASTQ path is looked up from the launch directory, as before, and
// failing that next to the samplesheet, so a samplesheet can travel with its
// data (tests/samplesheet_test.csv).
def fastq_file(String path) {
    def f = file(path)
    if (f.exists() || new File(path).isAbsolute()) return f
    def beside = file(params.input).parent.resolve(path)
    return beside.exists() ? beside : f
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def create_fastq_channel(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id         = row.sample
    meta.single_end = row.single_end.toBoolean()

    // add path(s) of the fastq file(s) to the meta map
    def fastq_meta = []
    def fastq_1 = fastq_file(row.fastq_1)
    if (!fastq_1.exists()) {
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${row.fastq_1}"
    }
    if (meta.single_end) {
        fastq_meta = [ meta, [ fastq_1 ] ]
    } else {
        def fastq_2 = fastq_file(row.fastq_2)
        if (!fastq_2.exists()) {
            exit 1, "ERROR: Please check input samplesheet -> Read 2 FastQ file does not exist!\n${row.fastq_2}"
        }
        fastq_meta = [ meta, [ fastq_1, fastq_2 ] ]
    }
    return fastq_meta
}
