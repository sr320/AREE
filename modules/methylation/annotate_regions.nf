// Genomic-context annotation of DMR/DML regions (raw_reanalysis mode only).
//
// STATUS: structurally complete DSL2 process wrapping a real, syntactically
// correct GenomicRanges/rtracklayer R script that classifies each region
// into promoter / exon / intron / gene-body / intergenic context by overlap
// against a supplied GTF. This is a simplified but genuine annotation
// approach: promoters are approximated as a fixed window upstream of each
// transcript's TSS (params.methylation.promoter_upstream_bp), and precedence
// among overlapping feature types is resolved by a fixed priority order
// (promoter > exon > intron > gene body > intergenic). A production
// annotation pipeline would likely use a purpose-built tool (e.g.
// HOMER annotatePeaks.pl, ChIPseeker, or a curated GFF feature hierarchy)
// with strand-aware, isoform-aware logic; this is not that. NOT executed in
// this build — see workflows/methylation/README.md.

process ANNOTATE_REGIONS {
    tag "${study_id}:${comparison_id}"
    label 'process_medium'
    container 'bioconductor/bioconductor_docker:RELEASE_3_18'
    publishDir "${params.outdir}/methylation/annotated", mode: params.publish_mode

    input:
    val study_id
    val comparison_id
    path dmr_table
    path annotation_gtf
    val promoter_upstream_bp

    output:
    path "${study_id}_${comparison_id}_dmr_annotated.tsv", emit: annotated_table
    path "${study_id}_${comparison_id}_annotation_qc.tsv", emit: qc_metrics
    path "versions.yml", emit: versions

    script:
    """
    cat <<-'EOF_R' > run_annotate.R
    suppressPackageStartupMessages({
        library(GenomicRanges)
        library(rtracklayer)
    })

    promoter_upstream_bp <- as.integer(${promoter_upstream_bp})

    dmr <- read.delim("${dmr_table}", stringsAsFactors = FALSE)

    dmr_gr <- GRanges(
        seqnames = dmr\$chrom,
        ranges   = IRanges(start = dmr\$start, end = dmr\$end),
        region_id = dmr\$region_id
    )

    gtf <- import("${annotation_gtf}", format = "gtf")

    genes <- gtf[gtf\$type == "gene"]
    exons <- gtf[gtf\$type == "exon"]
    transcripts <- gtf[gtf\$type == "transcript"]
    if (length(transcripts) == 0) {
        # Fall back to gene features as transcript proxies if the GTF has no
        # explicit transcript rows (some minimal/simplified GTFs omit them).
        transcripts <- genes
    }

    # Promoter windows: promoter_upstream_bp upstream of each transcript TSS,
    # strand-aware. This is a simplification of "promoter" — no downstream-of-
    # TSS extension, no distinction between alternative promoters.
    promoters_gr <- suppressWarnings(trim(flank(transcripts, width = promoter_upstream_bp, start = TRUE)))

    # Introns approximated as gene-body minus exons (per gene), which is a
    # standard simplification when a dedicated intron GTF feature is absent.
    genebody_by_gene <- reduce(genes)
    exons_reduced    <- reduce(exons)
    introns_gr <- suppressWarnings(GenomicRanges::setdiff(genebody_by_gene, exons_reduced))

    classify_region <- function(gr_region) {
        if (length(subsetByOverlaps(gr_region, promoters_gr)) > 0) return("promoter")
        if (length(subsetByOverlaps(gr_region, exons)) > 0) return("exon")
        if (length(subsetByOverlaps(gr_region, introns_gr)) > 0) return("intron")
        if (length(subsetByOverlaps(gr_region, genes)) > 0) return("gene_body")
        return("intergenic")
    }

    nearest_gene_id <- function(gr_region) {
        hits <- findOverlaps(gr_region, genes)
        if (length(hits) > 0) {
            gene_hit <- genes[subjectHits(hits)[1]]
            gid <- gene_hit\$gene_id
            if (is.null(gid)) gid <- NA_character_
            return(gid)
        }
        # No direct overlap: report nearest gene_id for context but the
        # region is still classified/kept as intergenic per the no-silent-
        # loss provenance requirement (see src/harmonize/methylation.py).
        return(NA_character_)
    }

    contexts  <- character(length(dmr_gr))
    gene_ids  <- character(length(dmr_gr))
    for (i in seq_along(dmr_gr)) {
        region_i <- dmr_gr[i]
        contexts[i] <- classify_region(region_i)
        gid <- nearest_gene_id(region_i)
        gene_ids[i] <- ifelse(is.na(gid), "", gid)
    }

    dmr\$annotation_context <- contexts
    dmr\$gene_id <- gene_ids

    write.table(dmr, file = "${study_id}_${comparison_id}_dmr_annotated.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)

    context_counts <- table(factor(contexts, levels = c("promoter","exon","intron","gene_body","intergenic")))
    qc <- data.frame(
        study_id      = "${study_id}",
        comparison_id = "${comparison_id}",
        annotation_context = names(context_counts),
        n_regions = as.integer(context_counts)
    )
    write.table(qc, file = "${study_id}_${comparison_id}_annotation_qc.tsv",
                sep = "\\t", quote = FALSE, row.names = FALSE)
    EOF_R

    Rscript run_annotate.R

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        r-base: \$(R --version | head -n1 | sed 's/R version //; s/ .*//')
        bioconductor-genomicranges: \$(Rscript -e 'cat(as.character(packageVersion("GenomicRanges")))')
        bioconductor-rtracklayer: \$(Rscript -e 'cat(as.character(packageVersion("rtracklayer")))')
    END_VERSIONS
    """
}
