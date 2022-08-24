.libPaths('/home/umcg-agulyaeva/SOFTWARE/R_LIB')
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
sessionInfo()


d <- DiagrammeR::grViz("
    digraph {
        graph []
        node [width = 2, fontname = Arial]
            N1 [label = '58,776 virus-like contigs\\nfrom the 4 Dutch cohort', shape = box]
            N2 [label = '41,284 genomes\\nfrom the databases', shape = box]
            e1 [label = 'Dereplication', shape = plaintext, height = 0.1]
            N3 [label = '30,461 vOTUs', shape = box]
            e2 [label = 'Read mapping', shape = plaintext, height = 0.1]
            N4 [label = '15,196 detected vOTUs', shape = box]
            e3 [label = <<I>Caudoviricetes</I>  genomes with terminal repeats>, shape = plaintext, height = 0.1]
            N5 [label = '1,899 selected vOTUs', shape = box]
            e4 [label = 'Detected in > 5% samples of a Dutch cohort', shape = plaintext, height = 0.1]
            N6 [label = '54 selected vOTUs', shape = box]
        edge []
            N1->e1 [arrowhead = none]
            N2->e1 [arrowhead = none]
            e1->N3 [arrowhead = vee]
            N3->e2 [arrowhead = none]
            e2->N4 [arrowhead = vee]
            N4->e3 [arrowhead = none]
            e3->N5 [arrowhead = vee]
            N5->e4 [arrowhead = none]
            e4->N6 [arrowhead = vee]
        { rank = same; N1; N2 }
    }

")


d <- DiagrammeRsvg::export_svg(d)
d <- charToRaw(d)
rsvg::rsvg_pdf(d, file = 'workflow_diagram.pdf')
