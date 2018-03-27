from neomodel import (StructuredNode, StructuredRel, StringProperty, IntegerProperty, ArrayProperty, BooleanProperty)
from neomodel.relationship_manager import RelationshipTo, RelationshipFrom
from neomodel.properties import FloatProperty

class IN_COUPLE(StructuredRel):
    position = IntegerProperty()
    
class WITH_OTHER_TRANSCRIPT(StructuredRel):
    position = IntegerProperty()   
    
class WITH_VIRUSES(StructuredRel):
    count_of_mapping_reads = IntegerProperty()

#NODES
class CellLine(StructuredNode):
    cell_line = StringProperty()
    name = StringProperty()
    disease = StringProperty()
    disease_name = StringProperty()

    with_viruses = RelationshipTo('Virus',"WITH_VIRUSES",model=WITH_VIRUSES)
    happen = RelationshipFrom('FusionCell', "FUSIONCELL_TO_CELLINE")

class GeneCouple(StructuredNode):
    gene_couple_id = StringProperty()
    gene1 = StringProperty()
    gene2 = StringProperty()
    
    cosmic = StringProperty()
    
    fusion_cells = RelationshipTo('FusionCell', "COUPLE_TO_FUSIONCELL")
    gene_five_prime = RelationshipFrom('Gene', "HAD")
    gene_three_prime = RelationshipTo('Gene', "WITH")

class FusionCell(StructuredNode):
    fusion_cell_id = StringProperty()
    num_algo = IntegerProperty()
    
    gene1 = StringProperty()
    gene2 = StringProperty()
    cell_line = StringProperty()
    
    quality = StringProperty()
    
    algos = ArrayProperty()
    annotations = ArrayProperty()
    
    gene_couple = RelationshipFrom('GeneCouple', "COUPLE_TO_FUSIONCELL")
    cell_lines = RelationshipTo('CellLine', "FUSIONCELL_TO_CELLINE")
    fusion_points = RelationshipTo('FusionPoint', "FUSIONCELL_TO_FUSIONPOINT")
    
class FusionPoint(StructuredNode):
    fusion_point_id = StringProperty()
    chromosome1 = StringProperty()
    chromosome2 = StringProperty()
    breakpoint1 = StringProperty()
    breakpoint2 = StringProperty()
    
    num_algo = IntegerProperty()
    fusion_cells = RelationshipFrom('FusionCell', "FUSIONCELL_TO_FUSIONPOINT")
#     algorithms = RelationshipTo('Algorithm', "FUSIONPOINT_TO_ALGORITHM")
    
    fusioncatcher_events = RelationshipTo('FusionCatcher', "WITH_FC_SCRIPT")
    ericscript_events = RelationshipTo('EricScript', "WITH_ERIC_SCRIPT")
    tophat_events = RelationshipTo('Tophat', "WITH_TOPHAT_SCRIPT")
    jaffa_events = RelationshipTo('Jaffa', "WITH_JAFFA_SCRIPT")
    
class Algorithm(StructuredNode):
    algorithm_id = StringProperty()
    algorithm = StringProperty()
    
    fusion_points = RelationshipFrom('FusionPoint', "FUSIONPOINT_TO_ALGORITHM")
    fusioncatcher_events = RelationshipTo('FusionCatcher', "WITH_FC_SCRIPT")
    ericscript_events = RelationshipTo('EricScript', "WITH_ERIC_SCRIPT")
    tophat_events = RelationshipTo('Tophat', "WITH_TOPHAT_SCRIPT")

class Chromosome(StructuredNode):
    chromosome = StringProperty()

    of_gene = RelationshipTo('Gene',"OF_GENE")
    fromFusiontoChromosome = RelationshipFrom('FusionPoint', "AT_CHROMOSOME")
    fromChromosomeToFusionPoint = RelationshipTo('FusionPoint', "FROM_CHROMOSOME")
    
class Couple(StructuredNode):
    couple = IntegerProperty() 

    with_other_transcript = RelationshipTo('Transcript',"WITH_OTHER_TRANSCRIPT", model=WITH_OTHER_TRANSCRIPT)
    with_protein = RelationshipTo('Protein',"WITH_PROTEIN")
    fromTranscriptToCouple = RelationshipFrom('Transcript',"IN_COUPLE", model=IN_COUPLE)
    fusioncatcher_events = RelationshipFrom('FusionCatcher',"WITH_TRANS_COUPLE")
    
class Fusion(StructuredNode):
    fusion_id = StringProperty()
    num_supporting_algorithms = IntegerProperty()

    with_fc_script = RelationshipTo('FusionCatcher',"WITH_FC_SCRIPT")
    with_eric_script =  RelationshipTo('EricScript',"WITH_ERIC_SCRIPT")
    with_tophat_script = RelationshipTo('Tophat',"WITH_TOPHAT_SCRIPT")
    with_gene = RelationshipTo('Gene',"WITH")
    at_chromosome = RelationshipTo('Chromosome',"AT_CHROMOSOME")
    fromCellLineToFusion = RelationshipFrom('CellLine',"HAPPEN")
    fromGeneToFusion = RelationshipFrom('Gene',"HAD")
    
class FusionCatcher(StructuredNode):
    fusioncatcher_id = StringProperty()
    
    quality = StringProperty()
    description = ArrayProperty()
    common_mapping_reads = IntegerProperty()
    spanning_pairs = IntegerProperty()
    spanning_unique_reads = IntegerProperty()
    longest_anchor_found = IntegerProperty()
    fusion_finding_method = ArrayProperty()
    fusion_sequence = StringProperty()
    fusion_point_1 = IntegerProperty()
    fusion_point_2 = IntegerProperty()
    strand_1 = StringProperty()
    strand_2 = StringProperty()
    predicted_effect_1 = StringProperty()
    predicted_effect_2 = StringProperty()
    cell_line = StringProperty()

    from_exon = RelationshipFrom('Exon', "FROM_EXON")
    at_exon = RelationshipTo('Exon', "AT_EXON")
    
    with_trans_couple = RelationshipTo('Couple', "WITH_TRANS_COUPLE")
    with_gene = RelationshipTo('Gene',"WITH")
    
    fusion_point = RelationshipFrom('FusionPoint', "WITH_FC_SCRIPT")
    
class EricScript(StructuredNode):
    ericscript_id = StringProperty()
    breakpoint_1 = StringProperty()
    strand_1 = StringProperty()
    breakpoint_2 = StringProperty()
    strand_2 = StringProperty()
    crossing_reads = IntegerProperty()
    spanning_reads = IntegerProperty()
    mean_intersize = FloatProperty()
    homology = StringProperty()
    fusion_type = StringProperty()
    junction_sequence = StringProperty()
    gene_expr_1 = FloatProperty()
    gene_expr_2 = FloatProperty()
    gene_expr_fused = FloatProperty()
    es = FloatProperty()
    gjs = StringProperty()
    us = FloatProperty()
    eric_score = FloatProperty()
    cell_line = StringProperty()
    
#     algorithm = RelationshipFrom('Algorithm', "WITH_ERIC_SCRIPT")
    fusion_point = RelationshipFrom('FusionPoint', "WITH_ERIC_SCRIPT")

class Tophat(StructuredNode):
    tophat_id = IntegerProperty()
    spanning_reads = StringProperty()
    spanning_mate_pairs = StringProperty()
    spanning_mate_pairs_end = StringProperty()
    cell_line = StringProperty()
    
#     algorithm = RelationshipFrom('Algorithm', "WITH_TOPHAT_SCRIPT")
    fusion_point = RelationshipFrom('FusionPoint', "WITH_TOPHAT_SCRIPT")

class Jaffa(StructuredNode):
    jaffa_id = IntegerProperty()
    transcript = StringProperty()
    spanning_pairs = IntegerProperty()
    spanning_reads = IntegerProperty()
    chrom1 = StringProperty()
    base1 = IntegerProperty()
    strand1 = StringProperty()
    chrom2 = StringProperty()
    base2 = IntegerProperty()
    strand2 = StringProperty()
    inframe = BooleanProperty()
    fusion_genes = StringProperty()
    known = BooleanProperty()
    classification = StringProperty()
    sequence = StringProperty()
    cell_line = StringProperty()
    
    fusion_point = RelationshipFrom('FusionPoint', "WITH_JAFFA_SCRIPT")

class Gene(StructuredNode):
    ensid = StringProperty()
    symbol = StringProperty()
    description = StringProperty()
    
    census = BooleanProperty()
    hallmark = BooleanProperty()

    fromExonToGene = RelationshipFrom('Exon',"IN_GENE")
    fromChromosomeToGene = RelationshipFrom('Chromosome',"OF_GENE")
    
    had = RelationshipTo('GeneCouple', "HAD")
    fromFusionToGene = RelationshipFrom('GeneCouple','WITH')
    
class Protein(StructuredNode):
    protein = StringProperty()

    fromCoupleToProtein = RelationshipFrom('Couple',"WITH_PROTEIN")
    
class Transcript(StructuredNode):
    transcript = StringProperty()
    
    to_couple = RelationshipTo('Couple',"IN_COUPLE", model=IN_COUPLE)
    from_couple = RelationshipFrom('Couple',"WITH_OTHER_TRANSCRIPT", model=WITH_OTHER_TRANSCRIPT)

class Exon(StructuredNode):
    exon = StringProperty()

    in_gene = RelationshipTo('Gene',"IN_GENE")
    fromExonToFusion = RelationshipTo('FusionCatcher', "FROM_EXON")
    fromFusionToExon = RelationshipFrom('FusionCatcher', "AT_EXON")

class Virus(StructuredNode):
    name = StringProperty()
    gi = StringProperty()
    ref = StringProperty()
    
    fromCellLineToVirus = RelationshipFrom('CellLine',"WITH_VIRUSES", model=WITH_VIRUSES)
    
    