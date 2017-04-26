from fusion.models import Gene, Chromosome, CellLine, Fusion
import json
import os
import csv
from sets import Set
from neomodel import db
import glob

# Create your views here.
from django.http import HttpResponse

def fusion_events(request, chromosome):
    return HttpResponse("You asked for fusion events of chromosome " + chromosome)

def genes(request, howmany):
    return HttpResponse("Genes: " + str(Gene.nodes.all()[0:int(howmany)]))

def count_genes(request, prefix):
    return HttpResponse(len(Gene.nodes.filter(symbol__startswith=prefix)))

def print_file(request):
    f = open("test.txt", "r")
    return HttpResponse(f.read())

def chromosomes(request):
    return HttpResponse(json.dumps(open(os.path.dirname(__file__) + "/" + "chromosomes.txt").read().splitlines()))

# def cell_lines(request):
#     response = []
#     for cell_line in CellLine.nodes.all():
#         response.append(cell_line.cell_line)
#         
#     return HttpResponse(json.dumps(response))

def search_indels_by_region(request, chromosome, start, end):
    
    print("Chromosome: " + chromosome + ", start="+start + ", end="+end)
    
    details = []
    response = {}
    response['details'] = {"header": ["id", "fusion point 1", "fusion point 2"], "items": details}
    
    for chrom in Chromosome.nodes.filter(chromosome__exact=chromosome):
        for res in chrom.fusion_events.filter(fusion_point_1__startswith=start):
            details.append([{"value": res.id}, {"value": res.fusion_point_1}, {"value": res.fusion_point_2}])
        
    return HttpResponse(json.dumps(response))

def get_chromosomes_cell_lines(request, cell_line):
    response = []
    
    c = CellLine.nodes.get(cell_line=cell_line)
    print(vars(c.fusion_events))
#     for fusion in c.fusion_events:
#         response.append(str(fusion.fusion_point_1) + "-" + str(fusion.fusion_point_2))
    
    return HttpResponse(json.dumps(response))

def show_info(request, filename):
    response = {}
    details = []
    header = []
    
    line_no = 0
    for line in open(os.path.dirname(__file__) + "/" + filename):
        line_no += 1
        
        fields = line.split("\t")
        if line_no <= 2:
            if line_no == 1:
                header = fields
                print(header)
            else:
                for i in range(0, len(header)):
                    header[i] = header[i] + " " + fields[i]
                    
        else:
            details.append(fields)
            
    response['details'] = {"header": header, "items": details}
    
    return HttpResponse(json.dumps(response))

# def search_for_cell_line(request,c_line):
#     response = {}
#     
#     header = get_header()
#     
#     # recupero fusioni nella linea cellulare
#     fusions = []
#     for fusion in CellLine.nodes.get(cell_line = c_line).happen:
#         fusions.append(fusion)
#     
#     rows = build_rows(fusions)
#     
#     response['rows'] = {"header": header, "items": rows}
#     
#     return HttpResponse(json.dumps(response))

def search_for_cell_line(request, c_line):
    response = {}
    header = get_header()
    #get_cell_line_from_disease("Colon adenocarcinoma")    
    # recupero fusioni nella linea cellulare
    fusions = []
    #se per tutte le linee cellulari mostra tutte le fusioni
    if c_line == "ALL":
        for c_l in CellLine.nodes.all():
            for fusion in c_l.happen:
                fusions.append(fusion)
                #print(fusion)
    else:
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)
    #print(fusions)
    rows = build_rows(fusions)
    #print(rows)
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

# def search_for_chromosome(request,c_line,chromos,start_point,end_point):
#     response = {}
#     rows = []
#     header = get_header()
#     
#     print(request)
#     
#     # recupero fusioni nella linea cellulare
#     fusions = []
#     for fusion in CellLine.nodes.get(cell_line = c_line).happen:
#         print(fusion)
#         if fusion.at_chromosome.match(fusion_point__gte=start_point) and fusion.at_chromosome.match(fusion_point__lte=end_point):
#             if fusion.at_chromosome.filter(id__exact=chromos):
#                 fusions.append(fusion)
#     
#     rows = build_rows(fusions)
# #     print(rows)
#     
#     response['rows'] = {"header": header, "items": rows}
#     return HttpResponse(json.dumps(response))

def search_for_chromosome(request, c_line, chromos):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    #se ho specificato solo il cromosoma, cerco tutte le fusioni in tutte le linee cellulari che coinvolgono il cromosoma
    if c_line == "ALL" and chromos != "ALL":
        c = Chromosome.nodes.get(chromosome = chromos)
        for fusion in c.fromFusiontoChromosome:
            fusions.append(fusion)
    #se ho specificato solo la linea cellulare, cerco tutte le fusioni che avvengono in uqella linea cellulare nel determinato intervallo
    elif c_line != "ALL" and chromos == "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)
    #se ho sia linea cellulare che cromosoma specificati, cerco tutte le fusioni nella linea cellulare che coinvolgono il cromosoma
    elif c_line != "ALL" and chromos != "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.at_chromosome.filter(chromosome__exact=chromos):
                fusions.append(fusion)
    
    rows = build_rows(fusions)
    #print(rows)
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_gene(request,c_line,gene_one,gene_two):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    
    #se ho specificato il gene ma non la linea cellulare, mostro tutte le fusioni che coinvolgono il gene
    if c_line == "ALL" and gene_one != "ALL" and gene_two == "ALL":
        if "ENSG" in gene_one:
            g = Gene.nodes.get(ensid = gene_one)
            if g.had:
                for fusion in g.had:
                    fusions.append(fusion)
            if g.fromFusionToGene:
                for fusion in g.fromFusionToGene:
                    fusions.append(fusion)
        else:
            g = Gene.nodes.get(symbol = gene_one)
            if g.had:
                for fusion in g.had:
                    fusions.append(fusion)
            if g.fromFusionToGene:
                for fusion in g.fromFusionToGene:
                    fusions.append(fusion)
    elif c_line == "ALL" and gene_one != "ALL" and gene_two != "ALL":
        if "ENSG" in gene_one and "ENSG" in gene_two:
            g = Gene.nodes.get(ensid = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].ensid == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(ensid = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].ensid == gene_one:
                        fusions.append(fusion)
        elif "ENSG" in gene_one and "ENSG" not in gene_two:
            g = Gene.nodes.get(ensid = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(symbol = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].ensid == gene_one:
                        fusions.append(fusion)
        elif "ENSG" not in gene_one and "ENSG" in gene_two:
            g = Gene.nodes.get(symbol = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].ensid == gene_two:
                        fusions.append(fusion)
        
            g = Gene.nodes.get(ensid = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_one:
                        fusions.append(fusion)
        else:
            g = Gene.nodes.get(symbol = gene_one)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_two:
                        fusions.append(fusion)
                        
            g = Gene.nodes.get(symbol = gene_two)
            if g.had:
                for fusion in g.had: 
                    if fusion.with_gene.all()[0].symbol == gene_one:
                        fusions.append(fusion)    
    #se ho specificato la linea cellulare ma non il gene, mostro tutte le fusioni per quella linea cellulare (ANALOGO A SEARCH FOR CELL_LINE)
    elif c_line != "ALL" and gene_one == "ALL" and gene_two == "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            fusions.append(fusion)     
    #se ho specificato sia la linea cellulare che il gene, cerco tutte le fusioni che coinvolgono quel gene nella linea cellulare
    elif c_line != "ALL" and gene_one != "ALL" and gene_two == "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.fromGeneToFusion.filter(symbol__exact=gene_one) or fusion.fromGeneToFusion.filter(ensid__exact=gene_one) or fusion.with_gene.filter(symbol__exact=gene_one) or fusion.with_gene.filter(ensid__exact=gene_one): #PROBLEMA
                fusions.append(fusion)
    elif c_line != "ALL" and gene_one != "ALL" and gene_two != "ALL":
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if (fusion.fromGeneToFusion.filter(symbol__exact=gene_one) and fusion.with_gene.filter(symbol__exact=gene_two)) or (fusion.fromGeneToFusion.filter(symbol__exact=gene_two) and fusion.with_gene.filter(symbol__exact=gene_one)) or (fusion.fromGeneToFusion.filter(ensid__exact=gene_one) and fusion.with_gene.filter(ensid__exact=gene_two)) or (fusion.fromGeneToFusion.filter(ensid__exact=gene_two) and fusion.with_gene.filter(ensid__exact=gene_one)): #PROBLEMA
                fusions.append(fusion)
            
    rows = build_rows(fusions)
    #print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_single_gene(request, gene_one, c_line):
    response = {}
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []

    if c_line == "ALL":
        for fusion in Gene.nodes.get(symbol = gene_one).had:
#             print(fusion)
            fusions.append(fusion)
    else:
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if fusion.fromGeneToFusion.filter(symbol__exact=gene_one) or fusion.fromGeneToFusion.filter(gene_id__exact=gene_one) or fusion.with_gene.filter(symbol__exact=gene_one) or fusion.with_gene.filter(gene_id__exact=gene_one):
                fusions.append(fusion)
            
    rows = build_rows(fusions)
#     print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_pair_gene(request, gene_one, gene_two, c_line):
    print("search for pair gene")
    response = {}
    header = get_header()
    
    fusions = []

    # recupero fusioni nella linea cellulare
    if c_line == "ALL":
        for fusion in Gene.nodes.get(symbol = gene_one).had:
            if fusion.with_gene.filter(symbol__exact=gene_two):
#                 print(fusion)
                fusions.append(fusion)
    else:
        for fusion in CellLine.nodes.get(cell_line = c_line).happen:
            if (fusion.fromGeneToFusion.filter(symbol__exact=gene_one) and fusion.with_gene.filter(symbol__exact=gene_two)) or (fusion.fromGeneToFusion.filter(symbol__exact=gene_two) and fusion.with_gene.filter(symbol__exact=gene_one)):
                fusions.append(fusion)

    rows = build_rows(fusions)
#     print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_single_exon(request,c_line,exon_one):
    response = {}
    rows = []
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    for fusion in CellLine.nodes.get(cell_line = c_line).happen:
        if fusion.at_exon.filter(exon__exact=exon_one):
            fusions.append(fusion)

    rows = build_rows(fusions)
#     print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_pair_exon(request,c_line,exon_one,exon_two):
    response = {}
    rows = []
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    for fusion in CellLine.nodes.get(cell_line = c_line).happen:
        if fusion.at_exon.filter(exon__exact=exon_one) and fusion.at_exon.filter(exon__exact=exon_two):
            fusions.append(fusion)

    rows = build_rows(fusions)
#     print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_single_transcript(request,c_line,transcript_one):
    response = {}
    rows = []
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    for fusion in CellLine.nodes.get(cell_line = c_line).happen:
        for couple in fusion.with_trans_couple:
            if couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one) or couple.with_other_transcript.filter(transcript__exact=transcript_one):
                fusions.append(fusion)

    rows = build_rows(fusions)
#     print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_pair_transcript(request,c_line,transcript_one,transcript_two):
    response = {}
    rows = []
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    for fusion in CellLine.nodes.get(cell_line = c_line).happen:
        for couple in fusion.with_trans_couple:
            if (couple.fromTranscriptToCouple.filter(transcript__exact=transcript_one) and couple.with_other_transcript.filter(transcript__exact=transcript_two)) or (couple.fromTranscriptToCouple.filter(transcript__exact=transcript_two) and couple.with_other_transcript.filter(transcript__exact=transcript_one)):
                fusions.append(fusion)

    rows = build_rows(fusions)
#     print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def search_for_fusion_information(request,c_line,algorithm,fusion_description,predicted_effect1,predicted_effect2):
    response = {}
    rows = []
    header = get_header()
    
    # recupero fusioni nella linea cellulare
    fusions = []
    for fusion in CellLine.nodes.get(cell_line = c_line).happen:
        predicted_effect_1 = fusion.fromGeneToFusion.relationship(fusion.fromGeneToFusion.all()[0]).predicted_effect
        predicted_effect_2 = fusion.with_gene.relationship(fusion.with_gene.all()[0]).predicted_effect
        if (algorithm in fusion.fusion_finding_method) and (fusion_description in fusion.description) and (predicted_effect1 == predicted_effect_1) and (predicted_effect2 == predicted_effect_2):
            fusions.append(fusion)

    rows = build_rows(fusions)
#     print(rows)

    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

# def build_rows(fusions):
#     rows = []
#     # ora che ho solo le fusioni interessate recupero le informazioni e mi costruisco la riga
#     for myfusion in fusions:
#         # recupero cell line
#         cellLine = myfusion.fromFusionToCellLine.all()[0].cell_line
#         
#         #recupero dati dai geni
#         gene1 = myfusion.fromGeneToFusion.all()[0]
#         strand_1 = myfusion.fromGeneToFusion.relationship(gene1).strand
#         predicted_effect_1 = myfusion.fromGeneToFusion.relationship(gene1).predicted_effect
#         gene2 = myfusion.with_gene.all()[0]
#         strand_2 = myfusion.with_gene.relationship(gene2).strand
#         predicted_effect_2 = myfusion.with_gene.relationship(gene2).predicted_effect
# 
#         #recupero cromosomi 
#         chromosome1 = []
#         chromosome2 = []
#         fusion_point_1 = ''
#         fusion_point_2 = ''
#         for chrom in myfusion.at_chromosome:
#             if chrom.of_gene.filter(symbol__exact=gene1.symbol):
#                 fusion_point_1 = myfusion.at_chromosome.relationship(chrom).fusion_point
#                 chromosome1 = chrom
#             if chrom.of_gene.filter(symbol__exact=gene2.symbol):
#                 fusion_point_2 = myfusion.at_chromosome.relationship(chrom).fusion_point
#                 chromosome2 = chrom
# 
#         #recupero esoni
#         exon1 = []
#         exon2 = []
# #         fusion_partner_1 = ''
# #         fusion_partner_2 = ''
#         for exon in myfusion.at_exon:
#             if exon.in_gene.filter(symbol__exact=gene1.symbol):
# #                 fusion_partner_1 = myfusion.at_exon.relationship(exon).fusion_partner
#                 exon1 = exon
#             if exon.in_gene.filter(symbol__exact=gene2.symbol):
# #                 fusion_partner_2 = myfusion.at_exon.relationship(exon).fusion_partner
#                 exon2 = exon
#                 
#         #recupero trascritti e proteine
#         transcript_couples = []
#         proteins = []
#         for couple in myfusion.with_trans_couple:
#             transcript1 = couple.fromTranscriptToCouple.all()[0]
# #             print(transcript1)
#             transcript1_position = couple.fromTranscriptToCouple.relationship(transcript1).position
#             transcript2 = couple.with_other_transcript.all()[0]
# #             print(transcript2)
#             transcript2_position = couple.with_other_transcript.relationship(transcript2).position
#             
#             transcript_couples.append(transcript1.transcript+":"+str(transcript1_position)+" - "+transcript2.transcript+":"+str(transcript2_position))
#             proteins.append(couple.with_protein.all()[0].protein)
#             
#         #costruisco la riga
#         row = []
#         row.append(cellLine)
#         row.append(gene1.symbol+" - "+gene2.symbol)
#         row.append(str(gene1.gene_id)+" - "+str(gene2.gene_id))
#         if exon1 or exon2:
#             row.append(str(exon1.exon)+" - "+str(exon2.exon))
#         else:
#             row.append("No exons")
#         row.append(str([str(chromosome1.chromosome)+":"+str(fusion_point_1)+":"+str(strand_1), str(chromosome2.chromosome)+":"+str(fusion_point_2)+":"+str(strand_2)]))
#         row.append(myfusion.description)
#         row.append(myfusion.common_mapping_reads)
#         row.append(myfusion.spanning_pairs)
#         row.append(myfusion.spanning_unique_reads)
#         row.append(myfusion.longest_anchor_found)
#         row.append(myfusion.fusion_finding_method)
#         row.append(myfusion.fusion_sequence)
#         if(predicted_effect_1!=predicted_effect_2):
#             row.append(predicted_effect_1+"/"+predicted_effect_2)
#         else:
#             row.append(predicted_effect_1)
#         row.append(transcript_couples)
#         row.append(proteins)
# 
#         rows.append(row)
#     return rows

def cell_lines(request):
    response = []

    response.append("ALL")
    
    for cell_line in CellLine.nodes.all():
        response.append(cell_line.cell_line)
        
    return HttpResponse(json.dumps(response))

def statistics_all(request):
    response = {}

    header = ['Fusion events', 'Transcripts', 'Genes', 'Predicted proteins', 'Exons', 'Viruses']
    rows = []
    response['details'] = {"header": header, "items": rows}
    fusion_events = 90319 # len(Fusion.nodes.all())
    transcripts = 25626 # len(Transcript.nodes.all())
    genes = 15601 # len(Gene.nodes.all())
    protein = 36828 # len(Protein.nodes.all())
    exon = 12976 # len(Exon.nodes.all())
    virus = 2793 # len(Virus.nodes.all())
    
    rows.append([fusion_events, transcripts, genes, protein, exon, virus])

    return HttpResponse(json.dumps(response))

def statistics_by_chromosome(request, chrom):
    response = {}
 
    header = ['Fusion events', 'Transcripts', 'Genes', 'Predicted proteins', 'Exons', 'Viruses']
    rows = []
    response['details'] = {"header": header, "items": rows}
    
    fusion_events = []
    transcripts = []
    genes = []
    proteins = []
    exons = []
    viruses = []
    
    node = Chromosome.nodes.get(chromosome = chrom)
    
    for fusion in node.fromFusiontoChromosome:
        
        fusion_events.append(fusion)
        
        for exon in fusion.at_exon:
            exons.append(exon)
            
        for gene in fusion.with_gene:
            genes.append(gene)
        
        for couple in fusion.with_trans_couple:
            for protein in couple.with_protein:
                proteins.append(protein)
            for transcript in couple.fromTranscriptToCouple:
                transcripts.append(transcript)
        
        for cell_line in fusion.fromFusionToCellLine:
            for virus in cell_line.with_viruses:
                viruses.append(virus)
                
    rows.append([len(fusion_events), len(transcripts), len(Set(genes)), len(proteins), len(exons), len(viruses)])
 
    return HttpResponse(json.dumps(response))

def sort_chromosomes(c):
    try:
        return int(c)
    except ValueError:
        return c
    
def fusion_by_chromosome(request):
    response = {}
    
    header = []
    rows = []
    
    response['details'] = {"header": header, "items": rows}
    
    chromosomes = Chromosome.nodes.all()
    sorted_chromosomes = sorted(chromosomes, key=lambda c: sort_chromosomes(c.chromosome))
    
    for chrom in sorted_chromosomes:
        
        header.append("chr " + chrom.chromosome)
        rows.append(len(chrom.fromFusiontoChromosome))
        
    return HttpResponse(json.dumps(response))

def download_data(request):
    
    response = {}

    header = ['Cell line', 'Fusion gene', 'Viruses', 'Summary']
    rows = []
    response['details'] = {"header": header, "items": rows}
    
    sizes = ['B', 'KB', 'MB', 'GB', 'TB']
    for cell in CellLine.nodes.all():
        cell_line_id = cell.cell_line
        
#         size_bytes = os.path.getsize(url)
#         size_converted = size_bytes
#         size_human_index = 0
#         
#         while size_converted/1024 > 1:
#             size_converted = size_converted / 1024
#             size_human_index += 1
#         
#         size_human = size_converted + " " + sizes[size_human_index]
        
        rows.append([
            cell_line_id,
            {"label": "candidate fusion genes",
             "url": "downloads/cells/" + cell_line_id + ".txt"},
            {"label": "viruses bacteria phages",
             "url": "downloads/viruses/" + cell_line_id + ".txt"},
            {"label": "summary",
             "url": "downloads/summary/" + cell_line_id + ".txt"}])

    return HttpResponse(json.dumps(response))

# def generate_statistics(request):
#     
#     infos = get_ccle_infos()
#     distribution = {}
#     for row in infos["items"]:
#         disease = row[2]
#         
#         if disease not in distribution:
#             distribution[disease] = 0
#             
#         distribution[disease] += 1
#         
#     filename = os.path.dirname(__file__) + "/statistics/Disease_CellLine.csv"
#     file = open(filename,'w') 
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Disease", "CellLine"])
#     for disease in distribution:
#         writer.writerow([disease, distribution[disease]])
#     file.close()
#     
#     return HttpResponse()

# def get_header():
#     return ["Cell line",
#         "Gene pair symbols",
#         "Gene pair EnsIDs",
#         "Exon pair",
#         "Chromosome : fusion point : strand",
#         "Description",
#         "Counts of common mapping reads",
#         "Spanning pairs",
#         "Spanning unique reads",
#         "Longest anchor found",
#         "Fusion finding method",
#         "Fusion sequence",
#         "Predicted effect",
#         "Predicted fused transcripts",
#         "Predicted fused proteins"]

def get_ccle_infos():
    header = ["ID","Cell Line","Disease","Disease name"]
    rows = []
    
    txt_file = open(os.path.dirname(__file__) + "/ccle_ids.txt", "r")
    next(txt_file)
    for line in txt_file:
        print(line)
        line = line.rstrip()
        words = line.split("\t")
        rows.append([words[0].rstrip(),words[1],words[2],words[3]])
        #print(words)
    
    response = {"header": header, "items": rows}

#     print(response)
    
    return response

    #prendo in input una stringa che <E8> il nome della malattia, mi ricavo le linee cellulari corrispondenti e mi ricavo la tabella relativa
def get_cell_line_from_disease(disease):
    ccle_infos = get_ccle_infos()["items"]
    
    cls = []
    for row in ccle_infos:
        if disease in row:
            cls.append(row[0])
    
    return cls

def search_for_disease(request, disease):
    cls = get_cell_line_from_disease(disease)
    fusions = []
    
    response = {}
    header = get_header()
    
    for cl in cls:
        
        for fusion in CellLine.nodes.get(cell_line = cl).happen:
            fusions.append(fusion)
        
    rows = build_rows(fusions)
    
    response['rows'] = {"header": header, "items": rows}    
        
    return HttpResponse(json.dumps(response))

def get_distribution(request, node1, node2, howmany, sorting):
    
    if not howmany: howmany = -1
    howmany = int(howmany)
    
    response = {}
    labels = []
    header = []
    items = []
    
    header_translation = get_ccle_infos()
    
    lines = open(os.path.dirname(__file__) + "/statistics/" + node1 + "_" + node2 + ".csv")
    line_no = 0
    for line in lines:
        line_no += 1
        line = line.rstrip()
        fields = line.split(",")
        if line_no == 1:
            labels = fields
            continue
        
        header.append(fields[0])
        items.append(fields[1])
    lines.close()
    
    if sorting == "DESC" or sorting == "ASC":
        header, items = zip(*sorted(zip(header, items), key=lambda pair: int(pair[1])))
        if sorting == "DESC":
            header = list(reversed(header))
            items = list(reversed(items))
    else:
        header, items = zip(*sorted(zip(header, items), key=lambda pair: (not pair[0].isdigit(), pair[0].zfill(3))))

    if howmany >= 0 and len(items) > howmany:
        items = items[:howmany]
        header = header[:howmany]
    
    header = list(header)
    
    if node1 == "CellLine":
        for i,item in enumerate(header):
            for ccle in header_translation["items"]:
                if ccle[0] == item:
                    header[i] = ccle[1]
                    break
    
    response['details'] = {"labels": labels, "header": header, "items": items}
    
    return HttpResponse(json.dumps(response))

def get_single_distribution(request, label, value):
    
    response = {}
    header = []
    items = []
    
    files = glob.glob(os.path.dirname(__file__) + "/statistics/" + label + "_*.csv")
    
    for filename in files:
        lines = open(filename)

        line_no = 0
        
        for line in lines:
            line_no += 1
            line = line.rstrip()
            fields = line.split(",")
            
            if line_no == 1:
                header_fields = fields
                continue
            
            if fields[0] == value:
                header.append(header_fields[1])
                items.append(fields[1])
                break
            
        lines.close()
    
#     header, items = zip(*sorted(zip(header, items), key=lambda pair: (not pair[0].isdigit(), pair[0].zfill(3))))
    
    response['details'] = {"header": header, "items": items}
    
    return HttpResponse(json.dumps(response))



def build_rows(fusions):
    rows = []
    
    ccle_infos = get_ccle_infos()
    for myfusion in fusions:
        # recupero dati cell line
        cellLine = myfusion.fromCellLineToFusion.all()[0].cell_line
        disease = ""
        acronym = ""
        
        for row in ccle_infos["items"]:
            if cellLine in row:
                disease = row[3]
                acronym = row[2]
                
        #recupero dati dai geni
        gene1 = myfusion.fromGeneToFusion.all()[0]
        gene2 = myfusion.with_gene.all()[0]
       
      
        #check
        fc_flag = "NO FC"
        es_flag = "NO ES"
        th_flag = "NO TH"
        
        if(myfusion.with_fc_script):
            fc_flag = '{ type: "button", action: "dialog", url: fusioncatcher/'+str(myfusion.fusion_id)+' }'
        if(myfusion.with_eric_script):
            es_flag = '{ type: "button", action: "dialog", url: ericscript/'+str(myfusion.fusion_id)+' }'
        if(myfusion.with_tophat_script):
            th_flag = '{ type: "button", action: "dialog", url: tophat/'+str(myfusion.fusion_id)+' }'
        
        rows.append([disease,acronym,cellLine,gene1.symbol,gene2.symbol,fc_flag,es_flag,th_flag])
        
        #print(gene1.symbol+" "+gene2.symbol)
        #print(fc_fusions)   
        #print(es_fusions)
        #print("\n\n\n")
        
    return rows


#BUILD_FC_TABLE E BUILD_ES_TABLE per costruire la tabella dei dettagli per ogni coppia. Quando premo il pulsante relativo all'algoritmo deve stamparmi tutte le fusioni relative

def build_fc_table(request,fus_id):
    
    response = {}
    header = get_fc_header()
    
    fusions = []
    
#     for fcFusion in Fusion.nodes.get(fusion_id = fus_id).with_fc_script:
#         fusions.append(fcFusion)
        
    for fusion in Fusion.nodes.all():
        if fusion.fusion_id == int(fus_id):
            if fusion.with_fc_script:
                fusions.append(fusion)
        
        
    rows = build_fc_rows(fusions)
    
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def build_es_table(request,fus_id):
    
    response = {}
    header = get_es_header()
    
    fusions = []
    
#     for esFusion in Fusion.nodes.get(fusion_id = fus_id).with_eric_script:
#         fusions.append(esFusion)

    for fusion in Fusion.nodes.all():
        if fusion.fusion_id == int(fus_id):
            if fusion.with_eric_script:
                fusions.append(fusion)
        
    rows = build_es_rows(fusions)
    
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def build_th_table(request,fus_id):
    
    response = {}
    header = get_tophat_header()
    
    fusions = []
    
#     for esFusion in Fusion.nodes.get(fusion_id = fus_id).with_eric_script:
#         fusions.append(esFusion)

    for fusion in Fusion.nodes.all():
        if fusion.fusion_id == int(fus_id):
            if fusion.with_tophat_script:
                fusions.append(fusion)
        
    rows = build_tophat_rows(fusions)
    
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def build_fc_rows(fusions):
    
    rows = []
    # ora che ho solo le fusioni FC interessate recupero le informazioni e mi costruisco la riga
    fc_fusions = []
    
    for fus in fusions:
        for fcfusion in fus.with_fc_script:
            fc_fusions.append(fcfusion)
    for myfusion in fc_fusions:
        gene1 = myfusion.fromFusionToFusionCatcher.all()[0].fromGeneToFusion.all()[0]
        gene2 = myfusion.fromFusionToFusionCatcher.all()[0].with_gene.all()[0]
        
        strand_1 = myfusion.strand_1
        predicted_effect_1 = myfusion.predicted_effect_1
        strand_2 = myfusion.strand_2
        predicted_effect_2 = myfusion.predicted_effect_2
        
        #recupero cromosomi 
        chromosome1 = gene1.fromChromosomeToGene.all()[0]
        chromosome2 = gene2.fromChromosomeToGene.all()[0]
        fusion_point_1 = myfusion.fusion_point_1
        fusion_point_2 = myfusion.fusion_point_2

        #recupero esoni        
        exon1 = []
        exon2 = []

        for exon in myfusion.at_exon:
            if exon.in_gene.filter(symbol__exact=gene1.symbol):
                exon1 = exon
            if exon.in_gene.filter(symbol__exact=gene2.symbol):
                exon2 = exon
        #recupero trascritti e proteine
        transcript_couples = []
        proteins = []
        for couple in myfusion.with_trans_couple:
            transcript1 = couple.fromTranscriptToCouple.all()[0]
            transcript1_position = couple.fromTranscriptToCouple.relationship(transcript1).position
            transcript2 = couple.with_other_transcript.all()[0]
            transcript2_position = couple.with_other_transcript.relationship(transcript2).position
            
            transcript_couples.append(transcript1.transcript+":"+str(transcript1_position)+" - "+transcript2.transcript+":"+str(transcript2_position))
            if couple.with_protein.all():
                proteins.append(couple.with_protein.all()[0].protein)
            
        #costruisco la riga
        row = []
        row.append(gene1.symbol+" - "+gene2.symbol)
        if exon1 or exon2:
            row.append(exon1.exon+" - "+exon2.exon)
        else:
            row.append("No exons")
        row.append([chromosome1.chromosome+":"+str(fusion_point_1)+":"+strand_1, chromosome2.chromosome+":"+str(fusion_point_2)+":"+strand_2])
        row.append(myfusion.description)
        row.append(myfusion.common_mapping_reads)
        row.append(myfusion.spanning_pairs)
        row.append(myfusion.spanning_unique_reads)
        row.append(myfusion.longest_anchor_found)
        row.append(myfusion.fusion_finding_method)
        row.append(myfusion.fusion_sequence)
        if(predicted_effect_1!=predicted_effect_2):
            row.append(predicted_effect_1+"/"+predicted_effect_2)
        else:
            row.append(predicted_effect_1)
        row.append(transcript_couples)
        row.append(proteins)   
        #print(row)
        rows.append(row)
    return rows

#SCRIVERE BUILD_ES_ROWS
def build_es_rows(fusions):
    rows = []
    #recupero dati degli eventi di fusione
        
    es_fusions = []
    
    for fus in fusions:
        for esfusion in fus.with_eric_script:
            es_fusions.append(esfusion)
    
    for myfusion in es_fusions:
        gene1 = myfusion.fromFusionToEricScript.all()[0].fromGeneToFusion.all()[0]
        gene2 = myfusion.fromFusionToEricScript.all()[0].with_gene.all()[0]
        
        #recupero cromosomi 
        chromosome1 = gene1.fromChromosomeToGene.all()[0]
        chromosome2 = gene2.fromChromosomeToGene.all()[0]    
        breakpoint1 = myfusion.breakpoint_1
        breakpoint2 = myfusion.breakpoint_2
        strand1 = myfusion.strand_1
        strand2 = myfusion.strand_2
        
        crossing_reads = myfusion.crossing_reads
        spanning_reads = myfusion.spanning_reads
        mean_intersize = myfusion.mean_intersize
        homology = myfusion.homology
        fusion_type = myfusion.fusion_type
        junction_sequence = myfusion.junction_sequence
        gene_expr_1 = myfusion.gene_expr_1
        gene_expr_2 = myfusion.gene_expr_2
        gene_expr_fused = myfusion.gene_expr_fused
        es = myfusion.es
        gjs = myfusion.gjs
        us = myfusion.us
        eric_score = myfusion.eric_score
        
        #costruisco la riga
        row = []
        row.append(gene1.symbol+" - "+gene2.symbol)
        row.append([chromosome1.chromosome+":"+str(breakpoint1)+":"+strand1, chromosome2.chromosome+":"+str(breakpoint2)+":"+strand2])
        row.append(crossing_reads)
        row.append(spanning_reads)
        row.append(mean_intersize)
        row.append(homology)
        row.append(fusion_type)
        row.append(junction_sequence)
        row.append(gene_expr_1)
        row.append(gene_expr_2)
        row.append(gene_expr_fused)
        row.append(es)
        row.append(gjs)
        row.append(us)
        row.append(eric_score)
        
        rows.append(row)

    return rows

def build_tophat_rows(fusions):
    rows = []
    #recupero dati degli eventi di fusione
        
    th_fusions = []
    
    for fus in fusions:
        for thfusion in fus.with_tophatscript:
            th_fusions.append(thfusion)
    
    for myfusion in th_fusions:
        gene1 = myfusion.fromFusionToTophatScript.all()[0].fromGeneToFusion.all()[0]
        gene2 = myfusion.fromFusionToTophatScript.all()[0].with_gene.all()[0]
        
        #recupero cromosomi 
        chromosome1 = gene1.fromChromosomeToGene.all()[0]
        chromosome2 = gene2.fromChromosomeToGene.all()[0]    
        left_coord = myfusion.left_coord
        right_coord = myfusion.right_coord
        
        spanning_reads = myfusion.spanning_reads
        spanning_mate_pairs = myfusion.spanning_mate_pairs
        spanning_mate_pairs_end = myfusion.spanning_mate_pairs_end
        
        #costruisco la riga
        row = []
        row.append(gene1.symbol+" - "+gene2.symbol)
        row.append([chromosome1.chromosome+":"+str(left_coord), chromosome2.chromosome+":"+str(right_coord)])
        row.append(spanning_reads)
        row.append(spanning_mate_pairs)
        row.append(spanning_mate_pairs_end)

        
        rows.append(row)

    return rows
    

def search_viruses(request,c_line,vir):
    response = {}
    rows = []
    header = ["Cell line",
        "Virus/bacteria/phage name",
        "GI",
        "NC"]
    
    #recupero virus nella linea cellulare
    #se ho specificato solo il virus, cerca per tutte le linee cellulari
    if c_line == "ALL" and vir != "":
        for c_l in CellLine.nodes.all():
            for virus in c_l.with_viruses:
                if vir in virus.name or vir == virus.gi or vir == virus.ref:
                    rows.append([c_l.cell_line, virus.name, virus.gi, virus.ref])
    #se ho specificato solo la linea cellulare, cerca tutti i virus per la linea cellulare
    elif c_line !="" and vir == "trolo":
        for virus in CellLine.nodes.get(cell_line = c_line).with_viruses:
            rows.append([c_line, virus.name, virus.gi, virus.ref])
    #se ho specificato sia virus che linea cellulare, cerca quel virus per quella linea cellulare
    elif c_line != "" and vir != "":
        for virus in CellLine.nodes.get(cell_line = c_line).with_viruses:
            if vir in virus.name or vir == virus.gi or vir == virus.ref:
                rows.append([c_line, virus.name, virus.gi, virus.ref])
                
    response['rows'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))
        
def generate_statistics(request):
    cellLine_gene = {(('CellLine','cell_line'),('Gene','symbol')):[[CellLine,'HAPPEN',Fusion,'HAD',Gene],[CellLine,'HAPPEN',Fusion,'WITH',Gene]]}
    cellLine_chromosome = {(('CellLine','cell_line'),('Chromosome','chromosome')):[[CellLine,'HAPPEN',Fusion,'AT_CHROMOSOME',Chromosome]]}
    cellLine_fusion = {(('CellLine','cell_line'),('Fusion','fusion_id')):[[CellLine,'HAPPEN',Fusion]]}
    chromosome_cellLine = {(('Chromosome','chromosome'),('CellLine','cell_line')):[[Chromosome,'AT_CHROMOSOME',Fusion,'HAPPEN',CellLine]]}
    chromosome_fusion = {(('Chromosome','chromosome'),('Fusion','fusion_id')):[[Chromosome,'AT_CHROMOSOME',Fusion]]}
    chromosome_gene = {(('Chromosome','chromosome'),('Gene','symbol')):[[Chromosome,'OF_GENE',Gene]]}
    gene_cellLine = {(('Gene','symbol'),('CellLine','cell_line')):[[Gene,'HAD',Fusion,'HAPPEN',CellLine],[Gene,'WITH',Fusion,'HAPPEN',CellLine]]}
    gene_chromosome = {(('Gene','symbol'),('Chromosome','chromosome')):[[Gene,'AT_CHROMOSOME',Chromosome]]}
    gene_fusion = {(('Gene','symbol'),('Fusion','fusion_id')):[[Gene,'HAD',Fusion],[Gene,'WITH',Fusion]]}
    
    stats = [cellLine_gene, cellLine_chromosome, cellLine_fusion, chromosome_cellLine,chromosome_fusion, chromosome_gene, gene_cellLine, gene_chromosome, gene_fusion]
    for pair in stats:
        gen_stat_file(pair)
    #gen_stat_file(gene_cellLine)
    
    return HttpResponse()    
    
def gen_stat_file(pairs):
    ids = list(pairs.keys())[0]
    paths = list(pairs.values())
    paths = paths[0]
    #print(ids)
    #print(paths)
    startnode = paths[0][0]
    #print(startnode)
    file =  open(ids[0][0]+'_'+ids[1][0]+'.csv','w')
    writer = csv.writer(file, lineterminator='\n')
    writer.writerow([ids[0][0],ids[1][0]])
    for x in startnode.nodes.all():
        print(x)
        nodes = {}
        for path in paths:
            #print(path)
            #print(len(path))
            query = "match (" + path[0].__name__ + "{"+ids[0][1]+":'"+eval("x."+ids[0][1])+"'})-"
            for i in range(1,len(path)-2,2):
                query = query + "[:" +  path[i] + "]-(" + path[i+1].__name__ + ")-" 
                #print(path[i+1])
            query = query + "[:" + path[len(path)-2] + "]-(a:" + path[(len(path)-1)].__name__+") return distinct a"
            #print(query)
            if(db.cypher_query(query)[0]):
                for row in db.cypher_query(query)[0]:
                    #print(row[0].properties[eval("'"+ids[1][1]+"'")])
                    if row[0].properties[eval("'"+ids[1][1]+"'")] not in nodes:
                        #print(row[1].properties[eval("'"+ids[1][1]+"'")])
                        nodes[row[0].properties[eval("'"+ids[1][1]+"'")]] = row[0].properties[eval("'"+ids[1][1]+"'")]
        #print(len(nodes))
        writer.writerow([eval("x."+ids[0][1]),len(nodes)])
    file.close()
    
def get_header():
    return ["Cancer",
            "Acronym",
            "CCLE",
            "Gene 1",
            "Gene 2",
            "Fusion Catcher",
            "Ericscript",
            "Tophat"]
    
def get_fc_header():
    return ["Gene pair symbols",
        "Exon pair",
        "Chromosome : fusion point : strand",
        "Description",
        "Counts of common mapping reads",
        "Spanning pairs",
        "Spanning unique reads",
        "Longest anchor found",
        "Fusion finding method",
        "Fusion sequence",
        "Predicted effect",
        "Predicted fused transcripts",
        "Predicted fused proteins"]
    
def get_es_header():
    return ["Gene pair symbols",
        "Chromosome : breakpoint : strand",
        "Crossing reads",
        "Spanning reads",
        "Mean intersize",
        "Homology",
        "Fusion type",
        "Junction sequence",
        "Gene expr 1",
        "Gene expr 2",
        "Gene expr fused",
        "Es",
        "Gjs",
        "Us",
        "Eric score"]
    
def get_tophat_header():
    return ["Gene pair symbols",
            "Chromosome : coordinate",
            "Spanning reads",
            "Spanning mate pairs",
            "Spanning mate pairs where one end spans a fusion"]        
