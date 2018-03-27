from fusion.models import Gene, Chromosome, CellLine, FusionCell, Protein, FusionPoint, Virus, Algorithm, GeneCouple, Exon, Transcript, FusionCatcher, EricScript
import json
import os
import csv
from sets import Set
from neomodel import db
import glob
import cgi
import time
import copy
import math
import re
from django.views.decorators.csrf import ensure_csrf_cookie
import datetime

from collections import Counter

# Create your views here.
from django.http import HttpResponse

MAX_RESULTS = 50

red_list = ["1000genome", "banned", "bodymap2", "cacg", "conjoing", "cta_gene", "ctb_gene", "ctc_gene", "ctd_gene", "gap<1K", "10K<gap<100K", "100K<gap<200K", "duplicates", "ensembl_fully_overlapping", "ensembl_same_strand_overlapping", "gtex", "hpa", "mt", "non_cancer_tissues", "non_tumor_cells", "pair_pseudo_genes", "paralogs", "readthrough", "refseq_fully_overlapping", "refseq_same_strand_overlapping", "rp11_gene", "rp_gene", "rrna", "short_distance", "similar_reads", "similar_symbols,ucsc_fully_overlapping", "ucsc_same_strand_overlapping", "gencode_fully_overlapping", "gencode_partially_overlapping", "gencode_same_strand_overlapping"]
orange_list = ["adjacent", "ambiguous", "ensembl_partially_overlapping", "fragments", "healthy", "refseq_partially_overlapping", "ucsc_partially_overlapping"]
green_list = ["18cancers", "cell_lines", "cgp", "chimerdb2", "chimerdb3kb", "chimerdb3pub", "chimerdb3seq", "cosmic", "gliomas", "oncogene", "prostates", "tcga", "ticdb", "known"]

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
    
    response = []
    for chromosome in ["ALL"] + open(os.path.dirname(__file__) + "/" + "chromosomes.txt").read().splitlines():
        response.append({"id": chromosome, "label": chromosome, "img": "chromosome-icon.png"})
    
    return HttpResponse(json.dumps(response))

def chromosomes_simple(request):
    
    response = []
    for chromosome in open(os.path.dirname(__file__) + "/" + "chromosomes.txt").read().splitlines():
        response.append({"id": chromosome, "label": chromosome, "img": "chromosome-icon.png"})
    
    return HttpResponse(json.dumps(response))

def exons(request):
    
    response = []
    
    response.append({"id": "ALL", "label": "Include any exon in results"})
    
    for exon in Exon.nodes:
        response.append({"id": exon.exon, "label": exon.exon})
    
    return HttpResponse(json.dumps(response))

def search_indels_by_region(request, chromosome, start, end):
    
    print("Chromosome: " + chromosome + ", start="+start + ", end="+end)
    
    details = []
    response = {}
    response['details'] = {"header": ["id", "fusion point 1", "fusion point 2"], "items": details}
    
    for chrom in Chromosome.nodes.filter(chromosome__iexact=chromosome):
        for res in chrom.fusion_events.filter(fusion_point_1__startswith=start):
            details.append([{"value": res.id}, {"value": res.fusion_point_1}, {"value": res.fusion_point_2}])
        
    return HttpResponse(json.dumps(response))

#def get_chromosomes_cell_lines(request, cell_line):
    #response = []
    
    #c = CellLine.nodes.get(cell_line=cell_line)
       
    #return HttpResponse(json.dumps(response))

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

def filter_table(response, c_line):
    if c_line == "ALL": return response
    else:
        simplified_header = response["header"][1:]
        simplified_rows = []
        for row in response["items"]:
            print(row)
            simplified_rows.append(row[1:])
            
        response["header"] = simplified_header
        response["items"] = simplified_rows
                
        return response

def search_for_cell_line_with_chromosome(request, c_line):
    response = {}
    header = ["Gene1", "Gene1 start", "Gene1 chromosome", "Gene2", "Gene2 start", "Gene2 chromosome"]

    rows = []
    
    info = get_gene_infos()
    
    #se per tutte le linee cellulari mostra tutte le fusioni
    if c_line != "ALL":
        for fusion in CellLine.nodes.get(name = c_line).happen.filter(num_algo__gte=2):
            event = {}
            
            genes = [fusion.gene1.upper(), fusion.gene2.upper()]
            if genes[0] not in info: continue
            if genes[1] not in info: continue
            
            single_info1 = [genes[0]] + info[genes[0]]
            single_info2 = [genes[1]] + info[genes[1]]
            
            event["g1"] = genes[0] 
            event["g1start"] = single_info1[2]
            event["g1chr"] = single_info1[1]
            
            event["g2"] = genes[1] 
            event["g2start"] = single_info2[2]
            event["g2chr"] = single_info2[1]
            
            rows.append(event)

    #print(rows)
    response = {"header": header, "items": rows}
    
    #if c_line != "ALL": response = filter_table(response, c_line)
    
    return HttpResponse(json.dumps(response))

def search_for_cell_line_events(request, c_line):
    response = {}
    header = ["5' gene", "3' gene", "5' chrom.", "3' chrom."]

    rows = []
    
    info = get_gene_infos()
    
    #se per tutte le linee cellulari mostra tutte le fusioni
    if c_line != "ALL":
        for fusion in CellLine.nodes.get(name = c_line).happen.filter(num_algo__gte=2):
            event = []
            
            genes = [fusion.gene1.upper(), fusion.gene2.upper()]
            if genes[0] not in info: continue
            if genes[1] not in info: continue
            
            single_info1 = [genes[0]] + info[genes[0]]
            single_info2 = [genes[1]] + info[genes[1]]
            
            event = [genes[0], genes[1], single_info1[1], single_info2[1]]

            rows.append(event)

    #print(rows)
    response['details'] = {"header": header, "items": rows}
    
    #if c_line != "ALL": response = filter_table(response, c_line)
    
    return HttpResponse(json.dumps(response))

def search_for_chromosome2(request):
    
    return HttpResponse(json.dumps(search_for_chromosome2_private(request)))

def search_for_chromosome2_private(request):
    
    global annotations
    
    data = json.loads(request.body)    
    c_line = data["cell_line"]
    chromos1 = data["chrom1"]
    chromos2 = data["chrom2"]
    num_algorithms = data["num_algorithms"] if "num_algorithms" in data else 1
    novel = data["novel"] if "novel" in data else False
    
    print(data)
    
    total = 0
    
    fusions = []
    
    if c_line == "ALL":
        
        if chromos1 != "ALL" and chromos2 == "ALL":
            c = Chromosome.nodes.get(chromosome = chromos1)
            genes = [gene.symbol for gene in c.of_gene]
            fusions = FusionCell.nodes.filter(gene1__in = genes).filter(num_algo__gte=num_algorithms)
            
            if novel:
                filtered_fusions = []
                for fusion in fusions:
                    if (not fusion.annotations or len(set(fusion.annotations).intersection(annotations.keys())) == 0) and not fusion.gene_couple[0].cosmic:
                        if (novel == "strict" and fusion.quality == "N/A") or (novel == "relaxed" and (fusion.quality == "N/A" or fusion.quality == "grey")):
                            filtered_fusions.append(fusion)
                fusions = filtered_fusions
                    
            total = len(fusions)
            
        elif chromos1 == "ALL" and chromos2 != "ALL" :
            c = Chromosome.nodes.get(chromosome = chromos2)
            genes = [gene.symbol for gene in c.of_gene]
            fusions = FusionCell.nodes.filter(gene2__in = genes).filter(num_algo__gte=num_algorithms)
            
            if novel:
                filtered_fusions = []
                for fusion in fusions:
                    if (not fusion.annotations or len(set(fusion.annotations).intersection(annotations.keys())) == 0) and not fusion.gene_couple[0].cosmic:
                        if (novel == "strict" and fusion.quality == "N/A") or (novel == "relaxed" and (fusion.quality == "N/A" or fusion.quality == "grey")):
                            filtered_fusions.append(fusion)
                fusions = filtered_fusions
            
            total = len(fusions)
            
        elif chromos1 != "ALL" and chromos2 != "ALL":
            
            # Non cosmic + non in database e non verde/giallo/rosso
            if novel == "strict":
                novel_condition = "NOT EXISTS(g.cosmic) AND f.quality = 'N/A' AND (NOT EXISTS(f.annotations) OR (" + " AND ".join(["NOT '"+x+"' IN f.annotations" for x in annotations.keys()]) + ")) AND"
                novel_path = "--(g:GeneCouple)"
            elif novel == "relaxed":
                novel_condition = "NOT EXISTS(g.cosmic) AND (f.quality = 'N/A' OR f.quality = 'grey') AND (NOT EXISTS(f.annotations) OR (" + " AND ".join(["NOT '"+x+"' IN f.annotations" for x in annotations.keys()]) + ")) AND"
                novel_path = "--(g:GeneCouple)"

            r, m = db.cypher_query("MATCH (p:FusionPoint)--(f:FusionCell)"+novel_path+" WHERE "+novel_condition+" f.num_algo >= "+str(num_algorithms)+" AND p.chromosome1 = '"+chromos1+"' AND p.chromosome2 = '"+chromos2+"' return count(DISTINCT f)")
            total = r[0][0]
             
            results, meta = db.cypher_query("MATCH (p:FusionPoint)--(f:FusionCell)"+novel_path+" WHERE "+novel_condition+" f.num_algo >= "+str(num_algorithms)+" AND p.chromosome1 = '"+chromos1+"' AND p.chromosome2 = '"+chromos2+"' return DISTINCT f ")
            fusions = [FusionCell.inflate(row[0]) for row in results]

#             for fusion_point in FusionPoint.nodes.filter(chromosome1__exact = chromos1).filter(chromosome2__exact = chromos2):
#                 for fusion in fusion_point.fusion_cells.filter(num_algo__gte=num_algorithms):
#                     fusions.append(fusion)
    else:
        if chromos1 == "ALL" and chromos2 == "ALL":
            for fusion in CellLine.nodes.get(cell_line = c_line).happen.filter(num_algo__gte=num_algorithms):
                fusions.append(fusion)
                
        elif chromos1 != "ALL" and chromos2 == "ALL":
            for fusion in CellLine.nodes.get(cell_line = c_line).happen.filter(num_algo__gte=num_algorithms):
                if len(fusion.fusion_points.filter(chromosome1__exact=chromos1)) > 0:
                    fusions.append(fusion)
                
        elif chromos1 == "ALL" and chromos2 != "ALL":
            for fusion in CellLine.nodes.get(cell_line = c_line).happen.filter(num_algo__gte=num_algorithms):
                if len(fusion.fusion_points.filter(chromosome2__exact=chromos2)) > 0:
                    fusions.append(fusion)
                    
        elif chromos1 != "ALL" and chromos2 != "ALL":
            for fusion in CellLine.nodes.get(cell_line = c_line).happen.filter(num_algo__gte=num_algorithms):
                if len(fusion.fusion_points.filter(chromosome1__exact=chromos1).filter(chromosome2__exact=chromos2)) > 0:
                    fusions.append(fusion)
        
        if novel:
            filtered_fusions = []
            for fusion in fusions:
                if (not fusion.annotations or len(set(fusion.annotations).intersection(annotations.keys())) == 0) and not fusion.gene_couple[0].cosmic:
                    if (novel == "strict" and fusion.quality == "N/A") or (novel == "relaxed" and (fusion.quality == "N/A" or fusion.quality == "grey")):
                        filtered_fusions.append(fusion)
            fusions = filtered_fusions
        
        total = len(fusions)
    
    fusions = sorted(fusions, key=lambda fusion: fusion.num_algo, reverse=True)
    if "offset" in data:
        fusions = fusions[data["offset"]:data["offset"]+data["limit"]]
    
    rows = build_rows2(fusions, c_line)
    
    response = {"structure": {"field_list": get_header2()}, "total": total, "hits": rows}
    
    return response

def search_for_gene2(request):
    
    return HttpResponse(json.dumps(search_for_gene2_private(request)))

def search_for_gene2_private(request):
    global annotations
    data = json.loads(request.body)
    
    c_line = data["cell_line"]
    gene_one = data["gene_1"]
    gene_two = data["gene_2"]
    num_algorithms = data["num_algorithms"] if "num_algorithms" in data else 1
    novel = data["novel"] if "novel" in data else False
    
    print(data)
    
    fusions = []
    
    gene1 = None
    if gene_one != "ALL":
        gene1 = Gene.nodes.get(ensid = gene_one) if "ENSG" in gene_one else Gene.nodes.get(symbol = gene_one)
    gene2 = None
    if gene_two != "ALL":
        gene2 = Gene.nodes.get(ensid = gene_two) if "ENSG" in gene_two else Gene.nodes.get(symbol = gene_two)
    
    print(gene1, gene2)
    
    if c_line == "ALL":
        fusions = FusionCell.nodes.filter(num_algo__gte=num_algorithms)
        print("FOUND", len(fusions))
        if gene1 is not None: fusions = fusions.filter(gene1__exact=gene1.symbol)
        if gene2 is not None: fusions = fusions.filter(gene2__exact=gene2.symbol)
    else:
        fusions = CellLine.nodes.get(cell_line = c_line).happen
        if gene1 is not None: fusions = fusions.filter(gene1__exact=gene1.symbol)
        if gene2 is not None: fusions = fusions.filter(gene2__exact=gene2.symbol)
        fusions = fusions.filter(num_algo__gte=num_algorithms)
    
    if novel == "strict":
        filtered_fusions = []
        for fusion in fusions:
            if fusion.quality == "N/A" and (not fusion.annotations or not set(fusion.annotations).intersection(annotations.keys())):
                filtered_fusions.append(fusion)
        fusions = filtered_fusions
    elif novel == "relaxed":
        filtered_fusions = []
        for fusion in fusions:
            if (fusion.quality == "grey" or fusion.quality == "N/A") and (not fusion.annotations or not set(fusion.annotations).intersection(annotations.keys())):
                filtered_fusions.append(fusion)
        fusions = filtered_fusions
    
    total = len(fusions)
    fusions = sorted(fusions, key=lambda fusion: fusion.num_algo, reverse=True)
    
    if "offset" in data:
        fusions = fusions[data["offset"]:data["offset"]+data["limit"]]
    
    rows = build_rows2(fusions, c_line)

    response = {"structure": {"field_list": get_header2()}, "total": total, "hits": rows}
    
    return response

def search_for_exon2(request):
    
    return HttpResponse(json.dumps(search_for_exon2_private(request)))

def search_for_exon2_private(request):
    response = {}
    
    data = json.loads(request.body)
    
    c_line = data["cell_line"]
    exon_one = data["exon_1"]
    exon_two = data["exon_2"]
    novel = data["novel"] if "novel" in data else False
    
    if c_line != "ALL" and exon_one == "ALL" and exon_two == "ALL":
        return search_for_cell_line2_private(request)
    
    print(data)
    
    if not data:
        data = {"offset": 0, "limit": 50}

    field_list = get_header2()
    
    fusions_dict = {}
    
    if exon_one != "ALL" and exon_two == "ALL":
        e = Exon.nodes.get(exon = exon_one)
        gene = e.in_gene[0]
        
        for fcfusion in e.fromExonToFusion:
            for algorithm in fcfusion.algorithm:
                for fusion_point in algorithm.fusion_points:
                    for fusion in fusion_point.fusion_cells:
                        if fusion.gene1 == gene.symbol:
                            if c_line == "ALL" or fusion.cell_line == c_line:
                                fusions_dict[fusion.fusion_cell_id] = fusion
        
    elif exon_one == "ALL" and exon_two != "ALL":
        e = Exon.nodes.get(exon = exon_two)
        gene = e.in_gene[0]
        for fcfusion in e.fromFusionToExon:
            for algorithm in fcfusion.algorithm:
                for fusion_point in algorithm.fusion_points:
                    for fusion in fusion_point.fusion_cells:
                        if fusion.gene2 == gene.symbol:
                            if c_line == "ALL" or fusion.cell_line == c_line:
                                fusions_dict[fusion.fusion_cell_id] = fusion
                
    elif exon_one != "ALL" and exon_two != "ALL":
        e = Exon.nodes.get(exon = exon_one)
        gene1 = e.in_gene[0]
        gene2 = Exon.nodes.get(exon = exon_two).in_gene[0]
        
        for fcfusion in e.fromExonToFusion:
                
            if not fcfusion.at_exon.filter(exon__iexact=exon_two): continue
            
            for algorithm in fcfusion.algorithm:
                for fusion_point in algorithm.fusion_points:
                    for fusion in fusion_point.fusion_cells:
                        if fusion.gene1 == gene1.symbol and fusion.gene2 == gene2.symbol:
                            if c_line == "ALL" or fusion.cell_line == c_line:
                                fusions_dict[fusion.fusion_cell_id] = fusion
    
    fusions = fusions_dict.values()
    
    if novel:
        filtered_fusions = []
        for fusion in fusions:
            if (not fusion.annotations or len(set(fusion.annotations).intersection(annotations.keys())) == 0) and not fusion.gene_couple[0].cosmic:
                if (novel == "strict" and fusion.quality == "N/A") or (novel == "relaxed" and (fusion.quality == "N/A" or fusion.quality == "grey")):
                    filtered_fusions.append(fusion)
        fusions = filtered_fusions
    
    total = len(fusions)
    fusions = sorted(fusions, key=lambda fusion: fusion.num_algo, reverse=True)
    
    if "offset" in data:
        fusions = fusions[data["offset"]:data["offset"]+data["limit"]]
        
    rows = build_rows2(fusions, c_line)

    response = {"structure": {"field_list": field_list}, "total": total, "hits": rows}
    
    return response

def search_for_transcript2(request):
    
    return HttpResponse(json.dumps(search_for_transcript2_private(request)))
    
    
def search_for_transcript2_private(request):
    
    data = json.loads(request.body)
    
    c_line = data["cell_line"]
    transcript_one = data["transcript_1"]
    transcript_two = data["transcript_2"]
    novel = data["novel"] if "novel" in data else False

    print(data)
    
    if transcript_one == "ALL" and transcript_two == "ALL":
        return search_for_cell_line2_private(request)

    response = {}
    
    if not data:
        data = {"offset": 0, "limit": 50}

    fusions_dict = {}
    
    if transcript_one != "ALL" and transcript_two == "ALL":
        t1 = None
        try:
            t1 = Transcript.nodes.get(transcript = transcript_one)
        except Transcript.DoesNotExist:
            pass
        
        if t1 is not None:            
            for couple in t1.to_couple:
                for fcfusion in couple.fusioncatcher_events:
                    for fusion_point in algorithm.fusion_point:
                        for fusion in fusion_point.fusion_cells:
                            if c_line == "ALL" or fusion.cell_line == c_line:
                                fusions_dict[fusion.fusion_cell_id] = fusion
                                    
    elif transcript_one == "ALL" and transcript_two != "ALL":
        
        t2 = None
        try:
            t2 = Transcript.nodes.get(transcript = transcript_two)
        except Transcript.DoesNotExist:
            pass
        
        if t2 is not None:   
            for couple in t2.from_couple:
                for fcfusion in couple.fusioncatcher_events:
                    for fusion_point in algorithm.fusion_point:
                        for fusion in fusion_point.fusion_cells:
                            if c_line == "ALL" or fusion.cell_line == c_line:
                                fusions_dict[fusion.fusion_cell_id] = fusion
    
    elif transcript_one != "ALL" and transcript_two != "ALL":
        t1 = None
        try:
            t1 = Transcript.nodes.get(transcript = transcript_one)
        except Transcript.DoesNotExist:
            pass
        
        if t1 is not None:   
            for couple in t1.to_couple:
                for fcfusion in couple.fusioncatcher_events:
                    for fusion_point in fcfusion.fusion_point:
                        for fusion in fusion_point.fusion_cells:
                            if c_line == "ALL" or fusion.cell_line == c_line:
                                fusions_dict[fusion.fusion_cell_id] = fusion
                                
    fusions = fusions_dict.values()
    
    if novel == "strict":
        filtered_fusions = []
        for fusion in fusions:
            #if (not fusion.annotations or len(set(fusion.annotations).intersection(annotations.keys())) == 0) and not fusion.gene_couple[0].cosmic:
            if fusion.quality == "N/A" or fusion.quality == "grey":
                filtered_fusions.append(fusion)
        fusions = filtered_fusions
    
    total = len(fusions)
    fusions = sorted(fusions, key=lambda fusion: fusion.num_algo, reverse=True)
    
    if "offset" in data:
        fusions = fusions[data["offset"]:data["offset"]+data["limit"]]
               
    rows = build_rows2(fusions, c_line)

    response = {"structure": {"field_list": get_header2()}, "total": total, "hits": rows}
    
    return response

import time
def cell_lines(request):
    response = []

    for cell_line in CellLine.nodes.order_by("cell_line"):
        response.append({"id": cell_line.cell_line, "label": cell_line.name, "extra": cell_line.disease_name, "img": "cell-icon.png"})
    
    response = sorted(response, key=lambda cell: cell["label"]);
    
    response.insert(0, {"id": "ALL", "label": "Include any cell line"})
        
    return HttpResponse(json.dumps(response))

def cell_lines_simple(request, prefix = ""):
    response = []

    for cell_line in CellLine.nodes.filter(name__icontains=prefix).order_by("cell_line")[0:50]:
        response.append({"id": cell_line.cell_line, "label": cell_line.name, "extra": cell_line.disease_name, "img": "cell-icon.png"})
        
    response = sorted(response, key=lambda cell: cell["label"]);
    
    if len(response) > 1:
        response.insert(0, {"id": "ALL", "label": "Include any cell line"})
        
    return HttpResponse(json.dumps(response))

def get_genes(request):
    response = []

    response.append({"id": "ALL", "label": "Include any gene in results"})
    
    three_prime_genes = set([line.rstrip() for line in open(os.path.dirname(__file__) + "/data/genes3'.txt").readlines()])
    five_prime_genes = set([line.rstrip() for line in open(os.path.dirname(__file__) + "/data/genes5'.txt").readlines()])
        
    for gene in Gene.nodes:
        response.append({"id": gene.symbol, "label": gene.symbol, "img": "gene-icon.png", "three_prime": gene.symbol in three_prime_genes, "five_prime": gene.symbol in five_prime_genes})
        
    return HttpResponse(json.dumps(response))

def get_genes_three_prime(request, prefix = None):
    response = []

    three_prime_genes = set([line.rstrip() for line in open(os.path.dirname(__file__) + "/data/genes3'.txt").readlines()])
    
    for gene in Gene.nodes:
        disabled = gene.symbol not in three_prime_genes
        item = {"id": gene.symbol, "label": gene.symbol, "img": "gene-icon.png", "disabled": disabled}
        if disabled:
            item["extra"] = "No fusion events are known for this gene"
            
        if len(response) > 50:
            break
            
        if prefix is None or prefix.lower() in gene.symbol.lower():
            response.append(item)
        
    response = sorted(response, key=lambda gene: gene["label"])
    
    if len(response) > 1:
        response.insert(0, {"id": "ALL", "label": "Include any gene in results"})
        
    return HttpResponse(json.dumps(response))

def get_genes_five_prime(request, prefix = None):
    response = []

    five_prime_genes = set([line.rstrip() for line in open(os.path.dirname(__file__) + "/data/genes5'.txt").readlines()])
        
    for gene in Gene.nodes:
        disabled = gene.symbol not in five_prime_genes
        item = {"id": gene.symbol, "label": gene.symbol, "img": "gene-icon.png", "disabled": disabled}
        if disabled:
            item["extra"] = "No fusion events are known for this gene"
            
        if len(response) > 50:
            break
        
        if prefix is None or prefix.lower() in gene.symbol.lower(): 
            response.append(item)
        
    response = sorted(response, key=lambda gene: gene["label"])
    
    if len(response) > 1:
        response.insert(0, {"id": "ALL", "label": "Include any gene in results"})  
      
    return HttpResponse(json.dumps(response))

def transcripts(request):
    
    response = []
    
    response.append({"id": "ALL", "label": "Include any transcript in results"})
    
    for transcript in Transcript.nodes:
        response.append({"id": transcript.transcript, "label": transcript.transcript, "img": "transcript-icon.png"})
    
    return HttpResponse(json.dumps(response))

def get_transcripts_three_prime(request, prefix = None):
    
    filename = os.path.dirname(__file__) + "/" + "/data/transcripts3'.txt"
    if not os.path.exists(filename):
        print("Precalculating all the three prime transcripts")
        writer = open(filename, "w")
        for transcript in Transcript.nodes:
            if len(transcript.from_couple) > 0:
                writer.write(transcript.transcript + "\n")
        writer.close()
        print("Finished precalculating all the three prime transcripts")
    
    response = []

    response.append({"id": "ALL", "label": "Include any transcript in results"})
    
    three_prime_transcripts = set([line.rstrip() for line in open(filename).readlines()])
        
    for transcript in three_prime_transcripts:
        if len(response) > 50: break
        if prefix is None or prefix in transcript: 
            response.append({"id": transcript, "label": transcript, "img": "transcript-icon.png"})
        
    return HttpResponse(json.dumps(response))

def get_transcripts_five_prime(request, prefix = None):
    
    filename = os.path.dirname(__file__) + "/" + "/data/transcripts5'.txt"
    if not os.path.exists(filename):
        print("Precalculating all the five prime transcripts")
        writer = open(filename, "w")
        for transcript in Transcript.nodes:
            if len(transcript.to_couple) > 0:
                writer.write(transcript.transcript + "\n")
        writer.close()
        print("Finished precalculating all the five prime transcripts")
    
    response = []

    response.append({"id": "ALL", "label": "Include any transcript in results"})
    
    five_prime_transcripts = set([line.rstrip() for line in open(filename).readlines()])
        
    for transcript in five_prime_transcripts:
        if len(response) > 50: break
        if prefix is None or prefix in transcript:
            response.append({"id": transcript, "label": transcript, "img": "transcript-icon.png"})
        
    return HttpResponse(json.dumps(response))

def get_exons_three_prime(request, prefix = None):
    
    filename = os.path.dirname(__file__) + "/" + "/data/exons3'.txt"
    if not os.path.exists(filename):
        print("Precalculating all the three prime exons")
        writer = open(filename, "w")
        for exon in Exon.nodes:
            if len(exon.fromFusionToExon) > 0:
                writer.write(exon.exon + "\n")
        writer.close()
        print("Finished precalculating all the three prime exons")
    
    response = []

    three_prime_exons = set([line.rstrip() for line in open(filename).readlines()])
        
    for exon in three_prime_exons:
        if len(response) > 50: break
        if prefix is None or prefix in transcript: 
            response.append({"id": exon, "label": exon, "img": "exon-icon.png"})
    
    response.insert(0, {"id": "ALL", "label": "Include any exon in results"})
    
    return HttpResponse(json.dumps(response))

def get_exons_five_prime(request, prefix = None):
    
    filename = os.path.dirname(__file__) + "/" + "/data/exons5'.txt"
    if not os.path.exists(filename):
        print("Precalculating all the five prime exons")
        writer = open(filename, "w")
        for exon in Exon.nodes:
            if len(exon.fromExonToFusion) > 0:
                writer.write(exon.exon + "\n")
        writer.close()
        print("Finished precalculating all the five prime exons")
    
    response = []

    five_prime_exons = set([line.rstrip() for line in open(filename).readlines()])
        
    for exon in five_prime_exons:
        if len(response) > 50: break
        if prefix is None or prefix in transcript: 
            response.append({"id": exon, "label": exon, "img": "exon-icon.png"})

    response.insert(0, {"id": "ALL", "label": "Include any exon in results"})
        
    return HttpResponse(json.dumps(response))

def viruses(request):
    response = []

#     virus_info = {}
#     for line in open(os.path.dirname(__file__)+ "/statistics/Virus_FusionCell.tsv", "r"):
#         fields = line.rstrip().split("\t")
#         virus_info[fields[0]] = fields[1]
    
    for virus in Virus.nodes:
        
        total = len(virus.fromCellLineToVirus)
        # total = int(virus_info[virus.name])
        item = {
            "id": virus.name,
            "label": virus.name,
            "img": "virus-icon.png",
            "cell_lines": total,
            "extra": str(total) + " cell line" + ("s" if total > 1 else "")
            }
        
        response.append(item)
    
    response.sort(key=lambda x: x["label"])
    response.insert(0, {"id": "ALL", "label": "Include any virus in results"})
    
    return HttpResponse(json.dumps(response))

# def precompute_transcripts(request):
# 
#     writer = open(os.path.dirname(__file__)+ "/statistics/Transcript_FusionCell.tsv", "w")
#     writer.write("Transcript\tFusion events\n")
# 
#     for virus in Transcript.nodes.all():
#         
#         r, m = db.cypher_query("MATCH (f:FusionCell)--(c:CellLine)--(v:Virus) WHERE v.name='"+virus.name+"' return count(distinct(f))")
#         total = r[0][0]
#         
#         writer.write(virus.name +"\t"+ str(total) + "\n")
#     
#     writer.close()
#     return HttpResponse()

def statistics_all(request):
    response = {}

    header = ['# of fusion events', '# of involved gene couples', '# of involved transcripts', '# of involved genes', '# of predicted proteins', '# of involved exons', '# of detected viruses']
    rows = []
    response['details'] = {"header": header, "items": rows}
    fusion_events = len(FusionCell.nodes)
    gene_couples = len(GeneCouple.nodes)
    transcripts = len(Transcript.nodes)
    genes = len(Gene.nodes)
    protein = len(Protein.nodes)
    exon = len(Exon.nodes)
    virus = len(Virus.nodes)
    
#     fusion_events = 990627 # len(Fusion.nodes.all())
#     transcripts = 25626 # len(Transcript.nodes.all())
#     genes = 20260 # len(Gene.nodes.all())
#     protein = 36828 # len(Protein.nodes.all())
#     exon = 12663 # len(Exon.nodes.all())
#     virus = 465 # len(Virus.nodes.all())
    
    rows.append([fusion_events, gene_couples, transcripts, genes, protein, exon, virus])

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

def get_human_size(url):
    size_bytes = os.path.getsize(url)
    size_converted = size_bytes
    size_human_index = 0
    
    while float(size_converted)/1024 > 1:
        size_converted = size_converted / 1024
        size_human_index += 1
    
    sizes = ['B', 'KB', 'MB', 'GB', 'TB']
    size_human = str(size_converted) + " " + sizes[size_human_index]
    
    return size_human
    
# def download_data(request):
#     
#     response = {}
# 
#     header = [
#         'Cell line',
#         'Disease',
#         'FusionCatcher (FC) predicted events',
#         'EricScript (ES) predicted events',
#         'Tophat-Fusion (TF) predicted events',
#         'Viruses information (only by FusionCatcher)',
#         'Summary information (only by FusionCatcher)',
#         'Oncofuse - FusionCatcher (FC)',
#         'Oncofuse - EricScript (ES)',
#         'Oncofuse - Tophat-Fusion (TF)',
#         ]
#     rows = []
#     response['details'] = {"header": header, "items": rows}
#     
#     ccle_map = {}
#     ccle_infos = get_ccle_infos()
#     for cell_info in ccle_infos:
#         ccle_map[cell_info["ID"]] = cell_info
#     
#     for cell in CellLine.nodes.order_by('name'):
#         cell_line_id = cell.cell_line
#         cell_line_name = ccle_map[cell_line_id]["Cell Line"]
#         cell_line_disease = ccle_map[cell_line_id]["Disease name"]
# 
#         download_type = "download"
#         
#         fusion_gene_object = { "type": download_type, "action": "download", "title": "size: " + get_human_size(os.path.dirname(__file__)+ "/data/downloads/fusioncatcher/" + cell_line_id + ".txt"), "label": "FusionCatcher", "url": "downloads/fusioncatcher/" + cell_line_id + ".txt"}
#         
#         virus_object = { "type": download_type, "action": "download", "title": "size: " + get_human_size(os.path.dirname(__file__)+ "/data/downloads/viruses/" + cell_line_id + ".txt"), "label": "Viruses", "url": "downloads/viruses/" + cell_line_id + ".txt"}
#         
#         summary_object = { "type": download_type, "action": "download", "title": "size: " + get_human_size(os.path.dirname(__file__)+ "/data/downloads/summary/" + cell_line_id + ".txt"), "label": "Summary", "url": "downloads/summary/" + cell_line_id + ".txt"}
#         
#         ericscript_object = { "type": download_type, "action": "download", "title": "size: " + get_human_size(os.path.dirname(__file__)+ "/data/downloads/ericscript/" + cell_line_id + ".txt"), "label": "EricScript", "url": "downloads/ericscript/" + cell_line_id + ".txt"}
#         
#         tophat_object = { "type": download_type, "action": "download", "title": "size: " + get_human_size(os.path.dirname(__file__)+ "/data/downloads/tophat/" + cell_line_id + ".txt"), "label": "Tophat-Fusion", "url": "downloads/tophat/" + cell_line_id + ".txt"}
# 
#         oncofuse_fc_file = glob.glob(os.path.dirname(__file__)+ "/data/downloads/oncofuse/oncofuse.fusioncatcher/*/" + cell_line_id + "*")
#         if len(oncofuse_fc_file) > 0: oncofuse_fc_file = oncofuse_fc_file[0]
#         if oncofuse_fc_file:
#             oncofuse_fc = { "type": download_type, "action": "download", "title": "size: " + get_human_size(oncofuse_fc_file), "label": "FusionCatcher", "url": oncofuse_fc_file.replace(os.path.dirname(__file__)+ "/data/", "")}
#         else:
#             oncofuse_fc = { "type": "text", "title": "N/A", "label": "N/A"}
#          
#         oncofuse_es_file = glob.glob(os.path.dirname(__file__)+ "/data/downloads/oncofuse/oncofuse.ericscript/*/*" + cell_line_id + "*")
#         if len(oncofuse_es_file) > 0: oncofuse_es_file = oncofuse_es_file[0]
#         if oncofuse_es_file:
#             oncofuse_es = { "type": download_type, "action": "download", "title": "size: " + get_human_size(oncofuse_es_file), "label": "EricScript", "url": oncofuse_es_file.replace(os.path.dirname(__file__) + "/data/", "")}
#         else:
#             oncofuse_es = { "type": "text", "title": "N/A", "label": "N/A"}
#         
#         oncofuse_th_file = glob.glob(os.path.dirname(__file__)+ "/data/downloads/oncofuse/oncofuse.Tophat-fusion/*/" + cell_line_id + "*")
#         if len(oncofuse_th_file) > 0: oncofuse_th_file = oncofuse_th_file[0]
#         if oncofuse_th_file:
#             oncofuse_th = { "type": download_type, "action": "download", "title": "size: " + get_human_size(oncofuse_th_file), "label": "Tophat-Fusion", "url": oncofuse_th_file.replace(os.path.dirname(__file__) + "/data/", "")}
#         else:
#             oncofuse_th = { "type": "text", "title": "N/A", "label": "N/A"}
#         
#         rows.append([
#             cell_line_name,
#             cell_line_disease,
#             fusion_gene_object,
#             ericscript_object,
#             tophat_object,
#             virus_object,
#             summary_object,
#             oncofuse_fc,
#             oncofuse_es,
#             oncofuse_th
#             ])
#     
#     return HttpResponse(json.dumps(response))

def download_data2(request):

    data = json.loads(request.body)
    print(data)
    
    offset = data["offset"]
    limit = data["limit"]
    
    response = {}
    
    header = [
        {
            "label": "cell_line",
            "title": "Cell line",
            "tooltip": "Cell Line",
            "filters": {
                "title": "Cell line filters:",
                "list": [
                    {
                        "type": "select",
                        "key": "cell_line",
                        "title": "Select a cell line:",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]
            }
        },
        {
            "label": "disease",
            "title": "Disease",
            "tooltip": "Cell line disease",
            "filters": {
                "title": "Disease filters:",
                "list": []
            }
        },
        {
            "label": "fc",
            "title": "FC predicted events",
            "tooltip": "FusionCatcher (FC) predicted events",
            "filters": {
                "title": "FusionCatcher filters:",
                "list": []
            }
        },
        {
            "label": "es",
            "title": "ES predicted events",
            "tooltip": "EricScript (ES) predicted events",
            "filters": {
                "title": "EricScript filters:",
                "list": []
            }
        },
        {
            "label": "th",
            "title": "TF predicted events",
            "tooltip": "Tophat-Fusion (TF) predicted events",
            "filters": {
                "title": "Tophat-Fusion filters:",
                "list": []
            }
        },
#         {
#             "label": "summary",
#             "title": "Summary information",
#             "tooltip": "Summary information",
#             "filters": {
#                 "title": "Summary filters:",
#                 "list": []
#             }
#         },
        {
            "label": "virus",
            "title": "Viruses information",
            "tooltip": "Viruses information (only by FusionCatcher)",
            "filters": {
                "title": "Virus filters:",
                "list": []
            }
        },
        {
            "label": "oncofuse_fc",
            "title": "Oncofuse - FC",
            "tooltip": "Oncofuse - FusionCatcher (FC)",
            "filters": {
                "title": "Oncofuse - FusionCatcher filters:",
                "list": []
            }
        },
        {
            "label": "oncofuse_es",
            "title": "Oncofuse - ES",
            "tooltip": "Oncofuse - EricScript (ES)",
            "filters": {
                "title": "Oncofuse - EricScript filters:",
                "list": []
            }
        },
        {
            "label": "oncofuse_th",
            "title": "Oncofuse - TF",
            "tooltip": "Oncofuse - Tophat-Fusion (TF)",
            "filters": {
                "title": "Oncofuse - Tophat-Fusion filters:",
                "list": []
            }
        }
    ]
    
    rows = []
    
    ccle_map = {}
    ccle_infos = get_ccle_infos()
    for cell_info in ccle_infos:
        ccle_map[cell_info["ID"]] = cell_info
    
    for cell in CellLine.nodes.order_by('name')[offset:offset+limit]:
        cell_line_id = cell.cell_line
        cell_line_name = ccle_map[cell_line_id]["Cell Line"]
        cell_line_disease = ccle_map[cell_line_id]["Disease name"]

        download_type = "download"
        
        fusion_gene_object = { "type": download_type, "action": "download", "title": "size: " + get_human_size(os.path.dirname(__file__)+ "/data/downloads/fusioncatcher/" + cell_line_id + ".txt"), "label": "FusionCatcher", "url": "download/fusioncatcher/" + cell_line_id + ".txt"}
        
        virus_object = { "type": download_type, "action": "download", "title": "size: " + get_human_size(os.path.dirname(__file__)+ "/data/downloads/viruses/" + cell_line_id + ".txt"), "label": "Viruses", "url": "download/viruses/" + cell_line_id + ".txt"}
        
        summary_object = { "type": download_type, "action": "download", "title": "size: " + get_human_size(os.path.dirname(__file__)+ "/data/downloads/summary/" + cell_line_id + ".txt"), "label": "Summary", "url": "download/summary/" + cell_line_id + ".txt"}
        
        ericscript_object = { "type": download_type, "action": "download", "title": "size: " + get_human_size(os.path.dirname(__file__)+ "/data/downloads/ericscript/" + cell_line_id + ".txt"), "label": "EricScript", "url": "download/ericscript/" + cell_line_id + ".txt"}
        
        tophat_object = { "type": download_type, "action": "download", "title": "size: " + get_human_size(os.path.dirname(__file__)+ "/data/downloads/tophat/" + cell_line_id + ".txt"), "label": "Tophat-Fusion", "url": "download/tophat/" + cell_line_id + ".txt"}

        oncofuse_fc_file = glob.glob(os.path.dirname(__file__)+ "/data/downloads/oncofuse/oncofuse.fusioncatcher/*/" + cell_line_id + "*")
        if len(oncofuse_fc_file) > 0: oncofuse_fc_file = oncofuse_fc_file[0]
        if oncofuse_fc_file:
            oncofuse_fc = { "type": download_type, "action": "download", "title": "size: " + get_human_size(oncofuse_fc_file), "label": "FusionCatcher", "url": oncofuse_fc_file.replace(os.path.dirname(__file__)+ "/data/downloads/", "downloads/")}
        else:
            oncofuse_fc = { "type": "text", "title": "N/A", "label": "N/A"}
         
        oncofuse_es_file = glob.glob(os.path.dirname(__file__)+ "/data/downloads/oncofuse/oncofuse.ericscript/*/*" + cell_line_id + "*")
        if len(oncofuse_es_file) > 0: oncofuse_es_file = oncofuse_es_file[0]
        if oncofuse_es_file:
            oncofuse_es = { "type": download_type, "action": "download", "title": "size: " + get_human_size(oncofuse_es_file), "label": "EricScript", "url": oncofuse_es_file.replace(os.path.dirname(__file__) + "/data/downloads/", "downloads/")}
        else:
            oncofuse_es = { "type": "text", "title": "N/A", "label": "N/A"}
        
        oncofuse_th_file = glob.glob(os.path.dirname(__file__)+ "/data/downloads/oncofuse/oncofuse.Tophat-fusion/*/" + cell_line_id + "*")
        if len(oncofuse_th_file) > 0: oncofuse_th_file = oncofuse_th_file[0]
        if oncofuse_th_file:
            oncofuse_th = { "type": download_type, "action": "download", "title": "size: " + get_human_size(oncofuse_th_file), "label": "Tophat-Fusion", "url": oncofuse_th_file.replace(os.path.dirname(__file__) + "/data/downloads/", "downloads/")}
        else:
            oncofuse_th = { "type": "text", "title": "N/A", "label": "N/A"}
        
        single_row = {
            "cell_line": [{
                "type": "text",
                "url": "cell_line/" + cell_line_id,
                "label": cell_line_name,
                "color": "black"
                }],
            "disease": [{
                "type": "text",
                "label": cell_line_disease,
                "color": "black"
                }],
            "fc": [fusion_gene_object],
            "es": [ericscript_object],
            "th": [tophat_object],
            "virus": [virus_object],
#             "summary": [summary_object],
            "oncofuse_fc": [oncofuse_fc],
            "oncofuse_es": [oncofuse_es],
            "oncofuse_th": [oncofuse_th]
        }
        
        rows.append(single_row)
    
    total = len(CellLine.nodes)
    
    response = {"structure": {"field_list": header}, "total": total, "hits": rows}
    
    return HttpResponse(json.dumps(response))

def get_gene_infos():
    info = {}
    
    txt_file = open(os.path.dirname(__file__) + "/gene_start.txt", "r")
    for line in txt_file:
        line = line.rstrip()
        fields = line.split("\t")
        info[fields[0]] = fields[1:]
        #print(fields)
    
    return info

def get_ccle_infos():
    # header = ["ID", "Cell Line", "Disease", "Disease name"]
    rows = []
    
    txt_file = open(os.path.dirname(__file__) + "/data/disease_links.tsv", "r")
    disease_map = {}
    txt_file.readline()
    for line in txt_file:
        line = line.rstrip()
        words = line.split("\t")
        words[0] = words[0].replace("_", " ")
        disease_map[words[0]] = words[1]
        
    #print(disease_map)
    
    map = {}
    txt_file = open(os.path.dirname(__file__) + "/ccle_ids.txt", "r")
    next(txt_file)
    for line in txt_file:
        line = line.rstrip()
        words = line.split("\t")
        object = { "ID": words[0].rstrip(), "Cell Line": words[1], "Disease": words[2], "Disease name": words[3]}
        
        object["Disease link"] = disease_map[words[3]] if words[3] in disease_map else None
        
        ccle_key = object["Cell Line"].lower()
        
        map[ccle_key] = object
        rows.append(object)
        
    txt_file = open(os.path.dirname(__file__) + "/data/Cell_Lines_Details-1.csv", "r")
    header2 = txt_file.readline().rstrip().split("\t")
#     header += header2
    for line in txt_file:
        line = line.rstrip()
        words = line.split("\t")
        words[0] = words[0].replace("-", "").lower()
        if words[0] not in map: continue
        object = map[words[0]]
        
        for i in range(1, len(header2)):
            object[header2[i]] = words[i]
        
    txt_file = open(os.path.dirname(__file__) + "/data/Cell_Lines_Details-2.csv", "r")
    header3 = txt_file.readline().rstrip().split("\t")
    del header3[1]
#     header += header3
    for line in txt_file:
        line = line.rstrip()
        words = line.split("\t")
        del words[1]
        words[0] = words[0].replace("-", "")
        if words[0] not in map: continue
        
        object = map[words[0]]
        for i in range(1, len(header3)-1):
            object[header3[i]] = words[i]
    
    for obj in rows:
        if "COSMIC identifier" in obj:
            obj["Drug resistance"] = obj["COSMIC identifier"]
            
    return rows

def get_cell_line_from_disease(disease):
    ccle_infos = get_ccle_infos()
    
    cls = []
    for ccle in ccle_infos:
        if ccle["Disease"] == disease:
            cls.append(ccle["ID"])
    
    return cls

def diseases_simple(request, prefix = None):
    response = []
    
    disease_dict = {}
    
    counter = Counter()
    for cell_line in CellLine.nodes:
        disease = cell_line.disease
        disease_name = cell_line.disease_name
        
        if prefix is None or prefix.lower() not in disease_name.lower(): continue
        
        disease_dict[disease] = {"id": disease, "label": disease_name, "img": "disease-icon.png"}
        counter[disease] += 1
        
    for disease in disease_dict.values():
        n = counter[disease["id"]]
        disease["cell_lines"] = n
        disease["extra"] = str(n) + " cell line" + ("s" if n>1 else "")
    
    diseases = disease_dict.values()
#     diseases.sort(key=lambda x: x["cell_lines"], reverse=True)
    diseases.sort(key=lambda x: x["label"])
    
    for disease in diseases:
        response.append(disease)
    
    return HttpResponse(json.dumps(response))

def get_distribution(request, node1, node2, howmany, sorting):
    
    if not howmany: howmany = -1
    howmany = int(howmany)
    
    response = {}
    labels = []
    header = []
    items = []
    
    labels_translation = get_ccle_infos()
    
    lines = open(os.path.dirname(__file__) + "/statistics/" + node1 + "_" + node2 + ".csv")
    line_no = 0
    for line in lines:
        line_no += 1
        line = line.rstrip()
        fields = line.split(",")
        if line_no == 1:
            header = fields[0:2]
            continue
        
        items.append([fields[0], fields[1]])
    lines.close()
    
    if sorting == "DESC" or sorting == "ASC":
        items = sorted(items, key=lambda item: int(item[1]))
        if sorting == "DESC":
            items = list(reversed(items))
    else:
        items = sorted(items, key=lambda item: (not item[0].isdigit(), item[0].zfill(3)))

    if howmany >= 0 and len(items) > howmany:
        items = items[:howmany]
    
    labels = list(labels)
    
    if node1 == "CellLine":
        for i,item in enumerate(items):
            for ccle in labels_translation:
                if ccle["ID"] == item[0]:
                    value = node2
                    if value == "FusionCell": value = "fusion events"
                    item[0] = {"id": ccle["ID"], "label": ccle["Cell Line"], "name": "Cell line: " + ccle["Cell Line"] + " / Disease: "+ccle["Disease name"], "value": "Cell line " + value}
                    break
    elif node1 == "Disease":
        for i, item in enumerate(items):
            for ccle in labels_translation:
                if ccle["Disease"] == item[0]:
                    value = node2
                    if value == "FusionCell": value = "Fusion events"
                    if value == "CellLine": value = "Cell lines"
                    item[0] = {"id": ccle["Disease"], "label": ccle["Disease"], "name": "Disease: "+ ccle["Disease name"] + " ("+ccle["Disease"]+")", "value": value + " per disease"}
                    break
    elif node1 == "Chromosome":
        for i,item in enumerate(items):
            value = node2
            if value == "FusionCellRatio": value = "Normalized GFE count"
            item[0] = {"id": item[0], "label": "Chromosome " + item[0], "name": "Chromosome " + item[0], "value": value}
    elif node1 == "PredictedEffect":
        for i,item in enumerate(items):
            item[0] = {"id": item[0], "name": "Predicted effect: " + item[0], "label": item[0], "value": "Number of breakpoints"}
    
    response['details'] = {"labels": labels, "header": header, "items": items}
    
    return HttpResponse(json.dumps(response))

from decimal import Decimal
def get_single_distribution(request, label, value):
    
    response = {}
    header = []
    items = []
    item = []
    
    files = sorted(glob.glob(os.path.dirname(__file__) + "/statistics/" + label + "_*.csv"))
    
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
                element = {"type": "text", "label": header_fields[1]}
                if len(header_fields) > 2:
                    element["tooltip"] = header_fields[2]
                    
                v = fields[1]
                if len(header_fields) > 3:
                    if header_fields[3] == "scientific":
                        v = '%.2E' % Decimal(v)
                    
                header.append(element)
                item.append(v)
                break

        lines.close()
        
    response['details'] = {"header": header, "items": [item]}
    
    return HttpResponse(json.dumps(response))

def cell_lines_of_fusion(request, fus_id):
    response = {}
    header = ["Cell line name", "Disease"]
    
    fusions = []
    
    fusions.append(Fusion.nodes.get(fusion_id = fus_id))
    
    ccle_map = {}
    ccle_infos = get_ccle_infos()
    for cell_info in ccle_infos:
        ccle_map[cell_info["ID"]] = cell_info
    
    rows = []
    met = set()
    for fusion in fusions:
        for cell_line in fusion.fromCellLineToFusion:
            c_line = cell_line.cell_line
            if c_line in met: continue
            met.add(c_line)
            
            cell_info = ccle_map[c_line]
            cell_line_name = cell_info["Cell Line"] + " ("+cell_info["ID"]+")"
            cell_line_disease = cell_info["Disease name"] + " ("+cell_info["Disease"]+")"
            item = [cell_line_name, cell_line_disease]
            rows.append(item)
    
    response['details'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

annotations = {
            "bodymap2": {
                "name": "Illumina Body Map 2.0",
                "url": "http://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-513/",
                "title": "These data consist of 16 RNA samples from 16 different organs from healthy people. If a candidate couple falls into this list, it is tagged as known false positive fusion genes",
                "img_src": "illumina.body.map.png"
            },
            "ticdb": {
                "name": "TICdb",
                "url": "http://www.unav.es/genetica/TICdb/",
                "title": "TICdb is a database of 1374 Translocation breakpoints found in human tumors, involving 431 different genes.",
                "img_src": "http://www.unav.es/genetica/logo.jpg"
            },
            "chimerdb": {
                "name": "ChimerDB",
                "url": "http://203.255.191.229:8080/chimerdbv31/mindex.cdb",
                "title": "Literature-based annotation",
                "img_src": "chimerdb.logo.png"
            },
            "cell_lines": {
                "name": "Cell lines",
                "url": "http://cancer.sanger.ac.uk/cell_lines",
                "title": "Known fusion gene from paper: C. Klijn et al., A comprehensive transcriptional portrait of human cancer cell lines, Nature Biotechnology, Dec. 2014, DOI:10.1038/nbt.3080",
                "img_src": "http://cancer.sanger.ac.uk/images/domain_logos/cell_lines_domain_logo_60x60.png"
            },
            "1000genomes": {
                "name": "1000 genomes",
                "url": "http://www.internationalgenome.org/",
                "title": "Fusion gene has been seen in a healthy sample. It has been found in RNA-seq data from some samples from 1000 genomes project.",
                "img_src": "https://www.genome.gov/images/content/1000genomes.jpg"
            },
#             "18cancers": {
#                 "name": "18 cancers",
#                 "url": "http://www.pnas.org/content/113/48/13768",
#                 "title": "Fusion gene found in a RNA-seq dataset of 18 types of cancers from 600 tumor samples.",
#                 "img_src": ""
#             },
            "gtex": {
                "name": "gtex",
                "url": "https://www.gtexportal.org/home/",
                "title": "Fusion gene has been seen in a healthy sample. It has been found in GTEx database of healthy tissues.",
                "img_src": "http://phoenixrising.me/wp-content/uploads/gtex.jpg"
            },
#             "hpa": {
#                 "name": "hpa",
#                 "url": "http://www.mcponline.org/content/13/2/397",
#                 "title": "Fusion gene has been seen in a healthy sample. It has been found in RNA-seq database of 27 healthy tissues.",
#                 "img_src": ""
#             },
            "known": {
                "name": "known",
                "url": "",
                "title": "Fusion gene which has been previously reported or published in scientific articles/reports/books/abstracts/databases",
                "img_src": "https://learn.nlm.nih.gov//images/pubmed_logo_large.gif"
            },
            "tcga": {
                "name": "TCGA",
                "url": "https://gdc.cancer.gov/",
                "title": "Known fusion gene from the TCGA database",
                "img_src": "http://ciriellolab.org/images/tcga_logo.gif"
            }
            }

def build_rows2(fusions, c_line = "ALL", include_FC_only = False):
    global annotations
    
    rows = []
    
    print("BUILDING ROWS (for cell_line(s)="+c_line+") include_FC_only:"+str(include_FC_only))
    gold_pairs = get_gold_pairs()
#     gold_genes = get_gold_genes()
    
    oncofuse = get_oncofuse_info_simple()
    oncofuse_set = set()
    for item in oncofuse:
        oncofuse_set.add("#".join(item))
        
    # cell_name = CellLine.nodes.get(cell_line=myfusion.cell_line).name
    names = {}
    for c in CellLine.nodes:
        names[c.cell_line] = c.name
    
    print("ANALYZING "+str(len(fusions))+" FUSIONS")
    for myfusion in fusions:
#         cell_lines = []
#         if c_line == "ALL":
#             cell_lines = []
#             met = set()
#             for c in myfusion.cell_lines:
#                 if c.cell_line in met: continue
#                 met.add(c.cell_line)
#                 cell_lines.append(c)
#         else:
#             cell_lines = [CellLine.nodes.get(cell_line__exact=c_line)]
#         
#         print(cell_lines)
        
        disease = ""
        acronym = ""
        
        gold_pair = None
        for pair in gold_pairs["items"]:
            if set([pair[0], pair[1]]) == set([myfusion.gene1, myfusion.gene2]):
                gold_pair = pair
                break
            
        gold_gene1 = Gene.nodes.get(symbol=myfusion.gene1).census
#         gold_gene1 = None
#         for gene in gold_genes["items"]:
#             if gene[0] == myfusion.gene1:
#                 gold_gene1 = Gene.nodes.get(symbol=myfusion.gene1)
#                 break

        gold_gene2 = Gene.nodes.get(symbol=myfusion.gene2).census
#         gold_gene2 = None
#         for gene in gold_genes["items"]:
#             if gene[0] == myfusion.gene2:
#                 gold_gene2 = Gene.nodes.get(symbol=myfusion.gene2)
#                 break
      
        NA_object = {
                    "type": "text",
                    "label": "-",
                    "color": "black",
                    "title": "-"
                }
        
        fc_flag = NA_object
        es_flag = NA_object
        th_flag = NA_object
        ja_flag = NA_object
        
        fusioncatcher_events = "FC" in myfusion.algos
        ericscript_events = "ES" in myfusion.algos
        tophat_events = "TH" in myfusion.algos
        jaffa_events = "JA" in myfusion.algos
        
#         fusioncatcher_events = []
#         ericscript_events = []
#         tophat_events = []
#         for fusion_point in myfusion.fusion_points:
#             print(fusion_point)
#             for algorithm in fusion_point.algorithms:
#                 for fc in algorithm.fusioncatcher_events:
#                     fusioncatcher_events.append(fc)
#                          
#                 for es in algorithm.ericscript_events:
#                     ericscript_events.append(es)
#                          
#                 for th in algorithm.tophat_events:
#                     tophat_events.append(th)
                    
        if fusioncatcher_events:
             
            fc_flag = {
                "type": "button",
                "color": "lightblue",
                "action": "window",
                "title": "This fusion event is supported by FusionCatcher",
                "items": [{
                    "type": "text",
                    "label": "FC"
                }],
                "card": {
                    "title": "Detailed information given by FusionCatcher for gene pair: " + str(myfusion.gene1) + "/" + str(myfusion.gene2),
                    "width": "100",
                    "elements": [
                            {
                                    "type": "table",
                                    "data": {
                                            "url": "/fusion_api/fusion/fusioncatcher/"+str(myfusion.fusion_cell_id) + "/" + str(myfusion.cell_line)
                                    }
                            }
                    ]
                    }
            }
                 
        if ericscript_events and not include_FC_only:
             
            es_flag = {
                "type": "button",
                "color": "lightgreen",
                "action": "window",
                "title": "This fusion event is supported by EricScript",
                "items": [{
                    "type": "text",
                    "label": "ES"
                }],
                "card": {
                    "title": "Detailed information given by EricScript for gene pair: " + str(myfusion.gene1) + "/" + str(myfusion.gene2),
                    "width": "100",
                    "elements": [
                            {
                                    "type": "table",
                                    "data": {
                                            "url": "/fusion_api/fusion/ericscript/"+str(myfusion.fusion_cell_id) + "/" + str(myfusion.cell_line)
                                    }
                            }
                    ]
                    }
            }
                 
        if tophat_events and not include_FC_only:
            th_flag = {
                "type": "button",
                "color": "lightcoral",
                "action": "window",
                "title": "This fusion event is supported by Tophat-Fusion",
                "items": [{
                    "type": "text",
                    "label": "TF"
                }],
                "card": {
                    "title": "Detailed information given by Tophat-Fusion for gene pair: " + str(myfusion.gene1) + "/" + str(myfusion.gene2),
                    "width": "100",
                    "elements": [
                            {
                                    "type": "table",
                                    "data": {
                                            "url": "/fusion_api/fusion/tophat/"+str(myfusion.fusion_cell_id) + "/" + str(myfusion.cell_line)
                                    }
                            }
                    ]
                    }
            }
            
        if jaffa_events and not include_FC_only:
            ja_flag = {
                "type": "button",
                "color": "chocolate",
                "action": "window",
                "title": "This fusion event is supported by Jaffa",
                "items": [{
                    "type": "text",
                    "label": "JA"
                }],
                "card": {
                    "title": "Detailed information given by Jaffa for gene pair: " + str(myfusion.gene1) + "/" + str(myfusion.gene2),
                    "width": "100",
                    "elements": [
                            {
                                    "type": "table",
                                    "data": {
                                            "url": "/fusion_api/fusion/jaffa/"+str(myfusion.fusion_cell_id) + "/" + str(myfusion.cell_line)
                                    }
                            }
                    ]
                    }
            }
        
        is_gold_pair_element = {
                    "type": "text",
                    "label": myfusion.gene_couple[0].cosmic,
                    "color": "black"
                }
        if gold_pair is not None:
            is_gold_pair_element = {
                "type": "image",
                "img_src": "cosmic_logo.png",
                "width": "30px",
                "title": "This pair of genes is known to be in the COSMIC repository.",
                "target": "_blank",
                "url": gold_pair[4]
            }
            
        gene1_object = {
                    "type": "text",
                    "label": myfusion.gene1,
                    "color": "black"
                }
        if gold_gene1:            
            gene1_object["pedix"] = [
                {
                    "type": "image",
                    "title": "The cancer Gene Census is an ongoing effort to catalogue those genes for which mutations have been causally implicated in cancer.",
                    "img_src": "census-icon.png",
                    "width": "25px",
                    "url": "http://cancer.sanger.ac.uk/cosmic/census",
                    "target": "_blank"
                }
            ]
                            
        gene2_object = {
                    "type": "text",
                    "label": myfusion.gene2,
                    "color": "black"
                }
        if gold_gene2:
            gene2_object["pedix"] = [
                {
                    "type": "image",
                    "title": "The cancer Gene Census is an ongoing effort to catalogue those genes for which mutations have been causally implicated in cancer.",
                    "img_src": "census-icon.png",
                    "width": "25px",
                    "url": "http://cancer.sanger.ac.uk/cosmic/census",
                    "target": "_blank"
                }
            ]
        
        color_counter = Counter()
        color_tags = {"red": set(), "orange": set(), "green": set()}
        color_messages = {"red": "false positive event with high probability", "orange": "false positive event with medium probability", "green": "true positive event with high probability"}
        color = None
        description = ""
        annotations_found = set()
        all_annotations = myfusion.annotations if myfusion.annotations is not None else []
#         for fc_description in all_annotations:
#             if fcfusion.description is None: continue
#             if len(fcfusion.description) == 0: continue
#             if fcfusion.description[0] == " ": continue
            
#             for tag in fcfusion.description:
#                 annotations_found.add(tag)
        in_mitelman = "mitelman" in all_annotations
        
        is_mitelman_element = {}
        if in_mitelman:
            is_mitelman_element = {
                "type": "image",
                "img_src": "mitelman_logo.png",
                "width": "30px",
                "title": "This pair of genes is known to be in the Mitelman repository.",
                "target": "_blank"
            }

        # Check red list
        descriptions_in_red_list = filter(lambda x: x in red_list, all_annotations)
#         print("\tRED:" + descriptions_in_red_list)
        if len(descriptions_in_red_list) > 0:
            color = "red"
            color_counter[color] += 1
            for c in descriptions_in_red_list: 
                color_tags[color].add(c)
            description = "false positive event with high probability (because of the following tags: " + ",".join(descriptions_in_red_list) + ")."
        
        # Check orange list
        descriptions_in_orange_list = filter(lambda x: x in orange_list, all_annotations)
#         print("\tORANGE:" + descriptions_in_orange_list)
        if len(descriptions_in_orange_list) > 0:
            color = "orange"
            color_counter[color] += 1
            for c in descriptions_in_orange_list: 
                color_tags[color].add(c)
            description = "false positive event with medium probability (because of the following tags: " + ",".join(descriptions_in_orange_list) + ")."
            
        # Check green list
        descriptions_in_green_list = filter(lambda x: x in green_list, all_annotations)
#         print("\tGREEN:" + descriptions_in_green_list)
        if len(descriptions_in_green_list) > 0:
            color = "green"
            color_counter[color] += 1
            for c in descriptions_in_green_list: 
                color_tags[color].add(c)
            description = "true positive event with high probability (because of the following tags: " + ",".join(descriptions_in_green_list) + ")."
        
        if len(color_counter) == 1:
            annotation = {"type": "icon", "icon_img": "fa-circle", "modifiers": "fa-2x", "color": color}
            fc_flag["title"] = (fc_flag["title"] + " and is a " + description if "title" in fc_flag else description.capitalize())
            if "items" not in fc_flag: fc_flag["items"] = []
            fc_flag["items"].append(annotation)
            
        elif len(color_counter) > 1:
            annotation = {"type": "icon", "icon_img": "fa-circle", "modifiers": "fa-2x", "color": "grey"}
            
            counting_messages = []
            for color in color_counter:
                counting_messages.append(color_messages[color] + " because of the following tags: [" + ",".join(x.replace("<", "&lt;").replace(">", "&gt;") for x in color_tags[color]) + "]")
            counting_message = ", ".join([str(color_counter[key]) + " " + key for key in color_counter])            
            fc_flag["title"] = (fc_flag["title"] + " and is ambiguously tagged (" + ", ".join(counting_messages) + " events)." if "title" in fc_flag else description.capitalize())
            if "items" not in fc_flag: fc_flag["items"] = []
            fc_flag["items"].append(annotation)
                    
        
        annotations_found_also_in_gold = []
        for tag in all_annotations:
            annotation_query = tag
            if annotation_query.startswith("chimerdb"): annotation_query = "chimerdb"
            if annotation_query not in annotations: continue
            
            annotations_found_also_in_gold.append(tag)
            annotation_metadata = annotations[annotation_query]
            
            if "img_src" not in annotation_metadata or annotation_metadata["img_src"] == "":
                annotation = {"type": "text", "label": annotation_metadata["name"], "title": annotation_metadata["title"]}
            else:
                annotation = {"type": "image", "width": "50px", "title": annotation_metadata["title"], "img_src": annotation_metadata["img_src"]}
            
            # fc_flag["title"] = (fc_flag["title"] + " and is a " + description if "title" in fc_flag else description.capitalize())
            if "items" not in fc_flag: fc_flag["items"] = []
            fc_flag["items"].append(annotation)
            
        if annotations_found_also_in_gold:
            fc_flag["title"] += " It has been found in the following databases: " + ", ".join(annotations_found_also_in_gold)
        
        # cell_name = CellLine.nodes.get(cell_line=myfusion.cell_line).name
        cell_name = names[myfusion.cell_line]
        
        fc_flag_final = copy.deepcopy(fc_flag)
        if fc_flag:
            if "#".join([myfusion.cell_line, myfusion.gene1, myfusion.gene2, "fusioncatcher"]) in oncofuse_set:
                #fc_flag_final["type"] = "button"
                fc_flag_final["title"] += ". It has also been validated by Oncofuse."
                if "items" not in fc_flag_final: fc_flag_final["items"] = []
                fc_flag_final["items"].append(
                    {
                        "type": "image",
                        "img_src": "oncofuse-icon.png",
                        "width": "50px"
                     })
        
        es_flag_final = copy.deepcopy(es_flag)
        if es_flag:
            if "#".join([myfusion.cell_line, myfusion.gene1, myfusion.gene2, "ericscript"]) in oncofuse_set:
                #es_flag_final["type"] = "button"
                es_flag_final["title"] += ". It has also been validated by Oncofuse."
                if "items" not in es_flag_final: es_flag_final["items"] = []
                es_flag_final["items"].append(
                    {
                        "type": "image",
                        "img_src": "oncofuse-icon.png",
                        "width": "50px"
                     })
                
        th_flag_final = copy.deepcopy(th_flag)
        if th_flag:
            if "#".join([myfusion.cell_line, myfusion.gene1, myfusion.gene2, "tophat"]) in oncofuse_set:
                #th_flag_final["type"] = "button"
                th_flag_final["title"] += ". It has also been validated by Oncofuse."
                if "items" not in th_flag_final: th_flag_final["items"] = []
                th_flag_final["items"].append(
                    {
                        "type": "image",
                        "img_src": "oncofuse-icon.png",
                        "width": "50px"
                     })
                
        ja_flag_final = copy.deepcopy(ja_flag)
        if ja_flag:
            if "#".join([myfusion.cell_line, myfusion.gene1, myfusion.gene2, "jaffa"]) in oncofuse_set:
                #th_flag_final["type"] = "button"
                ja_flag_final["title"] += ". It has also been validated by Oncofuse."
                if "items" not in ja_flag_final: ja_flag_final["items"] = []
                ja_flag_final["items"].append(
                    {
                        "type": "image",
                        "img_src": "oncofuse-icon.png",
                        "width": "50px"
                     })
        
#         fusion_id_object = {
#                 "type": "image",
#                 "img_src": "fusionevent-icon.png",
#                 "label": myfusion.fusion_cell_id,
#                 "color": "grey",
#                 "title": myfusion.fusion_cell_id,
#                 "url": "fusion_event/" + myfusion.fusion_cell_id + "/",
#                 "target": "_blank"
#             }
        fusion_id_object = [
                {
                    "type": "image",
                    "img_src": "fusionevent-icon.png",
                    "label": myfusion.fusion_cell_id,
                    "url": "fusion_event/" + myfusion.fusion_cell_id + "/",
                    "target": "_blank"
                },
                {
                    "type": "text",
                    "label": myfusion.fusion_cell_id,
                    "color": "grey",
                    "url": "fusion_event/" + myfusion.fusion_cell_id + "/",
                    "target": "_blank"
                }
            ]
        
        cell_line_object = {
                "type": "text",
                "label": cell_name,
                "color": "black",
                "url": "cell_line/" + myfusion.cell_line + "/",
                "target": "_blank"
            }
            
            
        single_row = {
            "fusion_id": fusion_id_object,
            "cell_line": [cell_line_object],
            "gene1": [gene1_object],
            "gene2": [gene2_object],
            "cosmic": [is_gold_pair_element],
            "mitelman": [is_mitelman_element],
            "fusioncatcher": [fc_flag_final]
            }
        
        if not include_FC_only:
            single_row["ericscript"] = [es_flag_final]
            single_row["tophat_fusion"] = [th_flag_final]
            single_row["jaffa"] = [ja_flag_final]
            
        rows.append(single_row)
        
    return rows


def build_fc_table(request, fus_id, c_line):
    
    response = {}
    # header = get_fc_header()
    header = ["FusionCatcher information"]
    
#     fusions.append(Fusion.nodes.get(fusion_id = fus_id))
    fusions = [FusionCell.nodes.get(fusion_cell_id = fus_id)]
    
    rows = build_fc_rows(fusions, c_line)
    
    response['details'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def build_es_table(request, fus_id, c_line):
    
    response = {}
    
#     header = get_es_header()
    header = ["EricScript information"]
    
    print("Building ES TABLE with FUS_ID="+fus_id)
    
#     fusions.append(Fusion.nodes.get(fusion_id = fus_id))
    fusions = [FusionCell.nodes.get(fusion_cell_id = fus_id)]

    rows = build_es_rows(fusions, c_line)
    
    response['details'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def build_th_table(request, fus_id, c_line):
    
    response = {}
    header = get_tophat_header()
    
    fusions = [FusionCell.nodes.get(fusion_cell_id = fus_id.replace("_SHARP_", "#"))]
    
    rows = build_tophat_rows(fusions, c_line)
    
    response['details'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))


def build_ja_table(request, fus_id, c_line):
    
    response = {}
    # header = get_fc_header()
    header = ["Jaffa information"]
    
#     fusions.append(Fusion.nodes.get(fusion_id = fus_id))
    fusions = [FusionCell.nodes.get(fusion_cell_id = fus_id)]
    
    rows = build_ja_rows(fusions, c_line)
    
    response['details'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def build_fc_rows(fusions, cellLine):
    
    rows = []
    fc_fusions = []
    
    print(cellLine)
    print(len(fusions))

#     for fus in fusions:
#         for fcfusion in fus.with_fc_script:
#             if fcfusion.toCellLine[0].cell_line == cellLine:
#                 fc_fusions.append(fcfusion)
            
#     print("TOTAL " + str(len(fc_fusions)) + " FUSION CATCHER FUSIONS FOUND")

    oncofuse = get_oncofuse_info_simple()
    print(len(oncofuse))
    oncofuse_set = set()
    for item in oncofuse:
        oncofuse_set.add("#".join(item))
    
    print(next(iter(oncofuse_set)))
    
    descriptions_counter = Counter()
    tags_counter = Counter()
    
    for fusion in fusions:
        gene1 = Gene.nodes.get(symbol=fusion.gene1)
        gene2 = Gene.nodes.get(symbol=fusion.gene2)
        
        for fusion_point in fusion.fusion_points:
            
            chromosome1 = Chromosome.nodes.get(chromosome__exact=fusion_point.chromosome1)
            chromosome2 = Chromosome.nodes.get(chromosome__exact=fusion_point.chromosome2)
            
#             gene1 = chromosome1.of_gene[0]
#             gene2 = chromosome2.of_gene[0]
            
            fusion_point_1 = fusion_point.breakpoint1
            fusion_point_2 = fusion_point.breakpoint2
            
            for myfusion in fusion_point.fusioncatcher_events.filter(cell_line__exact=cellLine):
                
                strand_1 = myfusion.strand_1
                predicted_effect_1 = myfusion.predicted_effect_1
                strand_2 = myfusion.strand_2
                predicted_effect_2 = myfusion.predicted_effect_2
                
                #recupero esoni        
                exon1 = []
                exon2 = []
        
                for exon in myfusion.at_exon:
                    if exon.in_gene.filter(symbol__iexact=gene1.symbol):
                        exon1 = exon
                    if exon.in_gene.filter(symbol__iexact=gene2.symbol):
                        exon2 = exon
                #recupero trascritti e proteine
                transcript_couples = []
                proteins = []
                proteins_object = {"type": "paragraph", "data": {"value": "<b>Proteins</b>: No proteins available for this event"}}
                if myfusion.with_trans_couple:
                    proteins_object = {
                        "type": "multi",
                        "subdata":[
                            {
                                "type": "image",
                                "data": {
                                    "url": "https://neogenomics.com/Portals/0/Images/Icons-Circular/MultiOmyx-Icon-Lav.png",
                                    "width": "50px",
                                    "margin": "10px"
                                }
                            },
                            {
                                "type": "paragraph",
                                "data": {"value": "<b>Proteins and transcripts</b>"}
                            },
                            {
                                "type": "button", "action": "window", "label": str(len(myfusion.with_trans_couple)) + " proteins and transcripts", "card": {
                                "title": "Detailed protein information",
                                "width": "100",
                                "elements": [
                                        {
                                                "type": "table",
                                                "data": {
                                                        "url": "/fusion_api/fusion/fusioncatcher_proteins/"+str(myfusion.fusioncatcher_id)
                                                }
                                        }
                                ]
                            }}
                        ]
                    }
                trans_couples_number = len(myfusion.with_trans_couple)
                
                in_oncofuse_fc = "#".join([cellLine, gene1.symbol, gene2.symbol, "fusioncatcher"]) in oncofuse_set
                oncofuse_object = {"type": "paragraph", "data": {"value": "<b>Oncofuse</b>: No Oncofuse data available for this event"}}
                if in_oncofuse_fc:
                    oncofuse_object = {
                        "type": "multi",
                        "subdata":[
                            {
                                "type": "paragraph",
                                "data": {"value": "<b>Oncofuse data</b>"}
                            },
                            { "type": "button", "color": "lightsteelblue", "action": "window", "label": "See Oncofuse data", "card": {
                                "title": "Detailed Oncofuse information",
                                "width": "100",
                                "elements": [
                                        {
                                                "type": "table",
                                                "data": {
                                                        "url": "/fusion_api/fusion/oncofuse/"+str(cellLine)+"/fusioncatcher/"+gene1.symbol+"|"+gene2.symbol
                                                }
                                        }
                                ]
                            }}
                        ]
                     }
                
                row = []
                row.append({
                        "type": "multi",
                        "subdata": [
                            {
                                "type": "image",
                                "data": {
                                        "url": "https://d30y9cdsu7xlg0.cloudfront.net/png/43386-200.png",
                                        "width": "50px",
                                    }
                            },
                            {"type": "paragraph", "data": {"value": "<b>5' position</b>: " + "chr"+ chromosome1.chromosome+":"+str(fusion_point_1)+":"+strand_1}}
                        ]
                        })
                row.append({"type": "spacer"})
                row.append({
                        "type": "multi",
                        "subdata": [
                            {
                                "type": "image",
                                "data": {
                                        "url": "https://d30y9cdsu7xlg0.cloudfront.net/png/43386-200.png",
                                        "width": "50px",
                                    }
                            },
                            {"type": "paragraph", "data": {"value": "<b>3' position</b>: " + "chr"+ chromosome2.chromosome+":"+str(fusion_point_2)+":"+strand_2}}
                        ]
                        })
                
                row.append({"type": "spacer"})
                row.append({"type": "spacer"})
                
                if(predicted_effect_1 != predicted_effect_2):
                    row.append({"type": "paragraph", "data": {"value": "<b>Predicted effect</b>: " + predicted_effect_1+"/"+predicted_effect_2}})
                else:
                    row.append({"type": "paragraph", "data": {"value": "<b>Predicted effect</b>: " + predicted_effect_1}})
                row.append({"type": "spacer"})
                    
                row.append(oncofuse_object)
                
                descriptions = []
                if type(myfusion.description) is list:
                    for description in myfusion.description:
                        if description != " ":
                            descriptions.append(description)
                else:
                    description = myfusion.description
                    if description is not None: description = description.strip()
                    
                    if description is not None and description != "":
                        descriptions = [description]
                        
                for description in descriptions:
                    tags_counter[description] += 1
                    descriptions_counter[",".join(sorted(descriptions))] += 1
                
                tags = []
                for description in descriptions:
                    
                    color = "white"
                    
                    if description in red_list: color = "red"
                    elif description in orange_list: color = "orange"
                    elif description in green_list: color = "green"
                    
                    print(description, color)
                    
                    tags.append({
                        "type": "chip",
                        "value": description,
                        "background_color": color,
                        "color": "black" if color is not "green" else "white",
                    })
                
                if not tags:
                    row.append({"type": "paragraph", "data": {"value": "<b>Tags</b>: No tags available for this event"}})
                else:
                    row.append({
                        "type": "multi",
                        "layout": "row",
                        "align": "start center",
                        "subdata":[
                            {"type": "paragraph", "data": {"value": "<b>Tags</b>:"}},
                            {"type": "chips", "items": tags}
                        ]
                    })
                row.append({"type": "spacer"})
                
                row.append({"type": "paragraph", "data": {"value": "<b>Common mapping reads</b>: " + str(myfusion.common_mapping_reads) }})
                row.append({"type": "paragraph", "data": {"value": "<b>Spanning pairs</b>: " + str(myfusion.spanning_pairs)}})
                row.append({"type": "paragraph", "data": {"value": "<b>Spanning pairs (unique)</b>: " + str(myfusion.spanning_unique_reads)}})
                row.append({"type": "paragraph", "data": {"value": "<b>Longest anchor found</b>: " + str(myfusion.longest_anchor_found)}})
                row.append({"type": "spacer"})
                
                algorithms = []
                for algo in myfusion.fusion_finding_method:
                    algorithms.append({
                        "type": "chip",
                        "value": algo,
                    })
                    
                row.append({
                    "type": "multi",
                    "layout": "row",
                    "align": "start center",
                    "subdata":[{"type": "paragraph", "data": {"value": "<b>Fusion finding method</b>:"}}, {"type": "chips", "items": algorithms}]})
                row.append({"type": "spacer"})
                 
                sequence_object = {
                    "type": "multi",
                    "layout": "row",
                    "align": "start center",
                    "subdata":[
                        {
                            "type": "image",
                            "data": {
                                "url": "http://www.bioscience.co.uk/userfiles/icons/Molecular%20Biology.png",
                                "width": "50px"
                            }
                        },
                        {
                            "type": "paragraph", "data": {"value": "<b>Fusion sequence</b>"}
                        },
                        { "type": "button", "action": "window", "label": "See sequence", "card": {
                            "title": "Sequence",
                            "width": "100",
                            "elements": [
                                    {
                                            "type": "table",
                                            "data": {
                                                    "url": "/fusion_api/fusion/sequence/"+str(myfusion.fusioncatcher_id)
                                            }
                                    }
                            ]
                            }
                         }
                    ]
                }
                row.append(sequence_object)
                row.append({"type": "spacer"})
                
                row.append(proteins_object)
                row.append({"type": "spacer"})
                
                rows.append([{
                        "type": "multi",
                        "subdata": row
                    }])
    
    print("FINISHED BUILD FC ROWS")
    return rows

def build_es_rows(fusions, cellLine):
    
    print("Building ES rows " + str(len(fusions)))
    
    rows = []
    es_fusions = []
    
    oncofuse = get_oncofuse_info_simple()
    oncofuse_set = set()
    for item in oncofuse:
        oncofuse_set.add("#".join(item))
    
#     for fus in fusions:
#         for esfusion in fus.with_eric_script:
#             es_fusions.append(esfusion)
    
    for fusion in fusions:
        print(fusion)
        
        gene1 = Gene.nodes.get(symbol=fusion.gene1)
        gene2 = Gene.nodes.get(symbol=fusion.gene2)
        
        for fusion_point in fusion.fusion_points:
            
            chromosome1 = Chromosome.nodes.get(chromosome__exact=fusion_point.chromosome1)
            chromosome2 = Chromosome.nodes.get(chromosome__exact=fusion_point.chromosome2)
            
#             gene1 = chromosome1.of_gene[0]
#             gene2 = chromosome2.of_gene[0]
            
            breakpoint1 = fusion_point.breakpoint1
            breakpoint2 = fusion_point.breakpoint2
            
            for ericscript_event in fusion_point.ericscript_events.filter(cell_line__exact=cellLine):
    
                print(ericscript_event)
                strand1 = ericscript_event.strand_1
                strand2 = ericscript_event.strand_2
    
                crossing_reads = ericscript_event.crossing_reads
                spanning_reads = ericscript_event.spanning_reads
                mean_intersize = ericscript_event.mean_intersize
                homology = ericscript_event.homology
                fusion_type = ericscript_event.fusion_type
                junction_sequence = {
                        "type": "multi",
                        "subdata":[
                            {
                                "type": "image",
                                "data": {
                                    "url": "http://www.bioscience.co.uk/userfiles/icons/Molecular%20Biology.png",
                                    "width": "50px"
                                }
                            },
                            {
                                "type": "paragraph",
                                "data": {"value": "<b>Sequence</b>"}
                            },
                            { "type": "button", "color": "lightsteelblue", "action": "window", "label": "See sequence", "card": {
                                "title": "Junction sequence",
                                "width": "100",
                                "elements": [
                                        {
                                                "type": "table",
                                                "data": {
                                                        "url": "/fusion_api/fusion/junction/"+str(ericscript_event.ericscript_id)
                                                }
                                        }
                                ]
                            }}
                        ]
                     }
                
                gene_expr_1 = ericscript_event.gene_expr_1
                gene_expr_2 = ericscript_event.gene_expr_2
                gene_expr_fused = ericscript_event.gene_expr_fused
                es = ericscript_event.es
                gjs = ericscript_event.gjs
                us = ericscript_event.us
                eric_score = ericscript_event.eric_score
                
                in_oncofuse_es = "#".join([cellLine, gene1.symbol, gene2.symbol, "ericscript"]) in oncofuse_set
                oncofuse_object = None
                if in_oncofuse_es:
                    oncofuse_object = {
                        "type": "multi",
                        "subdata":[
                            {
                                "type": "paragraph",
                                "data": { "value": "<b>Oncofuse data</b>" }
                            },
                            { "type": "button", "color": "lightsteelblue", "action": "window", "label": "See Oncofuse data", "card": {
                                "title": "Detailed Oncofuse information",
                                "width": "100",
                                "elements": [
                                        {
                                                "type": "table",
                                                "data": {
                                                        "url": "/fusion_api/fusion/oncofuse/"+str(cellLine)+"/ericscript/"+gene1.symbol+"|"+gene2.symbol
                                                }
                                        }
                                ]
                            }}
                        ]
                     }
                
                row = []
                row.append({
                        "type": "multi",
                        "subdata": [
                            {
                                "type": "image",
                                "data": {
                                        "url": "https://d30y9cdsu7xlg0.cloudfront.net/png/43386-200.png",
                                        "width": "50px",
                                    }
                            },
                            {"type": "paragraph", "data": {"value": "<b>5' position</b>: " + "chr"+ chromosome1.chromosome+":"+str(breakpoint1)+":"+strand1}}
                        ]
                        })
                row.append({"type": "spacer"})
                row.append({
                        "type": "multi",
                        "subdata": [
                            {
                                "type": "image",
                                "data": {
                                        "url": "https://d30y9cdsu7xlg0.cloudfront.net/png/43386-200.png",
                                        "width": "50px",
                                    }
                            },
                            {"type": "paragraph", "data": { "value": "<b>3' position</b>: " + "chr"+ chromosome2.chromosome+":"+str(breakpoint2)+":"+strand2}}
                        ]
                        })
                row.append({"type": "spacer"})
                
                row.append(oncofuse_object)
                
                row.append({"type": "paragraph", "data": { "value": "<b>Crossing reads</b>: " + str(crossing_reads)}})
                row.append({"type": "paragraph", "data": { "value": "<b>Spanning reads</b>: " + str(spanning_reads)}})
                row.append({"type": "paragraph", "data": { "value": "<b>Mean insert size</b>: " + str(mean_intersize)}})
                row.append({"type": "paragraph", "data": { "value": "<b>Homology</b>: " + str(homology)}})
                row.append({"type": "paragraph", "data": { "value": "<b>Fusion type</b>: " + str(fusion_type)}})
                
                row.append({"type": "spacer"})
                row.append(junction_sequence)
                row.append({"type": "spacer"})
                
                row.append({"type": "paragraph", "data": { "value": "<b>Gene expr 1</b>: " + str(gene_expr_1)}})
                row.append({"type": "paragraph", "data": { "value": "<b>Gene expr 2</b>: " + str(gene_expr_2)}})
                row.append({"type": "paragraph", "data": { "value": "<b>Gene expr fused</b>: " + ("Infinity" if math.isinf(gene_expr_fused) else str(gene_expr_fused))}})
                row.append({"type": "paragraph", "data": { "value": "<b>Es</b>: " + str(es)}})
                row.append({"type": "paragraph", "data": { "value": "<b>GJS</b>: " + str(gjs)}})
                row.append({"type": "paragraph", "data": { "value": "<b>US</b>: " + str(us)}})
                row.append({"type": "paragraph", "data": { "value": "<b>EricScore</b>: " + str(eric_score)}})
#                     row.append(gene_expr_1)
#                     row.append(gene_expr_2)
#                     row.append("Infinity" if math.isinf(gene_expr_fused) else gene_expr_fused)
#                     row.append(es)
#                     row.append(gjs)
#                     row.append(us)
#                     row.append(eric_score)
            
                rows.append([{
                        "type": "multi",
                        "subdata": row
                    }])
                    
# {"type": "text", "label": "Chromosome:breakpoint:strand", "tooltip": "It is the pair of chromosomal position of the 5' and 3' end of fusion junction (chromosome:position:strand); 1-based coordinate."},
#         {"type": "text", "label": "Oncofuse", "tooltip": "An indicator that shows whether this event is supported also by Oncofuse."},
#         {"type": "text", "label": "Crossing reads", "tooltip": "The number of paired end discordant reads."},
#         {"type": "text", "label": "Spanning reads", "tooltip": "The number of paired end reads spanning the junction."},
#         {"type": "text", "label": "Mean insert size", "tooltip": "Mean of insert sizes of crossing + spanning reads."},
#         {"type": "text", "label": "Homology", "tooltip": "If filled, all the homologies between the fusion junction and Ensembl genes."},
#         {"type": "text", "label": "Fusion type", "tooltip": "Intra-chromosomal, inter-chromosomal, read-through or CIS."},
#         {"type": "text", "label": "Junction sequence", "tooltip": "Predicted junction fusion sequence."},
#         {"type": "text", "label": "Gene expr 1", "tooltip": "Read count based estimation of the expression level of 5' gene."},
#         {"type": "text", "label": "Gene expr 2", "tooltip": "Read count based estimation of the expression level of 3' gene."},
#         {"type": "text", "label": "Gene expr fused", "tooltip": "Read count based estimation of the expression level of the predicted chimeric transcript."},
#         {"type": "text", "label": "Es", "tooltip": "Edge score."},
#         {"type": "text", "label": "GJS", "tooltip": "Genuine Junction score."},
#         {"type": "text", "label": "US", "tooltip": "Uniformity score."},
#         {"type": "text", "label": "EricScore", "tooltip": "EricScore score (adaboost classifier)."}
    print(len(rows))
    
    return rows

def build_tophat_rows(fusions, cellLine):
    rows = []
    th_fusions = []
    
#     for fus in fusions:
#         for thfusion in fus.with_tophat_script:
#             th_fusions.append(thfusion)

    oncofuse = get_oncofuse_info_simple()
    oncofuse_set = set()
    for item in oncofuse:
        oncofuse_set.add("#".join(item))
    
    for fusion in fusions:
        
        gene1 = Gene.nodes.get(symbol=fusion.gene1)
        gene2 = Gene.nodes.get(symbol=fusion.gene2)
        
        for fusion_point in fusion.fusion_points:
            
            chromosome1 = Chromosome.nodes.get(chromosome__exact=fusion_point.chromosome1)
            chromosome2 = Chromosome.nodes.get(chromosome__exact=fusion_point.chromosome2)
            
#             gene1 = chromosome1.of_gene[0]
#             gene2 = chromosome2.of_gene[0]
            
            left_coord = fusion_point.breakpoint1
            right_coord = fusion_point.breakpoint2
            
            for tophat_event in fusion_point.tophat_events.filter(cell_line__exact=cellLine):

                spanning_reads = tophat_event.spanning_reads
                spanning_mate_pairs = tophat_event.spanning_mate_pairs
                spanning_mate_pairs_end = tophat_event.spanning_mate_pairs_end
                
                in_oncofuse_th = "#".join([cellLine, gene1.symbol, gene2.symbol, "tophat"]) in oncofuse_set
                print("In oncofuse?", in_oncofuse_th, "#".join([cellLine, gene1.symbol, gene2.symbol, "tophat"]))
                oncofuse_object = {"type": "text", "label": "No Oncofuse data available for this event"}
                if in_oncofuse_th:
                    oncofuse_object = {
                        "type": "button", "color": "lightsteelblue", "action": "window", "label": "See Oncofuse data",
                        "card": {
                            "title": "Detailed Oncofuse information",
                            "width": "100",
                            "elements": [
                                {
                                        "type": "table",
                                        "data": {
                                                "url": "/fusion_api/fusion/oncofuse/"+str(cellLine)+"/tophat/"+gene1.symbol+"|"+gene2.symbol
                                        }
                                }
                            ]
                        }}
                    
                row = []
                row.append(",\n".join([chromosome1.chromosome+":"+str(left_coord), chromosome2.chromosome+":"+str(right_coord)]))
                row.append(oncofuse_object)
                row.append(spanning_reads)
                row.append(str(spanning_mate_pairs))
                row.append(str(spanning_mate_pairs_end))
                
                rows.append(row)

    return rows

def build_ja_rows(fusions, cellLine):
    
    rows = []
    ja_fusions = []
    
    print(cellLine)
    print(len(fusions))

#     for fus in fusions:
#         for fcfusion in fus.with_fc_script:
#             if fcfusion.toCellLine[0].cell_line == cellLine:
#                 fc_fusions.append(fcfusion)
            
#     print("TOTAL " + str(len(fc_fusions)) + " FUSION CATCHER FUSIONS FOUND")

    oncofuse = get_oncofuse_info_simple()
    print(len(oncofuse))
    oncofuse_set = set()
    for item in oncofuse:
        oncofuse_set.add("#".join(item))
    
    print(next(iter(oncofuse_set)))
    
    descriptions_counter = Counter()
    tags_counter = Counter()
    
    for fusion in fusions:
        gene1 = Gene.nodes.get(symbol=fusion.gene1)
        gene2 = Gene.nodes.get(symbol=fusion.gene2)
        
        for fusion_point in fusion.fusion_points:
            print(fusion_point)
            
            chromosome1 = Chromosome.nodes.get(chromosome__exact=fusion_point.chromosome1)
            chromosome2 = Chromosome.nodes.get(chromosome__exact=fusion_point.chromosome2)
            
#             gene1 = chromosome1.of_gene[0]
#             gene2 = chromosome2.of_gene[0]
            
            fusion_point_1 = fusion_point.breakpoint1
            fusion_point_2 = fusion_point.breakpoint2
            
            for myfusion in fusion_point.jaffa_events.filter(cell_line__exact=cellLine):
                
                strand_1 = myfusion.strand1
                strand_2 = myfusion.strand2
                
                in_oncofuse_ja = "#".join([cellLine, gene1.symbol, gene2.symbol, "jaffa"]) in oncofuse_set
                oncofuse_object = {"type": "paragraph", "data": {"value": "<b>Oncofuse</b>: No Oncofuse data available for this event"}}
                if in_oncofuse_ja:
                    oncofuse_object = {
                        "type": "multi",
                        "subdata":[
                            {
                                "type": "paragraph",
                                "data": {"value": "<b>Oncofuse data</b>"}
                            },
                            { "type": "button", "color": "lightsteelblue", "action": "window", "label": "See Oncofuse data", "card": {
                                "title": "Detailed Oncofuse information",
                                "width": "100",
                                "elements": [
                                        {
                                                "type": "table",
                                                "data": {
                                                        "url": "/fusion_api/fusion/oncofuse/"+str(cellLine)+"/jaffa/"+gene1.symbol+"|"+gene2.symbol
                                                }
                                        }
                                ]
                            }}
                        ]
                     }
                
                row = []
                row.append({
                        "type": "multi",
                        "subdata": [
                            {
                                "type": "image",
                                "data": {
                                        "url": "https://d30y9cdsu7xlg0.cloudfront.net/png/43386-200.png",
                                        "width": "50px",
                                    }
                            },
                            {"type": "paragraph", "data": {"value": "<b>5' position</b>: " + "chr"+ chromosome1.chromosome+":"+str(fusion_point_1)+":"+strand_1}}
                        ]
                        })
                row.append({"type": "spacer"})
                row.append({
                        "type": "multi",
                        "subdata": [
                            {
                                "type": "image",
                                "data": {
                                        "url": "https://d30y9cdsu7xlg0.cloudfront.net/png/43386-200.png",
                                        "width": "50px",
                                    }
                            },
                            {"type": "paragraph", "data": {"value": "<b>3' position</b>: " + "chr"+ chromosome2.chromosome+":"+str(fusion_point_2)+":"+strand_2}}
                        ]
                        })
                
                row.append({"type": "spacer"})
                row.append({"type": "spacer"})
                
                row.append(oncofuse_object)
                
                row.append({"type": "paragraph", "data": {"value": "<b>Spanning pairs</b>: " + str(myfusion.spanning_pairs)}})
                row.append({"type": "paragraph", "data": {"value": "<b>Spanning pairs (unique)</b>: " + str(myfusion.spanning_reads)}})
                row.append({"type": "spacer"})
                
                row.append({
                        "type": "multi",
                        "layout": "row",
                        "align": "start center",
                        "subdata":[
                            {"type": "paragraph", "data": {"value": "<b>Inframe</b>:"}},
                            {"type": "chips", "items": [{
                                "type": "chip",
                                "value": str(myfusion.inframe),
                                "color": "white",
                                "background_color": "green" if str(myfusion.inframe) == "True" else "red"
                            }]}
                        ]
                    })
                row.append({
                        "type": "multi",
                        "layout": "row",
                        "align": "start center",
                        "subdata":[
                            {"type": "paragraph", "data": {"value": "<b>Known (Mitelman)</b>:"}},
                            {"type": "chips", "items": [{
                                "type": "chip",
                                "value": str(myfusion.known),
                                "color": "white",
                                "background_color": "green" if str(myfusion.known) == "True" else "red"
                            }]}
                        ]
                    })
                row.append({"type": "paragraph", "data": {"value": "<b>Classification</b>: " + str(myfusion.classification)}})
                
                sequence_object = {
                            "type": "paragraph", "data": {"value": "<img width='50px' src='http://www.bioscience.co.uk/userfiles/icons/Molecular%20Biology.png'><b>Fusion sequence</b>:<span>" + myfusion.sequence + "</span>"}
                        }
                row.append(sequence_object)
                row.append({"type": "spacer"})
                
                rows.append([{
                        "type": "multi",
                        "subdata": row
                    }])
    
    print("FINISHED BUILD JAFFA ROWS")
    return rows

def search_viruses2(request):
    
    return HttpResponse(json.dumps(search_viruses2_private(request)))

def search_viruses2_private(request):
    response = {}
    
    data = json.loads(request.body)
    
    c_line = data["cell_line"]
    vir = data["virus"]
    
    print(data)
    
    field_list = [
        {
            "label": "cell_line",
            "title": "Cell line",
            "tooltip": "Cell Line",
            "filters": {
                "title": "Cell line filters:",
                "list": [
                    {
                        "type": "select",
                        "key": "cell_line",
                        "title": "Select a cell line:",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]
            }
        },
        {
            "label": "disease",
            "title": "Disease",
            "tooltip": "Disease",
            "filters": {
                "title": "Disease filters:",
                "list": [
                    {
                        "type": "select",
                        "key": "disease",
                        "title": "Select a disease:",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]
            }
        },
        {
            "label": "virus",
            "title": "Virus/bacteria/phage name",
            "tooltip": "Virus involved in the fusion",
            "filters": {
                "title": "Virus filters:",
                "list": [
                    {
                        "type": "select",
                        "key": "virus",
                        "title": "Select a virus:",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]
            }
            
        },
        {
            "label": "nc",
            "title": "NC (nucleotide)",
            "tooltip": "Nucleotide identifier",
            "filters": {
                "title": "NC filters:",
                "list": [
                    {
                        "type": "select",
                        "key": "nc",
                        "title": "Select a nucleotide identifier:",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]
            }
        }
    ]
                    
    rows = []
    
    #recupero virus nella linea cellulare
    #se ho specificato solo il virus, cerca per tutte le linee cellulari
    if c_line == "ALL" and vir != "":
        for c_l in CellLine.nodes:
            
            cell_line_object = {
                "type": "text",
                "label": c_l.name,
                "color": "black",
                "url": "cell_line/" + c_l.cell_line + "/",
                "target": "_blank"
            }
            
            disease_line_object = {
                    "type": "text",
                    "label": c_l.disease_name,
                    "color": "black"
                }
            
            for virus in c_l.with_viruses:
                
                virus_object = {
                    "type": "text",
                    "label": virus.name,
                    "color": "black"
                }
                
                if vir == virus.name or vir == virus.gi or vir == virus.ref:
                    rows.append(
                        {
                            "cell_line": [cell_line_object],
                            "disease": [disease_line_object],
                            "virus": [virus_object],
                                 #{"type": "link", "data": {"label": virus.gi, "url": "https://www.ncbi.nlm.nih.gov/gene/?term="+virus.gi}},
                            "nc": [{"type": "text", "label": virus.ref, "target": "_blank", "url": "https://www.ncbi.nlm.nih.gov/nuccore/" + virus.ref}]
                        }
                    )
                    
    #se ho specificato solo la linea cellulare, cerca tutti i virus per la linea cellulare
    elif c_line != "":
        
        c_l = CellLine.nodes.get(cell_line = c_line)
        
        cell_line_object = {
                "type": "text",
                "label": c_l.name,
                "color": "black",
                "url": "cell_line/" + c_l.cell_line + "/",
                "target": "_blank"
            }
        
        disease_line_object = {
                    "type": "text",
                    "label": c_l.disease_name,
                    "color": "black"
                }
        
        for virus in c_l.with_viruses:
            
            virus_object = {
                    "type": "text",
                    "label": virus.name,
                    "color": "black"
                }
            
            if vir != "" and (vir == virus.name or vir == virus.gi or vir == virus.ref):
                rows.append(
                        {
                            "cell_line": [cell_line_object],
                            "disease": [disease_line_object],
                            "virus": [virus_object],
                                 #{"type": "link", "data": {"label": virus.gi, "url": "https://www.ncbi.nlm.nih.gov/gene/?term="+virus.gi}},
                            "nc": [{"type": "text", "label": virus.ref, "target": "_blank", "url": "https://www.ncbi.nlm.nih.gov/nuccore/" + virus.ref}]
                        }
                    )
            elif vir == "ALL":
                rows.append(
                        {
                            "cell_line": [cell_line_object],
                            "disease": [disease_line_object],
                            "virus": [virus_object],
                                 #{"type": "link", "data": {"label": virus.gi, "url": "https://www.ncbi.nlm.nih.gov/gene/?term="+virus.gi}},
                            "nc": [{"type": "text", "label": virus.ref, "target": "_blank", "url": "https://www.ncbi.nlm.nih.gov/nuccore/" + virus.ref}]
                        }
                    )
    
    total = len(rows)
    if "offset" in data:
        rows = rows[data["offset"]:data["offset"]+data["limit"]]
    response = {"structure": {"field_list": field_list}, "total": total, "hits": rows}
    
    return response
        
from itertools import chain, combinations
def generate_statistics(request):
    
#     # Chromosome
#     print("Producing statistics about Chromosome/FusionCell")
#     file =  open(os.path.dirname(__file__) + "/statistics/Chromosome_FusionCell.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Chromosome", "Number of fusion events"])
#     for chromosome in Chromosome.nodes:
#         print(chromosome)
#         five_prime_events = len(chromosome.fromChromosomeToFusionPoint)
#         three_prime_events = len(chromosome.fromFusiontoChromosome)
#         total = five_prime_events + three_prime_events
#         writer.writerow([chromosome.chromosome, total])
#     file.close()
#  
#     # Chromosome/CellLine
#     print("Producing statistics about Chromosome/CellLine")
#     file =  open(os.path.dirname(__file__) + "/statistics/Chromosome_CellLine.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Chromosome", "Number of cell lines"])
#     for chromosome in Chromosome.nodes:
#         print(chromosome)
#         cell_lines = set()
#         for gene in chromosome.of_gene:
#             for couple in gene.had:
#                 for event in couple.fusion_cells:
#                     cell_lines.add(event.cell_line)
#         writer.writerow([chromosome.chromosome, len(cell_lines)])
#     file.close()
# 
#     print("Producing statistics about Chromosome/Gene")
#     file =  open(os.path.dirname(__file__) + "/statistics/Chromosome_Gene.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Chromosome", "Number of involved genes"])
#     for chromosome in Chromosome.nodes:
#         print(chromosome)
#         genes = set()
#         for gene in chromosome.of_gene: genes.add(gene.symbol)
#         writer.writerow([chromosome.chromosome, len(genes)])
#     file.close()
# 
#     # CellLine
#     print("Producing statistics about CellLine/FusionCell")
#     file =  open(os.path.dirname(__file__) + "/statistics/CellLine_FusionCell.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["CellLine", "Number of fusion events"])
#     for cell in CellLine.nodes:
#         print(cell)
#         writer.writerow([cell.cell_line, len(cell.happen)])
#     file.close()      

    # Disease/Virus
#     print("Producing statistics about Disease/Virus")
#     file =  open(os.path.dirname(__file__) + "/statistics/Disease_Virus.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Disease", "Number of viruses"])
#     disease_viruses = {}
#     for cell in CellLine.nodes:
#         print(cell)
#         if cell.disease not in disease_viruses:
#             disease_viruses[cell.disease] = set()
#         for virus in cell.with_viruses: disease_viruses[cell.disease].add(virus.name)
#     for disease in disease_viruses:
#         writer.writerow([disease, len(disease_viruses[disease])])
#     file.close()

    # Disease/FusionCell
#     print("Producing statistics about Disease/FusionCell")
#     file =  open(os.path.dirname(__file__) + "/statistics/Disease_FusionCell.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Disease", "Number of fusion events"])
#     disease_events = {}
#     for cell in CellLine.nodes:
#         print(cell)
#         if cell.disease not in disease_events:
#             disease_events[cell.disease] = set()
#         for event in cell.happen: disease_events[cell.disease].add(event)
#     for disease in disease_events:
#         writer.writerow([disease, len(disease_events[disease])])
#     file.close()

#     # Disease/FusionCellRatio
#     file =  open(os.path.dirname(__file__) + "/statistics/Disease_FusionCellRatio.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Disease", "Number of normalized fusion events"])
#     disease_events = {}
#     counter = Counter()
#     for cell in CellLine.nodes:
#         print(cell)
#         counter[cell.disease] += 1
#         if cell.disease not in disease_events:
#             disease_events[cell.disease] = set()
#         for event in cell.happen: disease_events[cell.disease].add(event)
#     for disease in disease_events:
#         writer.writerow([disease, len(disease_events[disease])/counter[disease]])
#     file.close()

#     # Disease/CellLine
#     print("Producing statistics about Disease/CellLine")
#     file =  open(os.path.dirname(__file__) + "/statistics/Disease_CellLine.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Disease", "Number of cell lines"])
#     disease_counter = Counter()
#     for cell in CellLine.nodes:
#         print(cell)
#         disease_counter[cell.disease] += 1
#     for disease in disease_counter:
#         writer.writerow([disease, disease_counter[disease]])
#     file.close()
    
    # Distribution of predicted effect
#     print("Producing statistics about PredictedEffect/FusionCell")
#     file =  open(os.path.dirname(__file__) + "/statistics/PredictedEffect_FusionCell.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Predicted effect", "Number of gene couples"])
#     effects_counter = Counter()
#     for fusioncatcher in FusionCatcher.nodes:
#         effect1 = fusioncatcher.predicted_effect_1
#         effect2 = fusioncatcher.predicted_effect_2
#         effect = effect1
#         if effect1 != effect2:
#             effect = effect1 + "/" + effect2
#         effects_counter[effect] += 1
#     for annotation in effects_counter:
#         writer.writerow([annotation, effects_counter[annotation]])
#     file.close()
    
    # Distribution of description
#     print("Producing statistics about Description/FusionCell")
#     file =  open(os.path.dirname(__file__) + "/statistics/Description_FusionCell.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Description", "Number of fusion events", "Color"])
#     tags_of_interest = ["chimerdb2", "cosmic", "tcga", "gtex", "oncogene", "ticdb", "mitelman"] 
#     annotation_counter = Counter()
#     step = 100000
#     subsets_colors = {}
#     total = len(FusionCell.nodes) # 990627
#     for x in xrange(0, total, step):
#         print("Loading all fusion cells between " + str(x) + " and " + str(x+step))
#         for event in FusionCell.nodes[x:x+step]:
#             annotations = event.annotations
#             if annotations is None: continue
#             if len(annotations) == 0: continue
#                  
#             intersection = [x for x in annotations if x in tags_of_interest]
#             if not intersection: continue
#                  
#             for subset in chain.from_iterable(combinations(intersection,n) for n in range(1, len(intersection)+1)):
#                 subset_string = "/".join(subset)
#                 annotation_counter[subset_string] += 1
#                 descriptions_in_green_list = filter(lambda x: x in green_list, subset)
#                 descriptions_in_orange_list = filter(lambda x: x in orange_list, subset)
#                 descriptions_in_red_list = filter(lambda x: x in red_list, subset)
#                    
#                 subset_colors = Counter()
#                    
#                 if len(descriptions_in_green_list) > 0:
#                     subset_colors["green"] += len(descriptions_in_green_list)
#                 if len(descriptions_in_orange_list) > 0:
#                     subset_colors["orange"] += len(descriptions_in_orange_list)
#                 if len(descriptions_in_red_list) > 0:
#                     subset_colors["red"] += len(descriptions_in_red_list)
#                    
#                 if len(subset_colors) == 1:
#                     subset_color = subset_colors.most_common()[0][0]
#                 elif len(subset_colors) > 1:
#                     subset_color = "grey";
#    
#                 subsets_colors[subset_string] = subset_color
#                    
#     for annotation in annotation_counter:
#         writer.writerow([annotation, annotation_counter[annotation], subsets_colors[annotation]])
#     file.close()
    
    
        # Distribution of description
#     print("Producing statistics about Description/FusionCell in the CCS")
#     file =  open(os.path.dirname(__file__) + "/statistics/Description_FusionCell_CCS.csv", 'w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow(["Description", "Number of fusion events in CCS", "Color"])
#     tags_of_interest = ["chimerdb2", "cosmic", "tcga", "gtex", "oncogene", "ticdb", "mitelman"] 
#     annotation_counter = Counter()
#     step = 100000
#     subsets_colors = {}
#     total = len(FusionCell.nodes)
#     for x in xrange(0, total, step):
#         print("Loading all fusion cells between " + str(x) + " and " + str(x+step))
#         for event in FusionCell.nodes[x:x+step]:
#             if event.num_algo < 3: continue
#              
#             annotations = event.annotations
#             if annotations is None: continue
#             if len(annotations) == 0: continue
#                  
#             intersection = [x for x in annotations if x in tags_of_interest]
#             if not intersection: continue
#                  
#             for subset in chain.from_iterable(combinations(intersection,n) for n in range(1, len(intersection)+1)):
#                 subset_string = "/".join(subset)
#                 annotation_counter[subset_string] += 1
#                 descriptions_in_green_list = filter(lambda x: x in green_list, subset)
#                 descriptions_in_orange_list = filter(lambda x: x in orange_list, subset)
#                 descriptions_in_red_list = filter(lambda x: x in red_list, subset)
#                    
#                 subset_colors = Counter()
#                    
#                 if len(descriptions_in_green_list) > 0:
#                     subset_colors["green"] += len(descriptions_in_green_list)
#                 if len(descriptions_in_orange_list) > 0:
#                     subset_colors["orange"] += len(descriptions_in_orange_list)
#                 if len(descriptions_in_red_list) > 0:
#                     subset_colors["red"] += len(descriptions_in_red_list)
#                    
#                 if len(subset_colors) == 1:
#                     subset_color = subset_colors.most_common()[0][0]
#                 elif len(subset_colors) > 1:
#                     subset_color = "grey";
#    
#                 subsets_colors[subset_string] = subset_color
#                    
#     for annotation in annotation_counter:
#         writer.writerow([annotation, annotation_counter[annotation], subsets_colors[annotation]])
#     file.close()
    
    return HttpResponse("OK")

def get_descriptions_statistics(request):
    
    header = ["Description", "ES", "FC", "TF"]
    rows = [
#         ["Cis", "39,814"],
        ["inter-chromosomal", "633,256", "42,679", "14,361"],
        ["intra-chromosomal", "112,124", "48,575", "31,634"]
#         ["Read-Through", "144,447", "25,534"]
        ]
    
    return HttpResponse(json.dumps({"header": header, "items": rows}))

def get_mapping_algorithms(request):
    response = ["BOWTIE", "BOWTIE+BLAT", "BOWTIE+BLAT;BOWTIE+STAR", "BOWTIE;BOWTIE+BLAT", "BOWTIE;BOWTIE+BLAT;BOWTIE+STAR", "BOWTIE;BOWTIE+STAR", "BOWTIE+STAR"]
    
    return HttpResponse(json.dumps(response))

def get_fusion_descriptions(request):
    
    response = []
    
    labels = ["100K<gap<200K", "10K<gap<100K", "1K<gap<10K", "adjacent", "ambiguous", "antisense", "banned", "bodymap2", "cacg", "cell_lines", "cgp", "chimerdb2", "conjoing", "cosmic", "duplicates", "ensembl_partially_overlapping", "gap<1K", "gencode_fully_overlapping", "gencode_partially_overlapping", "gencode_same_strand_overlapping", "gtex", "healthy", "hpa", "known", "lincrna", "no_protein", "oncogene", "pair_pseudo_genes", "paralogs", "prostates", "pseudogene", "readthrough", "refseq_partially_overlapping", "ribosomal_protein", "rp11_gene", "rp_gene", "similar_reads", "snrna", "tcga", "ticdb", "ucsc_fully_overlapping", "ucsc_partially_overlapping", "ucsc_same_strand_overlapping"]
    
    response.append({"id": "ALL", "label": "Include any description in results"})
    
    for label in labels:
        item = {"label": label, "id": label.replace("/", "|")}
        response.append(item)
    
    return HttpResponse(json.dumps(response))

def get_predicted_effects(request):

    response = []
    
    response.append({"id": "ALL", "label": "Include any predicted effect in results"})
    
    labels = ["CDS(truncated)/intronic", "CDS(truncated)/UTR", "exonic(no-known-CDS)/CDS(truncated)", "in-frame", "intronic/CDS(truncated)", "intronic/UTR", "out-of-frame", "UTR/CDS(truncated)", "UTR/exonic(no-known-CDS)", "UTR/intergenic", "UTR/intronic", "UTR/UTR"]
    
    for label in labels:
        item = {"label": label, "id": label.replace("/", "|")}
        response.append(item)
    
    return HttpResponse(json.dumps(response))

def get_algorithms(request):
    
    response = []
    
    labels = ["FusionCatcher", "Tophat-Fusion", "EricScript", "Jaffa"]
    
    for label in labels:
        item = {"label": label, "id": label}
        response.append(item)
    
    return HttpResponse(json.dumps(response))

def get_fusions2(request):
    return HttpResponse(json.dumps(get_fusions2_private(request)))

def get_fusions2_private(request):
    
    response = {}
    
    data = json.loads(request.body)
    
    c_line = data["cell_line"]
    fusion_descriptions_raw = data["fusion_descriptions_raw"]
    predicted_effects_raw = data["predicted_effects_raw"]
    num_algorithms = data["num_algorithms"] if "num_algorithms" in data else 1

    if fusion_descriptions_raw == "ALL" and predicted_effects_raw == "ALL":
        return search_for_cell_line2_private(request, True)
    
    print("get_fusions2_private", data)
    
    fusions = []
    
    fusion_descriptions = fusion_descriptions_raw
#     if not fusion_descriptions or fusion_descriptions[0] == "ALL": fusion_descriptions = []
    predicted_effects = predicted_effects_raw
#     if not predicted_effects or predicted_effects[0] == "ALL": predicted_effects = []
    
    predicted_effects_searched = "(fc.predicted_effect_1 = '" + predicted_effects+"' OR fc.predicted_effect_2 = '"+predicted_effects+"')" if predicted_effects != "ALL" else ""
    cell_lines_searched = "f.cell_line = " +"'"+c_line+"'" if c_line != "ALL" else ""
    # fusion_descriptions_searched = "'" + fusion_descriptions + "'"+" IN "+"fc.description" if fusion_descriptions != "ALL" else ""
    fusion_descriptions_searched = "'" + fusion_descriptions + "'"+" IN "+"f.annotations" if fusion_descriptions != "ALL" else ""
    at_least = "f.num_algo >= " + str(num_algorithms)

    novel = data["novel"] if "novel" in data else False
    if novel == "strict":
        novel_condition = "NOT EXISTS(g.cosmic)"
        novel_condition += " AND f.quality = 'N/A' AND (NOT EXISTS(f.annotations) OR (" + " AND ".join(["NOT '"+x+"' IN f.annotations" for x in annotations.keys()]) + "))"
#         novel_path = ",(g:GeneCouple)--(f)"
    elif novel == "relaxed":
        novel_condition = "NOT EXISTS(g.cosmic)"
#         novel_condition += " AND (f.quality = 'N/A' OR f.quality = 'grey') AND (NOT EXISTS(f.annotations) OR (" + " AND ".join(["NOT '"+x+"' IN f.annotations" for x in annotations.keys()]) + "))"
#         novel_path = ",(g:GeneCouple)--(f)"
    else:
        novel_condition = ""
#         novel_path = ""
    
    conditions = filter(lambda x: x != "", [novel_condition, at_least, cell_lines_searched, predicted_effects_searched, fusion_descriptions_searched])
    where_condition = "" if not conditions else " WHERE " + (" AND ".join(conditions))

    if predicted_effects != "ALL":
        path = "(fc:FusionCatcher)--(FusionPoint)--(f:FusionCell)--(g:GeneCouple)"
    else:
        path = "(f:FusionCell)--(g:GeneCouple)"
    
    count_query = "MATCH " + path + " " + where_condition + " return count(distinct(f))"
    print(count_query)
    r, m = db.cypher_query(count_query)
    total = r[0][0]
    
    query = "MATCH " + path + " " + where_condition + " return distinct(f) ORDER BY f.num_algo DESC " + ("SKIP "+ str(data["offset"])+ " LIMIT " + str(data["limit"]) if "offset" in data else "")
    print(query)
    results, meta = db.cypher_query(query)
    fusions = [FusionCell.inflate(row[0]) for row in results]
    for fusion in fusions:
        print(fusion)
    
    print("BUILDING ROWS FOR " + str(len(fusions)) + " *** Taking in consideration only FC ***")

    rows = build_rows2(fusions, c_line, True)
    header = get_header2()[:-3]
    
    response = {"structure": {"field_list": header}, "total": total, "hits": rows}
    
    return response

def search_by_algorithm(request):
    
    return HttpResponse(json.dumps(search_by_algorithm_private(request)))

def search_by_algorithm_private(request):
    
    response = {}
    
    data = json.loads(request.body)
    print(data)
    
    if "algorithms" not in data or not data["algorithms"]:
        return search_for_cell_line2_private(request)
    
    algorithms = data["algorithms"]
    c_line = data["cell_line"]
    strict = data["strict"] if "strict" in data else False
    
    field_list = get_header2()
    
    novel = data["novel"] if "novel" in data else False
    if novel == "strict":
        novel_condition = "NOT EXISTS(g.cosmic) AND f.quality = 'N/A' AND (NOT EXISTS(f.annotations) OR (" + " AND ".join(["NOT '"+x+"' IN f.annotations" for x in annotations.keys()]) + ")) AND"
        novel_path = "(g:GeneCouple)--"
    elif novel == "relaxed":
        novel_condition = "NOT EXISTS(g.cosmic) AND (f.quality = 'N/A' OR f.quality = 'grey') AND (NOT EXISTS(f.annotations) OR (" + " AND ".join(["NOT '"+x+"' IN f.annotations" for x in annotations.keys()]) + ")) AND"
        novel_path = "(g:GeneCouple)--"
    else:
        novel_condition = ""
        novel_path = ""
    
    where_conditions = []
    all_algorithms = ["FC", "ES", "TH", "JA"]
    for algo in all_algorithms:
        algorithm_name = "FusionCatcher" if algo == "FC" else "EricScript" if algo == "ES" else "Tophat-Fusion" if algo == "TH" else "Jaffa" if algo == "JA" else "Unknown"
        print(algorithms, algorithm_name, algorithm_name in algorithms)
        
        if algorithm_name in algorithms:            
            where_conditions.append(" '" + algo + "'" + " in f.algos ")
        else:
            if strict:
                where_conditions.append(" (NOT '" + algo + "'" + " in f.algos) ")
        
    if c_line != "ALL":
        where_conditions.append(" f.cell_line= '" + c_line + "' ")
        
    count_query = "MATCH "+novel_path+"(f:FusionCell) WHERE "+ novel_condition + " "+ "AND".join(where_conditions)+" return count(distinct(f))"
    print("COUNT QUERY", count_query)
    r, m = db.cypher_query(count_query)
    total = r[0][0]
    print(total)
    
    query = "MATCH "+novel_path+"(f:FusionCell) WHERE "+ novel_condition + " "+ "AND".join(where_conditions)+" return distinct(f) ORDER BY f.num_algo DESC " + ("SKIP "+ str(data["offset"])+ " LIMIT " + str(data["limit"]) if "offset" in data else "")
    print("QUERY", query)
    results, meta = db.cypher_query(query)
    fusions = [FusionCell.inflate(row[0]) for row in results]

#     query = FusionCell.nodes.filter(num_algo__exact=len(algorithms))
#     for algo in algorithms:
#         print(algo, type(algo))
#         query = query.filter(algos__contains=algo)
    
#     c_line = data["cell_line"] if "cell_line" in data else "ALL"
#     cell_lines_searched = "c.cell_line = " +"'"+c_line+"'" if c_line != "ALL" else ""
#     algorithms_searched = []
#     for i,algorithm in enumerate(algorithms):
#         algorithms_searched.append("(f)--(fp"+str(i)+":FusionPoint)--(algo"+str(i)+":"+"Algorithm)--" + "("+algorithm+")")
#     
#     conditions = filter(lambda x: x != "", [cell_lines_searched])
#     where_condition = "" if not conditions else " WHERE " + (" AND ".join(conditions))
#     
#     count_query = "MATCH "+",".join(algorithms_searched)+ ",(f:FusionCell)--(c:CellLine)" + where_condition + " return count(distinct(f))"
#     print(count_query)
#     r, m = db.cypher_query(count_query)
#     total = r[0][0]
#     
#     query = "MATCH "+",".join(algorithms_searched)+ ",(f:FusionCell)--(c:CellLine)" + where_condition + " return distinct(f) ORDER BY f.num_algo DESC " + ("SKIP "+ str(data["offset"])+ " LIMIT " + str(data["limit"]) if "offset" in data else "")
#     print(query)
#     results, meta = db.cypher_query(query)
#     fusions = [FusionCell.inflate(row[0]) for row in results]
    
    print("BUILDING ROWS FOR " + str(len(fusions)))

    rows = build_rows2(fusions, c_line)

    response = {"structure": {"field_list": field_list}, "total": total, "hits": rows}
    
    #if c_line != "ALL": response = filter_table(response, c_line)
    
    return response


# def gen_stat_file(pairs):
#     ids = list(pairs.keys())[0]
#     paths = list(pairs.values())
#     paths = paths[0]
#     #print(ids)
#     #print(paths)
#     startnode = paths[0][0]
#     #print(startnode)
#     file =  open(ids[0][0]+'_'+ids[1][0]+'.csv','w')
#     writer = csv.writer(file, lineterminator='\n')
#     writer.writerow([ids[0][0],ids[1][0]])
#     for x in startnode.nodes.all():
#         print(x)
#         nodes = {}
#         for path in paths:
#             #print(path)
#             #print(len(path))
#             query = "match (" + path[0].__name__ + "{"+ids[0][1]+":'"+eval("x."+ids[0][1])+"'})-"
#             for i in range(1,len(path)-2,2):
#                 query = query + "[:" +  path[i] + "]-(" + path[i+1].__name__ + ")-" 
#                 #print(path[i+1])
#             query = query + "[:" + path[len(path)-2] + "]-(a:" + path[(len(path)-1)].__name__+") return distinct a"
#             #print(query)
#             if(db.cypher_query(query)[0]):
#                 for row in db.cypher_query(query)[0]:
#                     #print(row[0].properties[eval("'"+ids[1][1]+"'")])
#                     if row[0].properties[eval("'"+ids[1][1]+"'")] not in nodes:
#                         #print(row[1].properties[eval("'"+ids[1][1]+"'")])
#                         nodes[row[0].properties[eval("'"+ids[1][1]+"'")]] = row[0].properties[eval("'"+ids[1][1]+"'")]
#         #print(len(nodes))
#         writer.writerow([eval("x."+ids[0][1]),len(nodes)])
#     file.close()
    
def get_header():
    return ["Cell line",
            "Gene 1",
            "Gene 2",
            "FusionCatcher",
            "EricScript",
            "Tophat-Fusion",
            "Cosmic",
            #"Annotation"
            #"In Oncofuse",
            ]
    
def get_fc_header():
    return [
        {"type": "text", "label": "Chromosome:breakpoint:strand", "tooltip": "It is the pair of chromosomal position of the 5' and 3' end of fusion junction (chromosome:position:strand); 1-based coordinate."},
        {"type": "text", "label": "Predicted effect", "tooltip": "Predicted effect of the candidate fusion gene using the annotation from Ensembl database."},
        {"type": "text", "label": "Oncofuse", "tooltip": "An indicator that shows whether this event is supported also by Oncofuse."},
        {"type": "text", "label": "Description", "tooltip": "Type of the fusion gene."},
        {"type": "text", "label": "Counts of common mapping reads", "tooltip": "Count of reads mapping simultaneously on both genes which form the fusion gene."},
        {"type": "text", "label": "Spanning pairs (of which unique)", "tooltip": "Count of pairs of reads supporting the fusion (unique reads shown in parentheses)."},
        {"type": "text", "label": "Longest anchor found", "tooltip": "Longest anchor (hangover) found among the unique reads mapping on the fusion junction."},
        {"type": "text", "label": "Fusion finding method(s)", "tooltip": "Aligning method used for mapping the reads and finding the fusion genes."},
        {"type": "text", "label": "Fusion sequence", "tooltip": "The inferred fusion junction (the asterisk sign marks the junction point)."},
        {"type": "text", "label": "Predicted fused transcript(s) and protein(s)", "tooltip": "All possible known fused transcripts and predicted amino acid sequences of all possible fused proteins."}
        ]
    
def get_es_header():
    return [
        #"Gene pair symbols",
        {"type": "text", "label": "Chromosome:breakpoint:strand", "tooltip": "It is the pair of chromosomal position of the 5' and 3' end of fusion junction (chromosome:position:strand); 1-based coordinate."},
        {"type": "text", "label": "Oncofuse", "tooltip": "An indicator that shows whether this event is supported also by Oncofuse."},
        {"type": "text", "label": "Crossing reads", "tooltip": "The number of paired end discordant reads."},
        {"type": "text", "label": "Spanning reads", "tooltip": "The number of paired end reads spanning the junction."},
        {"type": "text", "label": "Mean insert size", "tooltip": "Mean of insert sizes of crossing + spanning reads."},
        {"type": "text", "label": "Homology", "tooltip": "If filled, all the homologies between the fusion junction and Ensembl genes."},
        {"type": "text", "label": "Fusion type", "tooltip": "Intra-chromosomal, inter-chromosomal, read-through or CIS."},
        {"type": "text", "label": "Junction sequence", "tooltip": "Predicted junction fusion sequence."},
        {"type": "text", "label": "Gene expr 1", "tooltip": "Read count based estimation of the expression level of 5' gene."},
        {"type": "text", "label": "Gene expr 2", "tooltip": "Read count based estimation of the expression level of 3' gene."},
        {"type": "text", "label": "Gene expr fused", "tooltip": "Read count based estimation of the expression level of the predicted chimeric transcript."},
        {"type": "text", "label": "Es", "tooltip": "Edge score."},
        {"type": "text", "label": "GJS", "tooltip": "Genuine Junction score."},
        {"type": "text", "label": "US", "tooltip": "Uniformity score."},
        {"type": "text", "label": "EricScore", "tooltip": "EricScore score (adaboost classifier)."}
        ]
    
def get_tophat_header():
    return [
        #"Gene pair symbols",
            {"type": "text", "label": "Chromosome:breakpoint", "tooltip": "It is the pair of chromosomal position of the 5' and 3' end of fusion junction (chromosome:position); 1-based coordinate."},
            {"type": "text", "label": "Oncofuse", "tooltip": "An indicator that shows whether this event is supported also by Oncofuse."},
            {"type": "text", "label": "Spanning reads", "tooltip": "The number of paired end reads spanning the junction."},
            {"type": "text", "label": "Mate pairs", "tooltip": "Mate pairs that support the fusion."},
            {"type": "text", "label": "Spanning mate pairs", "tooltip": "Mate pairs that support the fusion whose one end spans the fusion."}
            ]

def get_gold_pairs():
    header = ["Gene1", "Gene2", "ID_Gene1", "ID_Gene2", "url"]
    rows = []
    
    txt_file = open(os.path.dirname(__file__) + "/data/fusion_pair_gold.txt", "r")
    for line in txt_file:
        line = line.rstrip()
        rows.append(line.split("\t"))
    
    return {"header": header, "items": rows}


# def get_gold_genes():
#     header = ["Gene", "url"]
#     rows = []
#     
#     txt_file = open(os.path.dirname(__file__) + "/data/genes_gold.txt", "r")
#     for line in txt_file:
#         line = line.rstrip()
#         rows.append(line.split("\t"))
#     
#     return {"header": header, "items": rows}

def get_oncofuse_info_simple():
    
    items = []
    for line in open(os.path.dirname(__file__) + "/data/downloads/oncofuse/oncofuse_info.txt"):
        line = line.rstrip()
        items.append(line.split("\t"))
    return items

def get_oncofuse_info(cell_line_list = [], algorithm = "", gene_list = []):
    
    header = []
    rows = []
        
    print(cell_line_list, algorithm, gene_list)
    
    ccle_map = {}
    ccle_infos = get_ccle_infos()
    for cell_info in ccle_infos:
        ccle_map[cell_info["ID"]] = cell_info["Cell Line"]
    
    for cell_line in cell_line_list:
        print(cell_line)
        regex = os.path.dirname(__file__) + "/data/downloads/oncofuse/*/*/" + "*"+cell_line+".*" + algorithm +".out"
        print(regex)
        for file in glob.glob(regex):
            print(file)
            txt_file = open(file, "r")
            line_no = 0
            for line in txt_file:
                
                line_no += 1
                
                fields = line.rstrip().split("\t")
                
                for i,field in enumerate(fields):
                    if ";" in field:
                        fields[i] = field.split(";")
                
                if line_no == 1:
                    fields[0] = "Cell line"
                    header = [field.replace("_", " ") for field in fields]
                else:
                    if gene_list is [] or gene_list[0] == fields[6] and gene_list[1] == fields[13]:
                        # fields[0] = "|".join([x for x in fields[0].split("/") if x.startswith("CCLE")])
                        start = fields[0].find("CCLE")
                        if start == -1:
                            print(file, fields[0])
                        end = fields[0].find("/", start)
                        if end == -1:
                            end = fields[0].find(".", start)                        
                            if end == -1:
                                end = len(fields[0])                    
                        
                        fields[0] = fields[0][start:end]
                        fields[0] = ccle_map[fields[0]]
                        
                        rows.append(fields)
            
    response = {"header": header, "items": rows}
        
    return response

def get_oncofuse(request, cell_line, algorithm = "", gene_list = ""):
    return HttpResponse(json.dumps(get_oncofuse_info(cell_line.split("|"), algorithm, gene_list.split("|"))))

def get_distribution_virus_by_disease(request):
    header = ["Disease", "Number of detected viruses"]
    rows = []
    
    ccle_map = {}
    disease_map = {}
    ccle_infos = get_ccle_infos()
    for cell_info in ccle_infos:
        ccle_map[cell_info["ID"]] = cell_info
        disease_map[cell_info["Disease"]] = {
            "label": cell_info["Disease"],
            "name": "Disease: " + cell_info["Disease name"] + " ("+cell_info["Disease"]+")",
            "value": "Number of detected viruses for this disease"
            }
    
    counter = Counter()
    for virus in Virus.nodes:
        diseases = set()
        for cell_line in virus.fromCellLineToVirus:
            cell_line_info = ccle_map[cell_line.cell_line]
            disease = cell_line_info["Disease"]
            diseases.add(disease)
            
        for disease in diseases:
            counter[disease] += 1
            
    for disease in counter.most_common():
        rows.append([disease_map[disease[0]], disease[1]])
        
    return HttpResponse(json.dumps({"header": header, "items": rows, "descriptions": header}))

def get_fusion_by_disease(request):
    header = ["Disease", "No. of fusion events for disease"]
    rows = []
    
    ccle_map = {}
    disease_map = {}
    ccle_infos = get_ccle_infos()
    for cell_info in ccle_infos:
        ccle_map[cell_info["ID"]] = cell_info
        disease_map[cell_info["Disease"]] = {
            "id": cell_info["Disease"],
            "label": cell_info["Disease"],
            "name": "Disease: " + cell_info["Disease name"],
            "value": "Normalized fusion events per disease"
        }
        
    disease2cells = Counter()
    counter = Counter()
    for cell_line in CellLine.nodes:
        cell_line_info = ccle_map[cell_line.cell_line]
        disease = cell_line_info["Disease"]
        counter[disease] += len(cell_line.happen)
        disease2cells[disease] += 1
    
#     print(counter)
#     print(disease2cells)
    for disease in counter:
        counter[disease] = counter[disease] / disease2cells[disease]
#     print(counter)
            
    for disease in counter.most_common():
        rows.append([disease_map[disease[0]], disease[1]])
        
    return HttpResponse(json.dumps({"header": header, "items": rows}))

def get_event_map(request):
    header = ["Type of event", "No. of fusion events"]
    rows = []
#     rows = [
#             [{"id": 0, "label": "False positives"}, 299],
#             [{"id": 1, "label": "True positives"}, 531],
#             [{"id": 2, "label": "False positives with medium probability"}, 32],
#             [{"id": 3, "label": "Novel candidates with ambiguous classification"}, 652],
#             [{"id": 4, "label": "Novel candidates"}, 1376]
#         ]

    chart = Counter()
    for fusion in FusionCell.nodes.filter(num_algo__exact = 3):
        if fusion.gene_couple[0].cosmic: continue
        
        chart[fusion.quality] += 1
    
    colors = ["#CC0000", "#6AA84F", "#FF9900", "#CCCCCC", "#990099"]
    
    quality_map = {
        'red': {"id": 0, "quality": "red", "label": "False positives"},
        'green': {"id": 1, "quality": "green","label": "True positives"},
        'orange': {"id": 2, "quality": "orange","label": "False positives with medium probability"},
        'grey': {"id": 3, "quality": "grey","label": "Novel candidates with ambiguous classification"},
        'N/A': {"id": 4, "quality": "N/A","label": "Novel candidates"}
    }
    
    for quality in sorted(quality_map.values(), key=lambda x: x["id"]):
        rows.append([quality, chart[quality["quality"]]])
    
    return HttpResponse(json.dumps({"header": header, "items": rows, "colors": colors}))

def cell_line_info(request, cell_line):
    header = []
    rows = []
    
    cell = CellLine.nodes.get(cell_line=cell_line)
    info = {
        "type": "multi",
        "layout": "row",
        "layout_align": "start center",
        "subdata":[
            {
                "type": "image",
                "data": {
                    "url": "cell-icon.png",
                    "width": "35px",
                    "margin": "10px"
                }
            },{
                "type": "paragraph",
                "data": {"value": "<b>Name</b>: " + cell.name}
            }]
        }
    rows.append([info])
    rows.append([{
        "type": "multi",
        "layout": "row",
        "layout_align": "start center",
        "subdata":[
            {
                "type": "image",
                "data": {
                    "url": "disease-icon.png",
                    "width": "35px",
                    "margin": "10px"
                }
            },{
                "type": "paragraph",
                "data": {"value": "<b>Disease</b>: " + cell.disease_name + " ("+cell.disease+")"}
            }]
          }])
    
    if cell.with_viruses:
        
        viruses = []
        for virus in cell.with_viruses.match(count_of_mapping_reads__gte=5):
            mapping_reads = virus.fromCellLineToVirus.relationship(cell).count_of_mapping_reads
            viruses.append({
                "type": "chip",
                "value": virus.name,
                "background_color": "grey",
                "color": "black",
                "number": mapping_reads,
                "title": "Number of supporting reads: " + str(mapping_reads)
            })
        viruses.sort(key=lambda x: x["number"], reverse=True)
        
        rows.append([{
        "type": "multi",
        "layout": "row",
        "layout_align": "start center",
        "subdata":[
            {
                "type": "image",
                "data": {
                    "url": "virus-icon.png",
                    "width": "35px",
                    "margin": "10px"
                }
            },{
                "type": "paragraph",
                "data": {"value": "<b>Virus(es)</b><br/>(" + str(len(viruses)) + " total)"}
            }, {"type": "chips", "limit": 5, "items": viruses}]
          }])
        
    return HttpResponse(json.dumps({"header": header, "items": rows}))

def cell_line_minor_info(request, cell_line):
    header = []
    rows = []
    
    cell = CellLine.nodes.get(cell_line=cell_line)
    infos = get_ccle_infos()
    ccle_info = None
    
    for item in infos:
        if item["Cell Line"] == cell.name:
            ccle_info = item
            
    print("INFO", ccle_info)
    ccle_info_minor = {}
    for key, value in ccle_info.iteritems():
        if key in ["Disease name", "Cell Line", "ID", "Disease"]: continue
        ccle_info_minor[key] = value
    
    print("MINOR", ccle_info_minor)
    
    header = []
    for key, value in ccle_info_minor.iteritems():
        image_name = None
        
        if key == "COSMIC identifier": image_name = "cosmic_logo.png"
        if key == "Methylation": image_name = "methylation-icon.png"; continue
        if key == "Site": image_name = "body/" + value + ".png"
        if key == "Drug": image_name = "drug-icon.png"; continue
        if key == "Histology": image_name = "histology-icon.png"
        if key == "Disease link": image_name = "disease-icon.png"
        if key == "Copy Number Alterations (CNA)": image_name = "cna-icon.png"; continue
        if key == "Whole Exome Sequencing (WES)": image_name = "wes-icon.png"; continue
        if key == "Gene Expression": image_name = "gene-icon.png"; continue
        if key == "Drug resistance": image_name = "drug-resistance-icon.png"
        
        if value == "Y": value = "Yes"
        if value == "N": value = "No"
        
        value = value.replace("_", " ").capitalize()
        
        if key == "Disease link":
            value = "<a target='_blank' href='"+value+"'>"+ccle_info["Disease name"]+"</a>"
            key = "Disease"
        if key == "COSMIC identifier": value = "<a target='_blank' href='http://cancer.sanger.ac.uk/cell_lines/sample/overview?id="+value+"'>"+value+"</a>"
        if key == "Drug resistance": value = "<a target='_blank' href='http://www.cancerrxgene.org/translation/CellLine/"+value+"'>"+"Explore"+"</a>"
        
        if image_name is None: image = {"type": "text", "label": "-"}
        else: image = {"type": "image", "data": {"width": "50px", "url": image_name}}
        
        if key == "Site": image = {"type": "image", "data": {"width": "75px", "url": image_name}}
        
        rows.append([
            image, 
            {"type": "text", "label": "<b>"+key+"</b>"},
            {"type": "text", "text_align": "left", "label": value}])
    
    if not rows:
        rows.append(["No additional information available for this cell line."])
    
    return HttpResponse(json.dumps({"header": header, "items": rows}))

def gene_couples(request, cell_line):
    start = datetime.datetime.now()
    
    counter = Counter()
    rows = []
    header = ["Pair of genes", "No. of fusion events"]
    fusion_cells = FusionCell.nodes.filter(cell_line = cell_line)
    for fusion in fusion_cells:
        if fusion.num_algo < 2: continue
        #if fusion.annotations is not None and len(filter(lambda x: x in red_list, fusion.annotations)) > 0: continue
        
        if len(counter) % 1000 == 0: print(len(counter))
        
        gene1 = fusion.gene1
        gene2 = fusion.gene2
        
        counter[gene1 + "#" + gene2] += fusion.num_algo
        
    end = datetime.datetime.now()
    
    total_time = end - start
    
    print("Building rows...")
    for pair in counter.most_common():
        rows.append([pair[0], pair[1]])
    
    return HttpResponse(json.dumps({"header": header, "items": rows}))

def cell_line_chromosome_info(request, cell_line):
    start = datetime.datetime.now()
    
    counter = Counter()
    rows = []
    header = ["Chromosomes", "No. of fusion events"]
    step = 100000
    fusion_cells = FusionCell.nodes.filter(cell_line = cell_line)
    for fusion in fusion_cells:
        #if fusion.num_algo < 3: continue
        #if fusion.annotations is not None and len(filter(lambda x: x in red_list, fusion.annotations)) > 0: continue
        
        if len(counter) % 1000 == 0: print(len(counter))
        
        for fusion_point in fusion.fusion_points:
            counter["chr"+fusion_point.chromosome1] += 1
            counter["chr"+fusion_point.chromosome2] += 1
        
    end = datetime.datetime.now()
    
    total_time = end - start
    
    print("Building rows...")
    for pair in counter.most_common():
        rows.append([{"title": "Chromosome " + pair[0].replace("chr", ""), "value": "Number of gene fusion events", "label": pair[0]}, pair[1]])
    
    return HttpResponse(json.dumps({"header": header, "items": rows, "descriptions": header}))


def fusion_event_info(request, fusion_id):
    
    header = []
    rows = []
    
    fusion_event = FusionCell.nodes.get(fusion_cell_id=fusion_id)
    cell = CellLine.nodes.get(cell_line = fusion_event.cell_line)
    print(fusion_event)
    
    infos = get_ccle_infos()
    tissue = "None"
    cell_info = None
    for el in infos:
        if el["ID"] == cell.cell_line:
            cell_info = el
            break
        
    header = []
    rows = []
    
    if "Site" in el:
        rows.append([{
            "type": "multi",
            "layout": "row",
            "layout_align": "start center",
            "subdata":[
                {
                    "type": "image",
                    "data": {
                        "url": "human-icon.png",
                        "width": "35px",
                        "margin": "10px"
                    }
                },{
                    "type": "paragraph",
                    "data": {"value": "<b>Primary tissue</b>: " + el["Site"].capitalize().replace("_", " ")}
                }]
            }])
    
    rows.append([{
        "type": "multi",
        "layout": "row",
        "layout_align": "start center",
        "subdata":[
            {
                "type": "image",
                "data": {
                    "url": "cell-icon.png",
                    "width": "35px",
                    "margin": "10px"
                }
            },{
                "type": "paragraph",
                "data": {"value": "<b>Involved cell line</b>: " + "<a href='cell_line/"+cell.cell_line+"'>"+cell.name+"</a>"}
            }]
        }])
    
    disease_object = el["Disease name"] + " ("+el["Disease"]+")"
    if el["Disease link"] is not None:
        disease_object = "<a target='_blank' href='"+ el["Disease link"] +"'>" + disease_object + "</a>"
    
    rows.append([{
        "type": "multi",
        "layout": "row",
        "layout_align": "start center",
        "subdata":[
            {
                "type": "image",
                "data": {
                    "url": "disease-icon.png",
                    "width": "35px",
                    "margin": "10px"
                }
            },{
                "type": "paragraph",
                "data": {"value": "<b>Involved disease</b>: " + disease_object}
            }]
        }])
    
    return HttpResponse(json.dumps({"header": header, "items": rows}))

def fusion_event_minor_info(request, fusion_id):
    
    header = []
    rows = []
    
    fusion_event = FusionCell.nodes.get(fusion_cell_id=fusion_id)
    cell = CellLine.nodes.get(cell_line = fusion_event.cell_line)
    print(fusion_event)
    
    header = []
    rows = []
    
    algocode2string = {"ES": "EricScript", "FC": "FusionCatcher", "TH": "Tophat-Fusion"}
    algo_strings = [algocode2string[x] for x in fusion_event.algos]
    
    fusioncatcher_button = {"type": "text", "data": {"value": "-"}}
    if "FC" in fusion_event.algos:
        fusioncatcher_button = {
                    "type": "button",
                    "color": "lightblue",
                    "action": "window",
                    "title": "This fusion event is supported by FusionCatcher",
                    "label": "FC",
                    "items": [{
                        "type": "text",
                        "label": "FC"
                    }],
                    "card": {
                        "title": "Detailed information given by FusionCatcher for gene pair: " + str(fusion_event.gene1) + "/" + str(fusion_event.gene2),
                        "width": "100",
                        "elements": [
                                {
                                        "type": "table",
                                        "data": {
                                                "url": "/fusion_api/fusion/fusioncatcher/"+str(fusion_event.fusion_cell_id) + "/" + str(fusion_event.cell_line)
                                        }
                                }
                        ]
                        }
                }
        
    ericscript_button = {"type": "text", "data": {"value": "-"}}
    if "ES" in fusion_event.algos:
        ericscript_button = {
                    "type": "button",
                    "color": "lightgreen",
                    "action": "window",
                    "title": "This fusion event is supported by EricScript",
                    "label": "ES",
                    "items": [{
                        "type": "text",
                        "label": "ES"
                    }],
                    "card": {
                        "title": "Detailed information given by EricScript for gene pair: " + str(fusion_event.gene1) + "/" + str(fusion_event.gene2),
                        "width": "100",
                        "elements": [
                                {
                                        "type": "table",
                                        "data": {
                                                "url": "/fusion_api/fusion/ericscript/"+str(fusion_event.fusion_cell_id) + "/" + str(fusion_event.cell_line)
                                        }
                                }
                        ]
                        }
                }
    
    tophat_button = {"type": "text", "data": {"value": "-"}}
    if "TH" in fusion_event.algos:
        tophat_button = {
                    "type": "button",
                    "color": "lightcoral",
                    "action": "window",
                    "title": "This fusion event is supported by Tophat-Fusion",
                    "label": "TH",
                    "items": [{
                        "type": "text",
                        "label": "TH"
                    }],
                    "card": {
                        "title": "Detailed information given by Tophat-Fusion for gene pair: " + str(fusion_event.gene1) + "/" + str(fusion_event.gene2),
                        "width": "100",
                        "elements": [
                                {
                                        "type": "table",
                                        "data": {
                                                "url": "/fusion_api/fusion/tophat/"+str(fusion_event.fusion_cell_id) + "/" + str(fusion_event.cell_line)
                                        }
                                }
                        ]
                        }
                }
    
    rows.append([{
        "type": "multi",
        "layout": "row",
        "layout_align": "start center",
        "subdata":[
            {
                "type": "image",
                "data": {
                    "url": "algorithm-icon.png",
                    "width": "35px",
                    "margin": "10px"
                }
            },{
                "type": "paragraph",
                "data": {"value": "<b>Number of supporting algorithms</b>: " + str(fusion_event.num_algo) + " ("+ ", ".join(algo_strings) +")"}
            }, fusioncatcher_button, ericscript_button, tophat_button]
        }])
    
    if fusion_event.annotations:
        rows.append([{
            "type": "multi",
            "layout": "row",
            "layout_align": "start center",
            "subdata":[
                {
                    "type": "image",
                    "data": {
                        "url": "information-icon.png",
                        "width": "35px",
                        "margin": "10px"
                    }
                },{
                    "type": "paragraph",
                    "data": {"value": "<b>Annotations</b>: " + str(len(fusion_event.annotations)) + " ("+ cgi.escape(", ".join(fusion_event.annotations)) +")"}
                }]
              }])
    
    return HttpResponse(json.dumps({"header": header, "items": rows}))

def fusion_event_link(request, fusion_id):
    
    header = []
    rows = []
    
    fusion_event = FusionCell.nodes.get(fusion_cell_id=fusion_id)
    event_annotations = fusion_event.annotations or []
    
    header = []
    rows = []
    
    gene_couple = fusion_event.gene_couple[0]
    
    if gene_couple.cosmic:
        rows.append([{
            "type": "multi",
            "layout": "row",
            "layout_align": "start center",
            "subdata":[
                {
                    "type": "image",
                    "data": {
                        "url": "cosmic_logo.png",
                        "width": "35px",
                        "margin": "10px"
                    }
                },{
                    "type": "paragraph",
                    "data": {"value": "<a target='_blank' href='"+str(gene_couple.cosmic)+"'><b>Cosmic</b></a>"}
                }]
              }])
    
    for annotation in event_annotations:
        if annotation in annotations:
            annotation_info = annotations[annotation]
            image = annotation_info["img_src"]
            width = "35px"
            
            if annotation in ["chimerdb", "known", "tcga"]: width = "80px"
             
            url = annotation_info["url"] if "url" in annotation_info else None
            
            if annotation == "ticdb":           
                    annotation_row = {
                        "type": "multi",
                        "layout": "row",
                        "layout_align": "start center",
                        "subdata":[
                            {
                                "type": "image",
                                "data": {
                                    "url": image,
                                    "width": width,
                                    "margin": "10px"
                                }
                            },{
                                "type": "paragraph",
                                "data": {"value": "<a target='_blank' href='"+url+"results.php?hgnc="+fusion_event.gene1+"'>"+"<b>" + fusion_event.gene1 + "</b>"+"</a>"}
                            },{
                                "type": "paragraph",
                                "data": {"value": ",&nbsp;"}
                            },{
                                "type": "paragraph",
                                "data": {"value": "<a target='_blank' href='"+url+"results.php?hgnc="+fusion_event.gene2+"'>"+"<b>" + fusion_event.gene2 + "</b>"+"</a>"}
                            }]
                      }
                    
            else:
                if annotation == "known":
                    annotation = "pubmed"
                    url = "https://www.ncbi.nlm.nih.gov/pubmed/?term=("+fusion_event.gene1+"[Title%2FAbstract])+AND+("+fusion_event.gene2+"[Title%2FAbstract])+AND+(gene%20fusion[Title%2FAbstract])"
                
                label = annotation.capitalize().replace("_", " ")
                
                value = "<b>" + label + "</b>"
                if url: value = "<a target='_blank' href='"+url+"'>"+value+"</a>"
    
                annotation_row = {
                    "type": "multi",
                    "layout": "row",
                    "layout_align": "start center",
                    "subdata":[
                        {
                            "type": "image",
                            "data": {
                                "url": image,
                                "width": width,
                                "margin": "10px"
                            }
                        },{
                            "type": "paragraph",
                            "data": {"value": value}
                        }]
                  }
                
            rows.append([annotation_row])
            
    if not rows:
        rows.append([
            {
                "type": "multi",
                "layout": "row",
                "layout_align": "start center",
                "subdata":[
                    {
                        "type": "image",
                        "data": {
                            "url": "generic-error.png",
                            "width": "40px",
                            "margin": "10px"
                        }
                    },{
                        "type": "paragraph",
                        "data": {"value": "No links available for this fusion event"}
                    }]
            }
        ])
            
    return HttpResponse(json.dumps({"header": header, "items": rows}))

def fusion_event_genes_info(request, fusion_id):
    
    header = []
    rows = []
    
    fusion_event = FusionCell.nodes.get(fusion_cell_id=fusion_id)
    gene_couple = fusion_event.gene_couple[0]
    five_prime_gene = gene_couple.gene_five_prime[0]
    three_prime_gene = gene_couple.gene_three_prime[0]
    five_prime_chromosome = "chr " + five_prime_gene.fromChromosomeToGene[0].chromosome
    three_prime_chromosome = "chr " + three_prime_gene.fromChromosomeToGene[0].chromosome
    
    header = []
    rows = []
    
    gene1 = Gene.nodes.get(symbol=gene_couple.gene1)
    gene2 = Gene.nodes.get(symbol=gene_couple.gene2)
    
    print(gene1)
    print(gene2)
    
    gene_card1 = {"type": "linkable_image", "data": {"title": "See GeneCards page for " + gene1.symbol, "width": "35px", "url": "gc-logo.png", "link": "http://www.genecards.org/cgi-bin/carddisp.pl?gene="+gene1.symbol}}
    gene_cosmic1 = {"type": "linkable_image", "data": {"title": "See Cosmic page for " + gene1.symbol, "width": "35px", "url": "cosmic-icon.png", "link": "http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln="+gene1.symbol}}
    gene_census1 = {"type": "linkable_image", "data": {"title": "See Census page for " + gene1.symbol, "width": "35px", "url": "census-icon.png", "link": "http://cancer.sanger.ac.uk/cosmic/census"}} if gene1.census else {}
    gene_hallmark1 = {"type": "linkable_image", "data": {"title": "See Hallmark page for " + gene1.symbol, "width": "35px", "url": "hallmark-icon.png", "link": "http://cancer.sanger.ac.uk/cosmic/census-page/"+gene1.symbol}} if gene1.hallmark else {}
    
    gene_card2 = {"type": "linkable_image", "data": {"title": "See GeneCards page for " + gene2.symbol, "width": "35px", "url": "gc-logo.png", "link": "http://www.genecards.org/cgi-bin/carddisp.pl?gene="+gene2.symbol}}
    gene_cosmic2 = {"type": "linkable_image", "data": {"title": "See Cosmic page for " + gene2.symbol, "width": "35px", "url": "cosmic-icon.png", "link": "http://cancer.sanger.ac.uk/cosmic/gene/analysis?ln="+gene2.symbol}}
    gene_census2 = {"type": "linkable_image", "data": {"title": "See Census page for " + gene2.symbol, "width": "35px", "url": "census-icon.png", "link": "http://cancer.sanger.ac.uk/cosmic/census/"}} if gene2.census else {}
    gene_hallmark2 = {"type": "linkable_image", "data": {"title": "See Hallmark page for " + gene2.symbol, "width": "35px", "url": "hallmark-icon.png", "link": "http://cancer.sanger.ac.uk/cosmic/census-page/"+gene2.symbol}} if gene2.hallmark else {}
    
    description1 = "<br/><b>Description</b>: " + (gene1.description if gene1.description is not None else "")
    description2 = "<br/><b>Description</b>: " + (gene2.description if gene2.description is not None else "")
    
    rows.append([{
        "type": "multi",
        "layout": "row",
        "layout_align": "start center",
        "subdata":[
            {
                "type": "image",
                "data": {
                    "url": "gene-icon.png",
                    "width": "35px",
                    "margin": "10px"
                }
            },{
                "type": "paragraph",
                "data": {"value": "<b>5' Gene</b>: " + "<a target='_blank' href='http://www.genecards.org/cgi-bin/carddisp.pl?gene="+fusion_event.gene1+"'>" + fusion_event.gene1 + "</a> ("+five_prime_chromosome+")" + description1}
            }, gene_card1, gene_cosmic1, gene_census1, gene_hallmark1]
          }])
    
    rows.append([{
        "type": "multi",
        "layout": "row",
        "layout_align": "start center",
        "subdata":[
            {
                "type": "image",
                "data": {
                    "url": "gene-icon.png",
                    "width": "35px",
                    "margin": "10px"
                }
            },{
                "type": "paragraph",
                "data": {"value": "<b>3' Gene</b>: " + "<a target='_blank' href='http://www.genecards.org/cgi-bin/carddisp.pl?gene="+fusion_event.gene2+"'>" + fusion_event.gene2 + "</a> ("+three_prime_chromosome+")" + description2}
            }, gene_card2, gene_cosmic2, gene_census2, gene_hallmark2]
          }])
    
    return HttpResponse(json.dumps({"header": header, "items": rows}))

def pairs_by_events(request):
    
    header = ["Pair of genes", "No. of fusion events"]
    rows = []
    
    filename = os.path.dirname(__file__) + "/statistics/Pairs_Fusion.txt"
    
    if os.path.isfile(filename):
        for line in open(filename, "r"):
            if len(rows) > 200: continue

            pair = line.rstrip().split("\t")
            rows.append(pair)
    else:
        print("Calculating this kind of statistics for the first time (then saving to file)")
        
        start = datetime.datetime.now()
        
        counter = Counter()
        
        step = 100000
        total = len(FusionCell.nodes) # 990627
        for x in xrange(0, total, step):
            print("Loading all fusion cells between " + str(x) + " and " + str(x+step))
            for fusion in FusionCell.nodes[x:x+step]:
                if fusion.num_algo < 3: continue
                if fusion.annotations is not None and len(filter(lambda x: x in red_list, fusion.annotations)) > 0: continue
                
                if len(counter) % 1000 == 0: print(len(counter))
                
                gene1 = fusion.gene1
                gene2 = fusion.gene2
                
                counter[gene1 + "#" + gene2] += 1
            
        end = datetime.datetime.now()
        
        total_time = end - start
        print("Total time to calculate the statistics (in seconds)", total_time.total_seconds())
        
        print("Building rows...")        
        for pair in counter.most_common():
            rows.append([pair[0], pair[1]])
        
        print("Saving statistics to file...")
        writer = open(filename, "w")
        for pair in counter.most_common():
            writer.write(str(pair[0]) + "\t" + str(pair[1]) + "\n")
        writer.close()

    return HttpResponse(json.dumps({"header": header, "items": rows}))

def fusioncatcher_proteins(request, fusioncatcher_fus_id):
    response = {}
    header = ["Transcript 1", "Transcript 2", "Protein"]
    
    fusions = []
    
    fusions.append(FusionCatcher.nodes.get(fusioncatcher_id = fusioncatcher_fus_id))
    
    met = ()
    rows = []
    for fusion in fusions:
        print(fusion)
        for couple in fusion.with_trans_couple:
            protein = {"type": "text", "label": couple.with_protein.all()[0].protein, "width": "700px"}
            transcript1 = couple.fromTranscriptToCouple[0].transcript
            transcript2 = couple.with_other_transcript[0].transcript
        
            item = [transcript1, transcript2, protein]
            rows.append(item)
    
    response['details'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))


def get_sequence(request, fus_id):
    response = {}
    header = ["Fusion sequence"]
    
    rows = []
    rows.append([FusionCatcher.nodes.get(fusioncatcher_id = fus_id).fusion_sequence])
    
    response['details'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

def get_junction(request, es_id):
    response = {}
    header = ["Junction sequence"]
    
    rows = []
    rows.append([EricScript.nodes.get(ericscript_id = es_id).junction_sequence])

    response['details'] = {"header": header, "items": rows}
    return HttpResponse(json.dumps(response))

INTERSECTION_SYMBOL = u"\u2229"
def get_algorithms_statistics(request, algorithms):
    
    requested_algorithms_raw = algorithms.split("|") if algorithms is not None else []
    print("REQ ALGOS:", requested_algorithms_raw)
    #print(requested_algorithms)
    
    map = {}
    map["FusionCatcher"] = "FC"
    map["EricScript"] = "ES"
    map["Tophat-Fusion"] = "TH"
    map["Jaffa"] = "JA"
    
    inv_map = {v: k for k, v in map.iteritems()}
    
    requested_algorithms = []
    for algo in requested_algorithms_raw:
        requested_algorithms.append(map[algo])

#     combinations = ["FusionCatcher",
#                     "EricScript",
#                     "Tophat-Fusion",
#                     "FusionCatcher"+INTERSECTION_SYMBOL+"EricScript",
#                     "FusionCatcher"+INTERSECTION_SYMBOL+"Tophat-Fusion",
#                     "EricScript"+INTERSECTION_SYMBOL+"Tophat-Fusion",
#                     "FusionCatcher"+INTERSECTION_SYMBOL+"EricScript"+INTERSECTION_SYMBOL+"Tophat-Fusion"]
#     
#     data = ["49032", "929638", "34202", "16475", "4843", "4221", "3294"]
    
    counter = Counter()
    header = []
    rows = []
    r, m = db.cypher_query("MATCH (n:FusionCell) RETURN n.algos,count(*)")
    print(r)
    
    for row in r:
        print("##################")
        set = Set(row[0])
        n = row[1]
        intersection = set.intersection(requested_algorithms)
        
        print(set, intersection, n)
        
        intersection_string_array = []
        for algo in intersection:
            intersection_string_array.append(inv_map[algo])
        
        for subset in reduce(lambda result, x: result + [sorted(subset + [x]) for subset in result], intersection_string_array, [[]]):
            if len(subset) == 0: continue
            counter[INTERSECTION_SYMBOL.join(subset)] += n
        
        print(counter)
        
    for x in counter.most_common():
        header.append(x[0])
        rows.append(x[1])
        
#     for i in range(len(combinations)):
#         combination = combinations[i]
#         print(combination)
#         print(combination, combination.split(INTERSECTION_SYMBOL))
#         if set(combination.split(INTERSECTION_SYMBOL)).issubset(set(requested_algorithms)):
#             header.append(combination)
#             row.append(data[i])
    
    rows = [rows]
    response = {"header": header, "items": rows}
    print("RESPONSE", response)

    return HttpResponse(json.dumps(response))

from distutils import util 
def get_external_database_statistics(request, databases, ccs):

    header = []
    row = []
    
    print(ccs)
    
    if ccs is None: ccs = "False"
    ccs = util.strtobool(ccs)
    
    requested_databases = databases.split("|")
    
    print(databases, requested_databases)
    
    filename = os.path.dirname(__file__) + "/statistics/Description_FusionCell.csv"
    if ccs: filename = os.path.dirname(__file__) + "/statistics/Description_FusionCell_CCS.csv"
    
    print("Reading from filename", filename)
    
    total = 0
    for line in open(filename, "r"):
        total += 1
        if total == 1: continue
        fields = line.rstrip().split(",")
        databases = fields[0].split("/")
        if set(databases).issubset(set(requested_databases)):        
            header.append(INTERSECTION_SYMBOL.join(fields[0].split("/")))
            row.append(fields[1])
    
    rows = [row]
    response = {"header": header, "items": rows}

    return HttpResponse(json.dumps(response))

def get_external_database_statistics_details(request, databases, ccs = False):

    header = ["Annotation", "Subset", "Cardinality"]
    rows = []
    
    requested_databases = databases.split("|")
    
    # print(databases, requested_databases)
    color_messages = {
        "red": "False positive event with high probability",
        "orange": "False positive event with medium probability",
        "green": "True positive event with high probability",
        "grey": "Ambiguous event because tagged with positive and negative evidences"}
    
    if ccs is None: ccs = "False"
    ccs = util.strtobool(ccs)
    
    filename = os.path.dirname(__file__) + "/statistics/Description_FusionCell.csv"
    if ccs: filename = os.path.dirname(__file__) + "/statistics/Description_FusionCell_CCS.csv"
    
    total = 0
    for line in open(filename, "r"):
        total += 1
        if total == 1: continue
        
        fields = line.rstrip().split(",")
        databases = fields[0].split("/")
        if set(databases).issubset(set(requested_databases)):
            element = INTERSECTION_SYMBOL.join(fields[0].split("/"))
            color = fields[2]        
            title = color_messages[color]
            
            icon_object = {
                    "type": "icon",
                    "text_align": "center",
                    "data":{
                        "color": color,
                        "icon": "fa-circle",
                        "title": title
                        }
                }
            rows.append([icon_object, element, int(fields[1])])
    
#     rows = [row]
    rows.sort(key=lambda x: x[2], reverse=True)

    response = {"header": header, "items": rows}

    return HttpResponse(json.dumps(response))

def get_header2():
    return [
        {
            "label": "fusion_id",
            "title": "Fusion ID",
            "tooltip": "Fusion event ID",
            "filters": {
                "title": "Fusion ID filters:",
                "list": [
                    {
                        "type": "select",
                        "key": "fusion_id",
                        "title": "Select a fusion ID:",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]
            }
        },
        {
            "label": "cell_line",
            "title": "Cell line",
            "tooltip": "Cell Line",
            "filters": {
                "title": "Cell line filters:",
                "list": [
                    {
                        "type": "select",
                        "key": "cell_line",
                        "title": "Select a cell line:",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]
            }
        },
        {
            "label": "gene1",
            "title": "5' gene",
            "tooltip": "5' partner gene involved in the fusion",
            "filters": {
                "title": "5' gene filters:",
                "list": [
                    {
                        "type": "select",
                        "key": "gene_five_prime",
                        "title": "Select a five prime gene:",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]
            }
            
        },
        {
            "label": "gene2",
            "title": "3' gene",
            "tooltip": "3' partner gene involved in the fusion",
            "filters": {
                "title": "3' gene filters:",
                "list": [
                    {
                        "type": "select",
                        "key": "gene_three_prime",
                        "title": "Select a three prime gene:",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]
            }
        },
        {
            "label": "cosmic",
            "title": "Cosmic",
            "width": "30px",
            "tooltip": "Known fusion couple in COSMIC Catalogue",
            "filters": {
                "title": "Cosmic filters",
                "list": [
                    {
                        "type": "checkbox",
                        "key": "cosmic",
                        "title": "Include results with Cosmic",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]                
            }
        },
        {
            "label": "mitelman",
            "title": "Mitelman",
            "width": "30px",
            "tooltip": "Known fusion couple in Mitelman Catalogue",
            "filters": {
                "title": "Mitelman filters",
                "list": [
                    {
                        "type": "checkbox",
                        "key": "mitelman",
                        "title": "Include results with Mitelman",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]                
            }
        },
        {
            "label": "fusioncatcher",
            "title": "FusionCatcher",
            "tooltip": "FusionCatcher algorithm",
            "filters": {
                "title": "FusionCatcher filters",
                "list": [
                    {
                        "type": "checkbox",
                        "key": "fusioncatcher",
                        "title": "Include results with FusionCatcher",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]                
            }
        },
        {
            "label": "ericscript",
            "title": "EricScript",
            "tooltip": "EricScript algorithm",
            "filters": {
                "title": "EricScript filters",
                "list": [
                    {
                        "type": "checkbox",
                        "key": "ericscript",
                        "title": "Include results with EricScript",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]                
            }
        },
        {
            "label": "tophat_fusion",
            "title": "Tophat-Fusion",
            "tooltip": "Tophat-Fusion algorithm",
            "filters": {
                "title": "Tophat-Fusion filters",
                "list": [
                    {
                        "type": "checkbox",
                        "key": "tophat_fusion",
                        "title": "Include results with Tophat-Fusion",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]                
            }
        },
        {
            "label": "jaffa",
            "title": "Jaffa",
            "tooltip": "Jaffa algorithm",
            "filters": {
                "title": "Jaffa filters",
                "list": [
                    {
                        "type": "checkbox",
                        "key": "jaffa",
                        "title": "Include results with Jaffa",
                        "placeholder": "",
                        "operators": "LIKE",
                        "chosen_value": ""
                    }
                ]                
            }
        }
    ]

@ensure_csrf_cookie
def search_for_cell_line2(request):
    
    return HttpResponse(json.dumps(search_for_cell_line2_private(request)))

def search_for_cell_line2_private(request, include_FC_only = False):
    data = json.loads(request.body)
    print(data)
    
    global annotations
    
    c_line = data["cell_line"] if "cell_line" in data else "ALL"
    num_algorithms = data["num_algorithms"] if "num_algorithms" in data else 1
    
    # TODO: safe-code
#     if c_line == "ALL" and num_algorithms < 3:
#         return {"structure": {"field_list": []}, "total": 0, "hits": []}
    
    field_list = get_header2()
    
    fusions = []
    
    # Non cosmic + non in database e non verde/giallo/rosso
    novel = data["novel"] if "novel" in data else False
    if novel == "strict":
        novel_condition = "NOT EXISTS(g.cosmic) AND f.quality = 'N/A' AND (NOT EXISTS(f.annotations) OR (" + " AND ".join(["NOT '"+x+"' IN f.annotations" for x in annotations.keys()]) + ")) AND"
        novel_path = "(g:GeneCouple)--"
    elif novel == "relaxed":
        novel_condition = "NOT EXISTS(g.cosmic) AND (f.quality = 'N/A' OR f.quality = 'grey') AND (NOT EXISTS(f.annotations) OR (" + " AND ".join(["NOT '"+x+"' IN f.annotations" for x in annotations.keys()]) + ")) AND"
        novel_path = "(g:GeneCouple)--"
    else:
        novel_condition = ""
        novel_path = ""
    
    if c_line == "ALL":
        # fusions = FusionCell.nodes.filter(num_algo__gte=num_algorithms)
        
        count_query = "MATCH "+novel_path+"(f:FusionCell) WHERE "+novel_condition+" f.num_algo >= "+str(num_algorithms) +" return count(distinct(f))"
        print(count_query)
        r, m = db.cypher_query(count_query)
        total = r[0][0]
        
        query = "MATCH "+novel_path+"(f:FusionCell) WHERE "+novel_condition+" f.num_algo >= "+str(num_algorithms) +" return distinct(f) ORDER BY f.num_algo DESC " + ("SKIP "+ str(data["offset"])+ " LIMIT " + str(data["limit"]) if "offset" in data else "")
        results, meta = db.cypher_query(query)
        fusions = [FusionCell.inflate(row[0]) for row in results]
        #total = len(fusions)
        #fusions = fusions[data["offset"]:data["offset"]+data["limit"]]
    else:
        count_query = "MATCH "+novel_path+"(f:FusionCell)--(c:CellLine) WHERE "+novel_condition+" f.num_algo >= "+str(num_algorithms) +" AND c.cell_line = '"+c_line+"' return count(distinct(f))"
        print(count_query)
        r, m = db.cypher_query(count_query)
        total = r[0][0]

        query = "MATCH "+novel_path+"(f:FusionCell)--(c:CellLine) WHERE "+novel_condition+" f.num_algo >= "+str(num_algorithms) +" AND c.cell_line = '"+c_line+"' return distinct(f) ORDER BY f.num_algo DESC " + ("SKIP "+ str(data["offset"])+ " LIMIT " + str(data["limit"]) if "offset" in data else "")
        print(query)
        results, meta = db.cypher_query(query)
        fusions = [FusionCell.inflate(row[0]) for row in results]

    rows = build_rows2(fusions, c_line, include_FC_only)
    if include_FC_only: field_list = field_list[:-2]
    
    print(len(rows))
    
    response = {"structure": {"field_list": field_list}, "total": total, "hits": rows}
    
    return response

def search_for_disease(request):
    return HttpResponse(json.dumps(search_for_disease_private(request)))

def search_for_disease_private(request):
    
    data = json.loads(request.body)
    print(data)
    
#     if "disease" not in data: return {"structure": {"field_list": get_header2()}, "total": 0, "hits": []}
    
    disease = data["disease"]
    if disease == "ALL":
        return search_for_cell_line2(request)
    
    num_algorithms = data["num_algorithms"] if "num_algorithms" in data else 1
    
    field_list = get_header2()
    
    fusions = []
    
    # Non cosmic + non in database e non verde/giallo/rosso
    novel_path = ""
    novel_condition = ""
    novel = data["novel"] if "novel" in data else False
    if novel == "strict":
        novel_condition = "NOT EXISTS(g.cosmic) AND f.quality = 'N/A' AND (NOT EXISTS(f.annotations) OR (" + " AND ".join(["NOT '"+x+"' IN f.annotations" for x in annotations.keys()]) + ")) AND"
        novel_path = "(g:GeneCouple)--"
    elif novel == "relaxed":
        novel_condition = "NOT EXISTS(g.cosmic) AND (f.quality = 'N/A' OR f.quality = 'grey') AND (NOT EXISTS(f.annotations) OR (" + " AND ".join(["NOT '"+x+"' IN f.annotations" for x in annotations.keys()]) + ")) AND"
        novel_path = "(g:GeneCouple)--"
    
#     if disease == "ALL":
#         fusions = FusionCell.nodes.filter(num_algo__gte=num_algorithms)
#         total = len(fusions)
#         fusions = fusions[data["offset"]:data["offset"]+data["limit"]]
#     else:
    count_query = "MATCH "+novel_path+"(f:FusionCell)--(c:CellLine) WHERE "+novel_condition+" f.num_algo >= "+str(num_algorithms) +" AND c.disease = '"+disease+"' return count(distinct(f))"
    print(count_query)
    r, m = db.cypher_query(count_query)
    total = r[0][0]
    
    query = "MATCH "+novel_path+"(f:FusionCell)--(c:CellLine) WHERE "+novel_condition+" f.num_algo >= "+str(num_algorithms) +" AND c.disease = '"+disease+"' return distinct(f) ORDER BY f.num_algo DESC " + ("SKIP "+ str(data["offset"])+ " LIMIT " + str(data["limit"]) if "offset" in data else "")
    print(query)
    results, meta = db.cypher_query(query)
    fusions = [FusionCell.inflate(row[0]) for row in results]

    print("Fusions", len(fusions))
    rows = build_rows2(fusions)
    print("Rows", len(rows))
    
    return {"structure": {"field_list": field_list}, "total": total, "hits": rows}

from django.utils.encoding import smart_str
from tempfile import NamedTemporaryFile
from wsgiref.util import FileWrapper

def download_disease(request):
    print(request)
    newfile = NamedTemporaryFile(mode="w+t", suffix='.txt')
    print("Creating temporary file at " + str(newfile))
    
    data = search_for_disease_private(request)
    
    print(len(data["hits"]), data["total"])
    for el in data["hits"]:
        newfile.write(str(el) + os.linesep)
    
    newfile.seek(0)
    data = {"filename": "disease_results.txt", "content": newfile.read()}
    
    response = HttpResponse(json.dumps(data))
    
    return response

def download_cell_lines(request):
    print(request)
    newfile = NamedTemporaryFile(mode="w+t", suffix='.txt')
    print("Creating temporary file at " + str(newfile))
    
    data = search_for_cell_line2_private(request)
    
    for el in data["hits"]:
        newfile.write(str(el) + os.linesep)
    
    newfile.seek(0)
    data = {"filename": "cell_line_results.txt", "content": newfile.read()}
    
    response = HttpResponse(json.dumps(data))
    
    return response

def download_genes(request):
    
    newfile = NamedTemporaryFile(suffix='.txt')
    
    data = search_for_gene2_private(request)
    for el in data["hits"]:
        newfile.write(str(el) + "\n")
    
    newfile.seek(0)    
    data = {"filename": "gene_results.txt", "content": newfile.read()}
    
    response = HttpResponse(json.dumps(data))
    
    return response

def download_chromosomes(request):
    
    newfile = NamedTemporaryFile(suffix='.txt')
    
    data = search_for_chromosome2_private(request)
    for el in data["hits"]:
        newfile.write(str(el) + "\n")
    
    newfile.seek(0)    
    data = {"filename": "chromosome_results.txt", "content": newfile.read()}
    
    response = HttpResponse(json.dumps(data))
    
    return response

def download_exons(request):
    
    newfile = NamedTemporaryFile(suffix='.txt')
    
    data = search_for_exon2_private(request)
    for el in data["hits"]:
        newfile.write(str(el) + "\n")
    
    newfile.seek(0)    
    data = {"filename": "exon_results.txt", "content": newfile.read()}
    
    response = HttpResponse(json.dumps(data))
    
    return response

def download_transcripts(request):
    
    newfile = NamedTemporaryFile(suffix='.txt')
    
    data = search_for_transcript2_private(request)
    for el in data["hits"]:
        newfile.write(str(el) + "\n")
    
    newfile.seek(0)    
    data = {"filename": "transcript_results.txt", "content": newfile.read()}
    
    response = HttpResponse(json.dumps(data))
    
    return response

def download_fusions(request):
    
    newfile = NamedTemporaryFile(suffix='.txt')
    
    data = get_fusions2_private(request)
    for el in data["hits"]:
        newfile.write(str(el) + "\n")
    
    newfile.seek(0)    
    data = {"filename": "fusion_results.txt", "content": newfile.read()}
    
    response = HttpResponse(json.dumps(data))
    
    return response

def download_algorithms(request):
    
    newfile = NamedTemporaryFile(suffix='.txt')
    
    data = search_by_algorithm_private(request)
    for el in data["hits"]:
        newfile.write(str(el) + "\n")
    
    newfile.seek(0)    
    data = {"filename": "algorithm_results.txt", "content": newfile.read()}
    
    response = HttpResponse(json.dumps(data))
    
    return response

def download_viruses(request):
    
    newfile = NamedTemporaryFile(suffix='.txt')
    
    data = search_viruses2_private(request)
    for el in data["hits"]:
        newfile.write(str(el) + "\n")
    
    newfile.seek(0)    
    data = {"filename": "virus_results.txt", "content": newfile.read()}
    
    response = HttpResponse(json.dumps(data))
    
    return response