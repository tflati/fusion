"""FusionCatcherDjango URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.10/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
from fusion import views

urlpatterns = [
    url(r'^cell_lines/', views.cell_lines, name='cell_lines'),
    url(r'^cell_lines_simple/?(.*)', views.cell_lines_simple, name='cell_lines_simple'),
    url(r'^diseases_simple/?(.*)', views.diseases_simple, name='diseases_simple'),
    url(r'^chromosomes2/', views.search_for_chromosome2, name='search_for_chromosome'),
    url(r'^disease/', views.search_for_disease, name='search_for_disease'),
    url(r'^chromosomes/', views.chromosomes, name='chromosomes'),
    url(r'^chromosomes_simple/', views.chromosomes_simple, name='chromosomes_simple'),
    url(r'^cell_line_events/(.+)/', views.search_for_cell_line_with_chromosome, name='cellline_chromosomes'),
    url(r'^cell_line_events_simple/(.+)/', views.search_for_cell_line_events, name='cellline_chromosomes_simple'),
    url(r'^cell_line2/', views.search_for_cell_line2, name='search_for_cell_line'),
    url(r'^gene2/', views.search_for_gene2, name='search_for_gene'),
    url(r'^genes/', views.get_genes, name='get_genes'),
    url(r'^genes_three_prime/(.+)?', views.get_genes_three_prime, name='get_genes_three_prime'),
    url(r'^genes_five_prime/(.+)?', views.get_genes_five_prime, name='get_genes_five_prime'),
    url(r'^transcripts_three_prime/(.+)?', views.get_transcripts_three_prime, name='get_transcripts_three_prime'),
    url(r'^transcripts_five_prime/(.+)?', views.get_transcripts_five_prime, name='get_transcripts_five_prime'),
    url(r'^exons_three_prime/(.+)?', views.get_exons_three_prime, name='get_exons_three_prime'),
    url(r'^exons_five_prime/(.+)?', views.get_exons_five_prime, name='get_exons_five_prime'),
    url(r'^exon2/', views.search_for_exon2, name='search_for_exon'),
    url(r'^exons/', views.exons, name='list_exons'),
    url(r'^transcripts/', views.transcripts, name='list_transcripts'),
    url(r'^transcript2/', views.search_for_transcript2, name='search_for_transcript'),
    url(r'^virus2/', views.search_viruses2, name='search_viruses'),
    url(r'^virus/', views.viruses, name='viruses'),
    url(r'^statistics_all/', views.statistics_all, name='statistics_all'),
    url(r'^statistics_by_chromosome/(.+)/', views.statistics_by_chromosome, name='statistics_by_chromosome'),
    url(r'^fusion_by_chromosome/', views.fusion_by_chromosome, name='statistics_all'),
    url(r'^download_data2/', views.download_data2, name='download_data2'),
    url(r'^get_distribution/(\w+)/(\w+)/(\d*)/?(\w*)/?$', views.get_distribution, name='get_distribution'),
    url(r'^get_single_distribution/(.+)/(.+)/', views.get_single_distribution, name='get_single_distribution'),
    url(r'^events/(.*)/', views.fusion_events, name='fusionEvents'),
    url(r'^count_genes/(.+)/', views.count_genes, name='count_genes'),
    url(r'^print_file/', views.print_file, name='print_file'),
    url(r'^cell_lines_of_fusion/(.*)/', views.cell_lines_of_fusion, name='cell_lines_of_fusion'),
    url(r'^search_indels_by_region/(.+)/(.+)/(.+)', views.search_indels_by_region, name='search_indels_by_region'),
    url(r'^show_info/(.+)', views.show_info, name='show_info'),
    url(r'^fusioncatcher/(.+)/(.+)/', views.build_fc_table, name='build_fc_table'),
    url(r'^ericscript/(.+)/(.+)/', views.build_es_table, name='build_es_table'),
    url(r'^tophat/(.+)/(.+)/', views.build_th_table, name='build_th_table'),
    url(r'^mapping_algorithms/', views.get_mapping_algorithms, name='mapping_algorithms'),
    url(r'^fusion_descriptions/', views.get_fusion_descriptions, name='fusion_descriptions'),
    url(r'^predicted_effects/', views.get_predicted_effects, name='predicted_effects'),
    url(r'^fusions2/', views.get_fusions2, name='fusions'),
    url(r'^algorithms/', views.get_algorithms, name='algorithm'),
    url(r'^search_by_algorithm/', views.search_by_algorithm, name='search_by_algorithm'),
    url(r'^oncofuse/(.+)/(.+)/(.+)/', views.get_oncofuse, name='get_oncofuse'),
    url(r'^oncofuse/(.+)/(.+)/', views.get_oncofuse, name='get_oncofuse'),
    url(r'^oncofuse/(.+)/', views.get_oncofuse, name='get_oncofuse'),
    url(r'^algorithms_statistics/(.+)?/', views.get_algorithms_statistics, name='get_algorithms_statistics'),
    url(r'^fusioncatcher_proteins/(.+)/', views.fusioncatcher_proteins, name='fusioncatcher_proteins'),
    url(r'^get_distribution_virus_by_disease/', views.get_distribution_virus_by_disease, name='get_distribution_virus_by_disease'),
    url(r'^pairs_by_events/', views.pairs_by_events, name='pairs_by_events'),
    url(r'^get_fusion_by_disease/', views.get_fusion_by_disease, name='get_fusion_by_disease'),
    url(r'^sequence/(.+)/', views.get_sequence, name='get_sequence'),
    url(r'^junction/(.+)/', views.get_junction, name='get_junction'),
    url(r'^download_disease/', views.download_disease, name='download_disease'),
    url(r'^download_cell_lines/', views.download_cell_lines, name='download_cell_lines'),
    url(r'^download_genes/', views.download_genes, name='download_genes'),
    url(r'^download_chromosomes/', views.download_chromosomes, name='download_chromosomes'),
    url(r'^download_exons/', views.download_exons, name='download_exons'),
    url(r'^download_transcripts/', views.download_transcripts, name='download_transcripts'),
    url(r'^download_fusions/', views.download_fusions, name='download_fusions'),
    url(r'^download_algorithms/', views.download_algorithms, name='download_algorithms'),
    url(r'^download_viruses/', views.download_viruses, name='download_viruses'),
    url(r'^generate_statistics/', views.generate_statistics, name='generate_statistics'),
    url(r'^get_external_database_statistics/(.+)/(.+)/', views.get_external_database_statistics, name='get_external_database_statistics'),
    url(r'^get_external_database_statistics_details/(.+)/(.+)/', views.get_external_database_statistics_details, name='get_external_database_statistics_details'),
    url(r'^get_descriptions_statistics/', views.get_descriptions_statistics, name='get_descriptions_statistics'),
    url(r'^get_event_map/', views.get_event_map, name='get_event_map'),
    url(r'^cell_line_info/(.+)/', views.cell_line_info, name='cell_line_info'),
    url(r'^cell_line_minor_info/(.+)/', views.cell_line_minor_info, name='cell_line_minor_info'),
    url(r'^gene_couples/(.+)/', views.gene_couples, name='gene_couples'),
    url(r'^cell_line_chromosome_info/(.+)/', views.cell_line_chromosome_info, name='chromosome_info'),
    url(r'^fusion_event_info/(.+)/', views.fusion_event_info, name='fusion_event_info'),
    url(r'^fusion_event_minor_info/(.+)/', views.fusion_event_minor_info, name='fusion_event_minor_info'),
    url(r'^fusion_event_genes_info/(.+)/', views.fusion_event_genes_info, name='fusion_event_genes_info'),
    url(r'^fusion_event_link/(.+)/', views.fusion_event_link, name='fusion_event_link'),
]
