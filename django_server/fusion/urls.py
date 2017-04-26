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
    url(r'^chromosomes/(.+)/(.+)/', views.search_for_chromosome, name='search_for_chromosome'),
#     url(r'^chromosomes/(.+)/(.+)/(.+)/(.+)/', views.search_for_chromosome, name='search_for_chromosome'),
    url(r'^cell_line/(.+)/', views.search_for_cell_line, name='search_for_cell_line'),
    url(r'^gene/(.+)/(.+)/(.+)/', views.search_for_gene, name='search_for_gene'),
#     url(r'^genes/(.+)/(.+)/(.+)/', views.search_for_pair_gene, name='search_for_pair_gene'),
#     url(r'^genes/(.+)/(.+)/', views.search_for_single_gene, name='search_for_single_gene'),
#     url(r'^genes/(.+)/', views.search_for_single_gene, name='search_for_single_gene'),
    url(r'^exon/single/(.+)/(.+)/', views.search_for_single_exon, name='search_for_single_exon'),
    url(r'^exon/pair/(.+)/(.+)/(.+)/', views.search_for_pair_exon, name='search_for_pair_exon'),
    url(r'^transcript/single/(.+)/(.+)/', views.search_for_single_transcript, name='search_for_single_transcript'),
    url(r'^transcript/pair/(.+)/(.+)/(.+)/', views.search_for_pair_transcript, name='search_for_pair_transcript'),
    url(r'^fusion_information/(.+)/(.+)/(.+)/(.+)/(.+)/', views.search_for_fusion_information, name='search_for_fusion_information'),
    url(r'^statistics_all/', views.statistics_all, name='statistics_all'),
    url(r'^statistics_by_chromosome/(.+)/', views.statistics_by_chromosome, name='statistics_by_chromosome'),
    url(r'^fusion_by_chromosome/', views.fusion_by_chromosome, name='statistics_all'),
    url(r'^download_data/', views.download_data, name='download_data'),
    url(r'^genstats/', views.generate_statistics, name='generate_statistics'),
    url(r'^get_distribution/(\w+)/(\w+)/(\d*)/?(\w*)/?$', views.get_distribution, name='get_distribution'),
    url(r'^get_single_distribution/(.+)/(.+)/', views.get_single_distribution, name='get_single_distribution'),
    url(r'^search_for_disease/(.+)/', views.search_for_disease, name='search_for_disease'),
    url(r'^events/(.*)/', views.fusion_events, name='fusionEvents'),
    url(r'^count_genes/(.+)/', views.count_genes, name='count_genes'),
    url(r'^print_file/', views.print_file, name='print_file'),
    url(r'^cell_lines/', views.cell_lines, name='cell_lines'),
    url(r'^get_chromosomes_cell_lines/(.+)', views.get_chromosomes_cell_lines, name='get_chromosomes_cell_lines'),
    url(r'^search_indels_by_region/(.+)/(.+)/(.+)', views.search_indels_by_region, name='search_indels_by_region'),
    url(r'^show_info/(.+)', views.show_info, name='show_info'),
    url(r'^fusioncatcher/(.+)/', views.build_fc_table, name='build_fc_table'),
    url(r'^ericscript/(.+)/', views.build_es_table, name='build_es_table'),
    url(r'^tophat/(.+)/', views.build_th_table, name='build_th_table'),
]
