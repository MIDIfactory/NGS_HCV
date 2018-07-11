import os
import sys
from Bio import SeqIO
import argparse
import pickle
from reads_selection_and_merging import reads_to_be_mapped
from blast_and_genotype_reads import get_genotype_of_the_reads_to_be_mapped
from datetime import datetime

############################## importo argomenti dal file principale #######################
scripts_folder = os.path.abspath("/".join(sys.argv[0].split("/")[:-1])) 
dictionary_args = pickle.load(open( "%s/dictionary_args" % (scripts_folder),"r"))
####################################### funzioni ############################################################

databaseHCV_folder=os.path.abspath("/".join((os.path.abspath("/".join(sys.argv[0].split("/")[:-1]))).split("/")[:-1])) + "/databases"

if dictionary_args["genotype"]=="unknown":
	if dictionary_args["user_ref"]==None: # in questa fase ci serve il genotipo solo se usiamo ref del software,se usiamo user_ref il genotipo per il momento non ci interessa.
		diz_genotype=get_genotype_of_the_reads_to_be_mapped(get_genotype=True,print_inf=False,remove_blastout=False,check_genotype=False)
		dictionary_args["genotype"]=diz_genotype['blast_genotype']
		
#print sys.argv[0]
# il mapping può essere eseguito considerando:
#	- single end mode: le reads forward e reverse che vengono da un paired end illumina come reads separate(quindi devi eseguire un mapping separatamente per ciascun set di reads)
#	- paired end mode: le reads forward e reverse sono trattate come un mate pair.
# nel caso del single end mode il comando deve semplicemente includere: il tipo di algoritmo(mem,aln ecc.), la reference e il file di reads di interesse
# nel caso del paired end mode il comando deve includere: l'algoritmo,la reference, il file di reads forward,il file di reads reverse(visto come mate.fq,cioè il "compagno" del forward))
# inoltre:
#- se l'opzione p non è specificata, significa che l'i-esima reads del forward reads file e l'i-esima reads del reverse reads file sono una coppia forward-reverse
#- se l'opzione p non è specificata, significa che l'2iesima read di un file reads e la 2iesima+1 reads dello stesso file reads sono una coppia forward-reverse,quindi significa che non devi mettere il secondo file di reads,perchè sia la forward che la reverse sono nello stesso file.

# funzione per  eseguire the single end bwa mapping(singola reads della coppia o risultato del merging) 
def bwa_mapping_single_end(os_system_action):
	name_reads_file_to_be_mapped=reads_to_be_mapped()
	# gli cambio il nome dell'ouput a seconda che la single end reads che sto mappando sia un elemento del mate pair o risultato del merging del mate pair
	# in caso di merging
	if dictionary_args["reads_processed"]==None and dictionary_args["paired_end_mapping"]==False:
		name_output_file=name_reads_file_to_be_mapped.split("./")[1].split(".")[0] + "_trim_contigs"
	# in caso di singolo elemento della coppia
	elif dictionary_args["reads_processed"]!=None and dictionary_args["paired_end_mapping"]==False:
		name_output_file=name_reads_file_to_be_mapped.split("./")[1].split(".")[0]
	# in caso lo user abbia deciso di non inserire nessuna sua reference si mappa contro le reference fornite all'interno del software
	if dictionary_args["user_ref"]==None:		
		os_input="bwa mem -t %s %s/databaseHCV_bwa/databaseHCV_reference_%s/HCV_%s_%s.TXT %s > %s.sam" % (dictionary_args["processors"],databaseHCV_folder,dictionary_args["genotype"],dictionary_args["genotype"],dictionary_args["gene"],name_reads_file_to_be_mapped,name_output_file)
	# in caso in cui lo user abbia scelto di mappare le reads contro una reference esterna alle reference fornite dal software
	elif dictionary_args["user_ref"]!=None and os_system_action==True: # dobbiamo indexare solo se dobbiamo performare il bwa mapping
		# indexiamo la nuova reference esterna in modo tale che possa essere usata dall'algoritmo bwa mem
		os_input_indexing="bwa index %s" % (dictionary_args["user_ref"]) # input per eseguire l'indexing
		os.system(os_input_indexing)
		os_input="bwa mem -t %s %s %s > %s.sam" % (dictionary_args["processors"],dictionary_args["user_ref"],name_reads_file_to_be_mapped,name_output_file)			
	if os_system_action==True:
		os.system(os_input)
	file_sam_to_be_processed="%s.sam" % (name_output_file)
	return file_sam_to_be_processed

# funzione per eseguire the pair-end end bwa mapping(coppia delle reads non mergiate).
def bwa_mapping_paired_end(os_system_action):
	name_reads_file_to_be_mapped_forward,name_reads_file_to_be_mapped_reverse=reads_to_be_mapped()
	name_output_file=name_reads_file_to_be_mapped_forward.split("./")[1].split("parsed")[0] + "R2_mate"
	# in caso lo user  abbia deciso di non inserire nessuna sua reference si mappa contro le reference fornite all'interno del software
	if dictionary_args["user_ref"]==None:		
		os_input="bwa mem -t %s %s/databaseHCV_bwa/databaseHCV_reference_%s/HCV_%s_%s.TXT %s %s > %s.sam" % (dictionary_args["processors"],databaseHCV_folder,dictionary_args["genotype"],dictionary_args["genotype"],dictionary_args["gene"],name_reads_file_to_be_mapped_forward,name_reads_file_to_be_mapped_reverse,name_output_file)
	# in caso in cui lo user abbia scelto di mappare le reads contro una reference esterna alle reference fornite dal software
	elif dictionary_args["user_ref"]!=None and os_system_action==True: # dobbiamo indexare solo se dobbiamo performare il bwa mapping
		# indexiamo la nuova reference esterna in modo tale che possa essere usata dall'algoritmo bwa mem
		os_input_indexing="bwa index %s" % (dictionary_args["user_ref"]) # input per eseguire l'indexing
		os.system(os_input_indexing)
		os_input="bwa mem -t %s %s %s %s > %s.sam" % (dictionary_args["processors"],dictionary_args["user_ref"],name_reads_file_to_be_mapped_forward,name_reads_file_to_be_mapped_reverse,name_output_file)			
	if os_system_action==True:
		os.system(os_input)
	file_sam_to_be_processed="%s.sam" % (name_output_file)
	return file_sam_to_be_processed

def perform_bwa_mapping():
	#------------------------- selezione single end reads ----------------------------------------------
	# nessun merging,solo un elemento della coppia subisce il mapping 
	if dictionary_args["reads_processed"]!=None and dictionary_args["paired_end_mapping"]==False:
		print "bwa_mapping_with_single_end_of_the_mate_pair_executed..."
		file_sam_to_be_processed=bwa_mapping_single_end(os_system_action=True)

	# --------------- selezione paired end reads con merging fra la forward e reverse ----------------------------------------------------------
	# il risultato del merging subisce il mapping
	elif dictionary_args["reads_processed"]==None and dictionary_args["paired_end_mapping"]==False:
		print "bwa_mapping_with_single_end_of_the_merged_reads_executed..."
		file_sam_to_be_processed=bwa_mapping_single_end(os_system_action=True)

	#----------------- selezione pair-end end reads senza merging --------------------------------------------------
	# sia la forward sia la reverse in contemporanea subiranno il mapping
	elif dictionary_args["reads_processed"]==None and dictionary_args["paired_end_mapping"]==True:
		print "bwa_mapping_with_paired_end_reads_executed..."
		file_sam_to_be_processed=bwa_mapping_paired_end(os_system_action=True)
	
# funzione per ottenere solo il nome del file sam che verrà processatto
def file_sam_to_be_processed():
	if dictionary_args["reads_processed"]!=None and dictionary_args["paired_end_mapping"]==False:
		file_sam_to_be_processed=bwa_mapping_single_end(os_system_action=False)

	# --------------- selezione paired end reads con merging fra la forward e reverse ----------------------------------------------------------
	# il risultato del merging subisce il mapping
	elif dictionary_args["reads_processed"]==None and dictionary_args["paired_end_mapping"]==False:
		file_sam_to_be_processed=bwa_mapping_single_end(os_system_action=False)

	#----------------- selezione pair-end end reads senza merging --------------------------------------------------
	# sia la forward sia la reverse in contemporanea subiranno il mapping
	elif dictionary_args["reads_processed"]==None and dictionary_args["paired_end_mapping"]==True:
		file_sam_to_be_processed=bwa_mapping_paired_end(os_system_action=False)
	return file_sam_to_be_processed    
##################################################################################################################################

if __name__=="__main__":


	# eseguiamo il mapping delle reads contro la reference corretta
	perform_bwa_mapping()
	# verifichiamo il file sam da processare
	file_sam_to_be_processed=file_sam_to_be_processed()
	print "file_sam_to_be_processed: " + file_sam_to_be_processed
	















