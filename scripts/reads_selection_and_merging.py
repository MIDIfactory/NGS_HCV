#coding: utf-8
import os
from Bio import SeqIO
import sys
import argparse
from samples_extraction_and_checking import check_folder_content 
from utility import batch_iterator # funzione utile che importo dallo script utility.
import itertools
import pickle
from datetime import datetime
import numpy as np
import random
############################## importo argomenti dal file principale #######################
scripts_folder = os.path.abspath("/".join(sys.argv[0].split("/")[:-1])) 
dictionary_args = pickle.load(open( "%s/dictionary_args" % (scripts_folder),"r"))


############################################### funzioni ########################################################################

# seleziono le reads che devo processare o su cui fare il merging (se voglio fare il merging). Le reads risultanti verranno usate per il mapping
# faccio il parsing come Seq_records  di queste reads usando SeqIO di bioPython.Devo creare un file fasta e  quindi perdo informazioni sulla qualità. Il parse trasforma la qualità in un formato non accettato da mothur.
# Nota: le sequenze derivano da Illumina. Illumina associa a ogni base una misura di qualità. Nell'ultima versione si usa il phred score, lo stesso del Sanger,ovvero Q=-10logP
# dove P è la prob. che la base sia incorretta. va da 0 a 40,40 significa bassa prob. che la base sia sbagliata.
def single_reads_selection_parsing_export(name_file_to_be_parsed,single_in_pair_index):
	
	#-------------------------------- reads_selection ---------------------------------------------
	lista_checked_files=check_folder_content(directory_type="reads_folder",extension_type=".fastq",extension_check=True)
	# facciamo in modo che questo messaggio venga printato solo quando si esegue lo script reads_selection_and_merging e non quando questa funzione viene utilizzata in altri script
	if sys.argv[0].split("/")[-1]=="reads_selection_and_merging.py":
		print "input_reads_file_processed_selected_by_the_user:","%s/%s" % (dictionary_args["reads_folder"],name_file_to_be_parsed) 
	try:
		index_reads_file=lista_checked_files.index("%s/%s" % (dictionary_args["reads_folder"],name_file_to_be_parsed.split(".")[0] + ".fastq"))
	except ValueError as InputError:	
		raise RuntimeError("InputError:The input reads_file_to_be_processed is not included in your reads_folder,please insert a valid reads_file input")
		 
	reads_file=lista_checked_files[index_reads_file]
	#------------------------- reads_parsing ---------------------------------------------
	# parsiamo i file (che ci dà un iterator) e questo iterator lo transformiamo in lista
	parse_reads_file=SeqIO.parse(reads_file, "fastq")
	if dictionary_args["batching"]!=None:
		# analizzo solo un batch di reads
		# uso itertool per prendere un certo batch(i primi 100,i secondi 100 ecc..). Altra soluzione è usare SeqIO.index(SeqIo.parse...)		
		batching_reads_file=batch_iterator(parse_reads_file,int(dictionary_args["batching"])) # crea un iterator del file di reads dove ogni elemento dell'iteratore è una lista di un certo numero di reads(un batch di reads)
		reads_parsed=list(itertools.islice(batching_reads_file,None)) # posso prendere solo le prime x sequenze,le prime x + le seconde x,oppure le prime x + seconde + le terze x ecc... In questo caso ho settato stop=None così prendo tutti i batch
		random_number= np.random.randint(low=0,high=len(reads_parsed),size=None)# genero un numero causale compreso fra l'indice minimo e massimo di reads_parsed che contiene un certo tot di batch del file di reads.
		index=random_number
		# nel caso in cui stiamo facendo un merging o un pair end mapping abbiamo bisogno che il batch di reads estratto sia uguale in entrambi i files
		# quindi nel secondo file non dobbiamo estratte un indice della lista di batch casuale,ma uno predeterminato,che corrisponde con quello del primo file.
		if single_in_pair_index!=None:
			index=single_in_pair_index
		reads_parsed=reads_parsed[index] # a questo punto posso selezionare quale gruppo di sequenze considerare,se le prime x(elemento 0 della lista reads_parsed),oppure le seconde x(elemento 1 della lista reads_parsed e così via). In questo caso prendo un batch a caso	
	elif dictionary_args["batching"]==None:
		reads_parsed=list(parse_reads_file)
	# ------------------ export reads_file ----------------------------------------
	with open('%s_parsed.fasta' % (name_file_to_be_parsed.split(".")[0]), 'w+') as file_reads: # name_file_to_be_parsed.split(".")[0] in caso ci siano punti nel nome del file non li considero
		for record in reads_parsed:
			SeqIO.write(record,file_reads,"fasta")
	return index
			
# funzione per seleziona una coppia di reads o per fare il mapping con il risultato del loro merging oppure per fare il mapping considerandole come una coppia
def paired_end_reads_selection(forward_reads_file,reverse_reads_file):
	# selezioniamo la coppia di reads su cui vogliamo fare il merging e la parsiamo
	index=single_reads_selection_parsing_export(forward_reads_file,None)
	single_reads_selection_parsing_export(reverse_reads_file,index)
	reads_input_forward="./%s_parsed.fasta" % (forward_reads_file.split(".")[0]) # forward_reads_file.split(".")[0]  in caso ci siano punti nel nome del file non li considero
	reads_input_reverse="./%s_parsed.fasta" % (reverse_reads_file.split(".")[0]) # stesso motivo
	return  reads_input_forward,reads_input_reverse
	 
	

#-------------------- funzione per selezionare i file di reads che verrà usato per fare il mapping sulla reference ---------------------
# se l'utente ha scelto di non fare il merging,significa che ha scelto un file di reads con cui vuole fare il mapping o cmq le operazioni successive 
# quindi in questo caso reads_processed!=None e il file di reads che da adesso in poi verrà usato la versione post-parsing(batchata oppure no) di questo file di reads scelto dall'utente
# se invece l'utente ha scelto di fare il merging quindi ha specificato rf e rrv e quindi dictionary_args["reads_processed"]==None,cioè r non è specificato,
# allora il file reads che adesso in po verrà processato sarà quello risultante dal merging.
def reads_to_be_mapped():
	#--------------------- selezione nel caso di assenza di merging: una solo file di reads come input ------------------------
	if dictionary_args["reads_processed"]!=None and dictionary_args["paired_end_mapping"]==False:
		name_reads_file_to_be_mapped="./%s_parsed.fasta" % ((dictionary_args["reads_processed"]).split(".")[0])
		return name_reads_file_to_be_mapped
	#------------------ selezione in caso di presenza di merging fra forward e reverse: due files di reads di input ---------
	elif dictionary_args["reads_processed"]==None and dictionary_args["paired_end_mapping"]==False:
		# Nota: mothur come output produce sempre un file con estesione fasta a prescindere che tu decida di mergiare due fasta e due fastq
		name_reads_file_to_be_mapped="./%s_parsed.trim.contigs.fasta" % ((dictionary_args["reads_forward"]).split(".")[0])
		return name_reads_file_to_be_mapped
	#------------------ selezione in caso di assenza di merging,come considerando le reads come una coppia da mappare: due files di reads di input ---------
	elif dictionary_args["reads_processed"]==None and  dictionary_args["paired_end_mapping"]==True:
		name_reads_file_to_be_mapped_forward,name_reads_file_to_be_mapped_reverse=paired_end_reads_selection(forward_reads_file=dictionary_args["reads_forward"],reverse_reads_file=dictionary_args["reads_reversed"])
		return name_reads_file_to_be_mapped_forward,name_reads_file_to_be_mapped_reverse
		
############################################################################################################################################
if __name__=="__main__":
	# se non ti interessa fare il merging,allora seleziona il singolo file di reads all'interno della cartella delle reads. 
	#------------------------- selezione single end reads ----------------------------------------------
	# nessun merging,solo un elemento della coppia subisce il mapping 
	if dictionary_args["reads_processed"]!=None and dictionary_args["paired_end_mapping"]==False:
		print "parsing_input_files..."
		index=single_reads_selection_parsing_export(name_file_to_be_parsed=dictionary_args["reads_processed"],single_in_pair_index=None)
		lista_checked_files=check_folder_content(directory_type="working_folder",extension_check=True,extension_type="_parsed.fasta")
		if len(lista_checked_files)==0:
			raise RuntimeError("no_parsed.fasta file foundend in the working directory,parsing failed")
		print "checking_valid_parsed_extension_files_in_the_working_directory:\n",lista_checked_files
		name_reads_file_to_be_mapped=reads_to_be_mapped()
		print "name_reads_file_that_will_be_mapped: ", name_reads_file_to_be_mapped
	# --------------- selezione paired end reads con merging fra la forward e reverse ----------------------------------------------------------
	# il risultato del merging subisce il mapping
	elif dictionary_args["reads_processed"]==None and dictionary_args["paired_end_mapping"]==False:

		name_reads_file_to_be_mapped_forward,name_reads_file_to_be_mapped_reverse=paired_end_reads_selection(forward_reads_file=dictionary_args["reads_forward"],reverse_reads_file=dictionary_args["reads_reversed"])
		os_input_mothur='mothur "#make.contigs(ffasta=%s,rfasta=%s, processors=%s)"' % (name_reads_file_to_be_mapped_forward,name_reads_file_to_be_mapped_reverse,dictionary_args["processors"])
		os.system(os_input_mothur)
		# eliminiamo file inutili
		remove_files="rm *.logfile; rm *.num.temp; rm *scrap*; rm *_parsed.fasta;rm *.qual"
		os.system(remove_files)
		# Nota: il formato di output del merging è sempre un fasta. a prescinde che dall'estensione dei due files che stai mergiando (2 fasta o 2 fastq)
		# Nota: se decidi di mergiare due fasta e non due fastq, il merging ovviamente non terrà in considerazione il quality score(phred quality score) delle basi delle sequenze 	
		lista_checked_files=check_folder_content(directory_type="working_folder",extension_check=True,extension_type=".trim.contigs.fasta")
		print "checking_valid_parsed_trim_contigs_fasta_extension_files_in_the_working_directory:\n",lista_checked_files
		if len(lista_checked_files)==0:
			raise RuntimeError(".trim.contigs.fasta file foundend in the working directory,parsing and merging failed")
		name_reads_file_to_be_mapped=reads_to_be_mapped()
		print "name_reads_file_that_will_be_mapped: ", name_reads_file_to_be_mapped
	#----------------- selezione pairend end reads senza merging --------------------------------------------------
	# sia la forward sia la reverse in contemporanea subiranno il mapping
	elif dictionary_args["reads_processed"]==None and dictionary_args["paired_end_mapping"]==True:
		print "parsing_input_files..."
		name_reads_file_to_be_mapped_forward,name_reads_file_to_be_mapped_reverse=reads_to_be_mapped()
		lista_checked_files=check_folder_content(directory_type="working_folder",extension_check=True,extension_type="_parsed.fasta")
		if len(lista_checked_files)==0:
			raise RuntimeError("no_parsed.fasta file foundend in the working directory,parsing failed")
		print "checking_valid_parsed_extension_files_in_the_working_directory:\n",lista_checked_files
		print "name_reads_pair_that_will_be_mapped:[%s,%s]" % (name_reads_file_to_be_mapped_forward,name_reads_file_to_be_mapped_reverse)

	#------------------- indichiamo qui che sta per iniziare il bwa mapping ------------------------------------
	if sys.argv[0].split("/")[-1]=="reads_selection_and_merging.py": # per evitare che se lo importo come modulo mi printi questa cosa anche da altre parti
		if dictionary_args["genotype"]!="unknown":
			print "______________________________________________________________"
			print "start of the bwa mapping..."
		elif dictionary_args["genotype"]=="unknown":
			print "the genotype of the reads is unknown: a blast will be performed before bwa mapping"

