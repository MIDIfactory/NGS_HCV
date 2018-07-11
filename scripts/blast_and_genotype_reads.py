from __future__ import division # necessario per assicurare il floating point
from collections import Counter
#coding:utf-8
import sys
import pickle
import argparse
import warnings
import os
from reads_selection_and_merging import reads_to_be_mapped
from datetime import datetime

############################## importo argomenti dal file principale #######################
scripts_folder = os.path.abspath("/".join(sys.argv[0].split("/")[:-1])) # + "/scripts"
dictionary_args = pickle.load(open( "%s/dictionary_args" % (scripts_folder),"r"))

 

############################################### funzioni ########################################################################
databaseHCV_folder=os.path.abspath("/".join((os.path.abspath("/".join(sys.argv[0].split("/")[:-1]))).split("/")[:-1])) + "/databases"



# funzione per eseguire il blast fra le reads da mappare e il database
def perform_blast_agaist_personalized_HCV_databases(os_system_action):
	if dictionary_args["reads_processed"]!=None and dictionary_args["paired_end_mapping"]==False:
   		name_reads_file_to_be_mapped=reads_to_be_mapped()  
		name_file_output=name_reads_file_to_be_mapped.split("./")[1].split(".")[0]   
		blast_imput="blastn -query %s -db %s/databaseHCV_blast/databaseHCV_blast -out blastout_%s.txt -num_alignments 1 -outfmt 6 -num_threads %s" % (name_reads_file_to_be_mapped,databaseHCV_folder,name_file_output,dictionary_args["processors"])
	elif dictionary_args["reads_processed"]==None and dictionary_args["paired_end_mapping"]==False:
  		name_reads_file_to_be_mapped=reads_to_be_mapped()
		name_file_output=name_reads_file_to_be_mapped.split("./")[1].split(".")[0]   
		blast_imput="blastn -query %s -db %s/databaseHCV_blast/databaseHCV_blast -out blastout_%s.txt -num_alignments 1 -outfmt 6 -num_threads %s" % (name_reads_file_to_be_mapped,databaseHCV_folder,name_file_output,dictionary_args["processors"])
	elif dictionary_args["reads_processed"]==None and dictionary_args["paired_end_mapping"]==True:
		# ne usiamo cmq solo una: se una appartiene a un genotipo anche all'altra appartiene allo stesso genotipo
  		name_reads_file_to_be_mapped_forward,name_reads_file_to_be_mapped_reverse=reads_to_be_mapped()
		name_file_output=name_reads_file_to_be_mapped_forward.split("./")[1].split(".")[0]  
		blast_imput="blastn -query %s -db %s/databaseHCV_blast/databaseHCV_blast -out blastout_%s.txt -num_alignments 1 -outfmt 6 -num_threads %s" % (name_reads_file_to_be_mapped_forward,databaseHCV_folder,name_file_output,dictionary_args["processors"])			
	if os_system_action==True:
		os.system(blast_imput)
	elif os_system_action==False:
		return name_file_output
# funzione per estrarre il genotipo delle reads da mappare dal file di output uscito dal blast.
def get_genotype_of_the_reads_to_be_mapped(get_genotype,print_inf,remove_blastout,check_genotype):
	name_file_output=perform_blast_agaist_personalized_HCV_databases(os_system_action=False)
	file_imput="./blastout_%s.txt" % (name_file_output)
	with open(file_imput,"r") as blastoutput:
		#cerco i genotipi che desidero nella colonna degli id del blastoutput e calcolo la prevalenza
		count_a=0
		count_b=0
		t=0
		lista_genotipi=[] # metto solo i tipi di genotipi presenti(quindi un rappresentante per ogni tipo)
		all_genotype=[] # metto tutti i genotipi che trova il for loop.
		for record in blastoutput:
			t+=1
			genotype = str(record).strip().split()[1][0:2] # genotipo trovato dal loop per una determinata record
			if genotype not in lista_genotipi:        
				lista_genotipi.append(genotype)
			all_genotype.append(genotype)
		diz_prevalence={}  
		for i in lista_genotipi:
			count = all_genotype.count(i) # metodo di liste per contare il numero di elementi "i" della lista.   
			prevalence = count/t * 100
			diz_prevalence[i]=prevalence # riporto prevalence in un dizionario. Le chiavi sono il tipo di genotipo,i cui valori rappresentano la prevalenza di quel genotipo nel file di records
		diz_genotype={} # dizionario dove inserisco il genotipo risultante dal blast che poi devo mettere ,se necessario, come imput nel file del processamento.
		for a in diz_prevalence:
			if diz_prevalence[a] == max(diz_prevalence.values()):
				genotype_sample=a # rappresenta la key del dizionario(e quindi il tipo di genotipo) con il massimo valore di prevalenza fra tutti i valori di prevalenza del dizionario, ovvero "diz_prevalence.values()"
				if "1" not in genotype_sample:
					diz_genotype["blast_genotype"]=str(genotype_sample)[0:1]
				else:
					diz_genotype["blast_genotype"]=genotype_sample
	#----------- printiamo varie informazioni sull'analisi del blastoutput ------------------------
	if print_inf==True:
		print "prevalance_genotype:",round(max(diz_prevalence.values()),3)
		if max(diz_prevalence.values())<90:
			print "NGS_HCV warning: the prevalence of the genotype extracted is below 90%"
			print "the genotype of the sample is: " + genotype_sample
		elif max(diz_prevalence.values())<=50:
			warnings.warn("impossible to determine a unique genotype for this sample,the genotype with higher prevalence will be chosen by default",RuntimeWarning)
			print "genotype chosen is: " + genotype_sample
		else:
			print "the genotype of the sample is: " + diz_genotype["blast_genotype"]
		#----------- controlliamo che il genotipo trovato sia fra i genotipi delle ref a disposizione -----------	
	if check_genotype==True:
		lista_valid_genotype=["1a","1b","2","3","4"]
		if diz_genotype["blast_genotype"] not in lista_valid_genotype:
			os_input="rm ./blastout_%s.txt" % (name_file_output)
			os.system(os_input) # eliminiamo il blastout dato che diventerebbe un file inutile se interrompiamo la sessione
			raise Exception("GenotypeExtractionError: the genotype of your input reads does not match with the genotype of the available references")
			sys.exit("execution interrupted")
		#---------- rimuoviamo il file di output del blast -----------------------------------------
	if remove_blastout==True:
		#rimuoviamo il blastoutput quando non ci serve piu
		os_input="rm ./blastout_%s.txt" % (name_file_output)
		os.system(os_input)	
	return diz_genotype






####################################################################################################################################  
if __name__=="__main__":
	if sys.argv[0].split("/")[-1]=="blast_and_genotype_reads.py":
		print "______________________________________________________________"
		print "start of the blast execution..."
	# ------------------  eseguiamo il blast contro il database di sequenze di HCV ---------------------------------------------
	perform_blast_agaist_personalized_HCV_databases(os_system_action=True)
	#-------------------- estraiamo il genotipo ------------------------------------------------------------------------------
	if sys.argv[0].split("/")[-1]=="blast_and_genotype_reads.py":
		print "start of the genotype extraction....."
	diz_genotype=get_genotype_of_the_reads_to_be_mapped(get_genotype=True,print_inf=True,remove_blastout=False,check_genotype=True)
	# ----------------------- scriviamo qui che ora possiamo iniziare il bwa mapping -----------------------------
	if sys.argv[0].split("/")[-1]=="blast_and_genotype_reads.py":
		print "______________________________________________________________"
		print "start of the bwa_mapping..."

