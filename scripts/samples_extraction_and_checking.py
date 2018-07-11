#estraggo tutti i file .gz da tutte le cartelle dove ci sono le mie reads(un file contiene le forward e l'altro file contiene le reverse)
#coding:utf-8
import os
import glob
import sys
import argparse
import numpy as np
import pickle
from datetime import datetime

############################## importo argomenti dal file principale #######################

scripts_folder = os.path.abspath("/".join(sys.argv[0].split("/")[:-1])) 
dictionary_args = pickle.load(open( "%s/dictionary_args" % (scripts_folder),"r"))

if sys.argv[0].split("/")[-1]=="samples_extraction_and_checking.py":
	print "______________________________________________________________"
	print "start of the reads_extraction and checking..."


#################################### funzioni ###################################################################
# funzione per selezionare i files dentro la cartella di interesse tramite glob
def files_selection(directory_type):
	# reads_folder cartella dove sono presenti le reads
	if directory_type=="reads_folder":
		glob_input="%s/*" % (dictionary_args["reads_folder"])
		glob_output=glob.iglob(glob_input)
	# working directory: cartella dove vengono processati i vari output degli scripts della pipeline.
	if directory_type=="working_folder":
		glob_input="./*" 
		glob_output=glob.iglob(glob_input)
	lista_glob_outputs=list(glob_output)
	return lista_glob_outputs


def files_decompression():
	os_input="gunzip -rkf %s" % (dictionary_args["reads_folder"]) # con r decomprimo molti file in maniera ricorsiva(agendo anche dentro altre cartelle),con k mantengo i file compressi senza cancellarli dopo la decompressione,con f costringo il comando a decomprimere senza che esca messaggi di input nel caso ci siano file con lo stesso nome da sovrascrivere
	os.system(os_input)

# funzione per controllare i files dentro le cartelle
def check_folder_content(directory_type,extension_check,extension_type):
	lista_glob_outputs=files_selection(directory_type=directory_type) # seleziona la directory su cui vuoi controllare i files
	#--- controlla presenza dei file in base a una corretta estensione -----------
	if extension_check==True:
		assert extension_type!=None, "please select an extension to be checked" 
		lista_checked_files_extension=[]
		for file_name in lista_glob_outputs:
			if file_name.endswith("%s" % (extension_type)):
				lista_checked_files_extension.append(file_name)
		assert len(lista_checked_files_extension)!=0, "InputError:%s extension file NOT founded in the current %s:impossible to proceed" % (extension_type,directory_type)	
		return lista_checked_files_extension
	#----- altri tipi di controlli ----------------------------------
	#.........
	#-----

##########################################################################################################################################
if __name__=="__main__":
	# se i files sono da decomprimere, perchè si trovano compressi nel gunzip format, allora lo script li decomprime. In caso contrario si prosegue al prossimo script della pipeline
	#----------------- eseguo la decompressione se ci sono file compressi ---------------------------------------
	files_decompression()

	#------ controllo il contenuto della cartella  ---------------------------------------------------------------
	# ricontrolliamo i files dentro la cartella dopo l'estrazione e verifichiamo che sia tutto a posto, per verificare che ci siano i fastq files necessari per proseguire ( es: i file con le reads forward R1 e reverse R2). 
	#Resistuisce la lista dei file fastq sui cui lavoro
	# controllo di files in base a estensione 
	lista_checked_files=check_folder_content(directory_type="reads_folder",extension_check=True,extension_type=".fastq")
	print "checking_valid_input_fastq_extension_files_in_reads_directory:\n",lista_checked_files
	# E' su questi fastq files che posso eseguire l'analisi
	# adesso inizia lo script di selezione delle reads,dove in caso verrà eseguito il merging.
	if sys.argv[0].split("/")[-1]=="samples_extraction_and_checking.py":
		print "______________________________________________________________"
		print "start of the reads_selection script..."
	
