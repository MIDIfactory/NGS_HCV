#coding: utf-8
### in questo file riporto semplicemente tutte le funzioni utili già pronte che importo in altri scripts
def batch_iterator(iterator, batch_size):
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True  # Make sure we loop once
    while entry:
        batch = []
        while len(batch) < batch_size:
            try:
                entry = iterator.next()
            except StopIteration:
                entry = None
            if entry is None:
                # End of file
                break
            batch.append(entry)
        if batch:
            yield batch

'''
una volta usata questa funzione su un "iteratore" posso usare:
list(itertools.islice(iteratore, 1)) # per selezionare solo il primo batch di elementi (es: i primi mille elementi) e metterli in una lista
list(itertools.islice(iteratore, 2)) # per selezionare oltre al primo anche il secondo batch di elementi(es: i primi elementi più i successivi mille elementi) 
e metterli in una lista dove il primo elemento della lista(elemento 0) è a sua volta una lista che rappresenta il primo batch(primi mille elementi) 
e il secondo elemento della lista (elemento 1) è a sua volta una lista che rappresenta il secondo batch(gli altri mille elementi)
# e così via...
'''
'''
# Creo la funzione per ottenere l'indexed database for the blast(it is necessary,only one time,because blastn accepts only database indexed)
def get_blast_indexed_databases():
	os_imput="makeblastdb -in %s/databaseHCV_blast/outdatabaseHCV.fasta -title databaseHCV -dbtype nucl -out databaseHCV_blast -parse_seqids;mv *databaseHCV_blast* %s/databaseHCV_blast" % (databaseHCV_blast_folder,databaseHCV_blast_folder)
	os.system(os_imput)
'''


