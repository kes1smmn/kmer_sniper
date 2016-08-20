## gawkmer - a quick kmer inspector 

<hr>

#### gawkmer takes input as fasta, fastq, or bed file [hg19] and returns the kmer counts for jellyfish database in a csv format.
#### 


<hr>

<hr>

1.	If not already on your system install python 3.5
	* I suggest using anaconda for fewer headaches

2.	Install Jellyfish with the python bindings.

    ``` bash
    wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz
    tar xvf jellyfish-2.2.6.tar.gz
    cd jellyfish-2.2.6
    ./configure --prefix=$HOME --enable-python-binding
    make
    make install
    ```
3.	Install gawkmer

    ``` bash
    git clone
    cd gawkmer
    python setup.py install
    ```

<hr>

Output:
``` python
<sequence_identifier>,<region_uniqueness_rank_score>,<percent_of_region_not_unique>,<kmer_1>,.,.,.,.,<kmer_[regions_size - kmer_size + 1]>
```

region_uniqueness_rank_score = 