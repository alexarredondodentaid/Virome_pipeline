#!/bin/bash

# --- Configuració del LOG ---
LOG_FILE="pipeline_$(date +%Y%m%d_%H%M%S).log"

# Funció per escriure al log i a la pantalla simultàniament
log_msg() {
    echo -e "$1" | tee -a "$LOG_FILE"
}

# --- Configuració de rutes i recursos ---
MAIN_DIR=$(pwd)
DB_KNEAD="/home/alexarredondo/Documents/Software/KneadData/HumanDB" 
TRIMMOMATIC="/home/alexarredondo/anaconda3/envs/kneaddata/bin/Trimmomatic-0.33"
SPADES="/home/alexarredondo/Documents/Software/SPAdes-4.2.0-Linux/bin/spades.py"
DVF_SCRIPT="/home/alexarredondo/Documents/Software/DeepVirFinder/dvf.py"
CHECKV_DB="/home/alexarredondo/Documents/Software/CheckV/checkv-db-v1.5"
GENOMAD_DB="/home/alexarredondo/Documents/2025_Viroma_test/genomad_db/genomad_db/"

THREADS=100
MEMORY=120

# Inicialització del fitxer log
echo "--- INICI DE LA PIPELINE: $(date) ---" > "$LOG_FILE"

# Creació de directoris finals (els que NO s'esborraran)
mkdir -p checkv_res genomad_res Mapped_res

# Carreguem conda
source ~/anaconda3/etc/profile.d/conda.sh

log_msg " "
log_msg "################################################################################################"
log_msg "#                                                                                              #"
log_msg "# Benvingudi a la pipeline de processament de dades metagenòmiques per a la detecció de virus! #"
log_msg "# Realitzada per Indispensable, PhD                                                            #"
log_msg "#                                                                                              #"
log_msg "# Comencem amb el trimatge de les lectures amb Cutadapt...                                     #"
log_msg "#                                                                                              #"
log_msg "################################################################################################"
log_msg " "

# --- BUCLE PRINCIPAL PER MOSTRA ---
for R1 in *_1.fastq.gz; do # revisar l'extensió dels arxius a utilitzar
    # Identifiquem la parella R2
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
    
    if [[ ! -f "$R2" ]]; then
        log_msg "ALERTA: No es troba la parella de $R1. Saltant mostra..."
        continue
    fi

    # Extraurem el nom de la mostra (ex: ESHE042)
    SAMPLE=$(basename "$R1" _1.fastq.gz)
    
    log_msg " "
    log_msg "========================================================================"
    log_msg " INICIANT PROCESSAMENT DE LA MOSTRA: $SAMPLE"
    log_msg "========================================================================"
    log_msg " "


    # 1. CUTADAPT (Trimatge)
    mkdir -p "temp_trim_${SAMPLE}"
    conda activate cutadapt_env
    log_msg " "
    log_msg "###################################"
    log_msg "# [$SAMPLE] Executant Cutadapt... #"
    log_msg "###################################"
    log_msg " "
    # -u 10 i -U 10 eliminen els primers 10 nucleòtids de R1 i R2, respectivament. -q 30,30 fa un trimatge de baixa qualitat a les lectures.
    # -o serveix per indicar el nom de sortida del fitxer trimat de R1, i -p per al fitxer trimat de R2. El paràmetre --cores=0 utilitza tots els nuclis disponibles.
    cutadapt -u 10 -U 10 -q 30,30 \
        -o "temp_trim_${SAMPLE}/${SAMPLE}_1.trim.fastq" \
        -p "temp_trim_${SAMPLE}/${SAMPLE}_2.trim.fastq" \
        "$R1" "$R2" --cores=0 --report=minimal >> "$LOG_FILE" 2>&1
    conda deactivate

    # 2. KNEADDATA (Eliminació d'ADN humà)
    mkdir -p "temp_knead_${SAMPLE}"
    conda activate kneaddata
    log_msg " "
    log_msg "############################################################"
    log_msg "# [$SAMPLE] Traiem les seqüències humanes amb KneadData... #"
    log_msg "############################################################"
    log_msg " "
    # -i1 i -i2 indiquen els fitxers trimats de R1 i R2, respectivament.
    # -db és la base de dades de referència per a eliminar seqüències humanes. 
    # -o és el directori de sortida. 
    # --trimmomatic especifica la ruta a l'executable de Trimmomatic que KneadData utilitzarà internament.
    kneaddata -i1 "temp_trim_${SAMPLE}/${SAMPLE}_1.trim.fastq" \
              -i2 "temp_trim_${SAMPLE}/${SAMPLE}_2.trim.fastq" \
              -db "$DB_KNEAD" \
              -o "temp_knead_${SAMPLE}" \
              --threads "$THREADS" \
              --trimmomatic "$TRIMMOMATIC" >> "$LOG_FILE" 2>&1
    conda deactivate

    # Busquem els fitxers resultants de kneaddata (paired)
    K1=$(find "temp_knead_${SAMPLE}" -name "*_paired_1.fastq")
    K2=$(find "temp_knead_${SAMPLE}" -name "*_paired_2.fastq")

    # 3. METASPADES (Assemblatge)
    log_msg " "
    log_msg "##########################################################"
    log_msg "# [$SAMPLE] Comencem amb l'assemblatge amb metaSPAdes... #"
    log_msg "##########################################################"
    log_msg " "
    # --meta indica que estem treballant amb dades metagenòmiques. -1 i -2 indiquen els fitxers de lectures parells resultants de kneaddata. 
    # -t i -m especifiquen el nombre de fils i la quantitat de memòria a utilitzar, respectivament. -o és el directori de sortida on es guardaran els resultats de l'assemblatge.
    python3 "$SPADES" --meta -1 "$K1" -2 "$K2" \
        -t "$THREADS" -m "$MEMORY" -o "temp_spades_${SAMPLE}" >> "$LOG_FILE" 2>&1

    CONTIGS="temp_spades_${SAMPLE}/contigs.fasta"
    if [[ ! -f "$CONTIGS" ]]; then
        log_msg "ERROR: No s'ha generat contigs.fasta per a $SAMPLE. Saltant a la següent."
        continue
    fi

    # 4. DETECCCIÓ VIRAL (VirSorter2 i DeepVirFinder)
    conda activate vs2
    log_msg " "
    log_msg "#####################################"
    log_msg "# [$SAMPLE] Executant VirSorter2... #"
    log_msg "#####################################"
    log_msg " "
    # -w és el directori de treball on VirSorter2 guardarà els resultats. 
    # -i és el fitxer d'entrada amb els contigs a analitzar. 
    # --min-length estableix la longitud mínima dels contigs a considerar. 
    # --include-groups especifica els grups virals a incloure en l'anàlisi. 
    # --keep-original-seq manté les seqüències originals dels contigs en els resultats. 
    # -j indica el nombre de fils a utilitzar.
    virsorter run -w "temp_vs2_${SAMPLE}" -i "$CONTIGS" --min-length 1000 \
        --include-groups dsDNAphage,ssDNA --keep-original-seq -j "$THREADS" >> "$LOG_FILE" 2>&1
    
    # VirSorter2 genera un fitxer final-viral-score.tsv on es poden filtrar els contigs virals segons el seu score i altres criteris. 
    # En aquest cas, es seleccionen els contigs que tenen un score de 0.9 o superior, almenys 1000 bases de longitud, i que compleixen altres criteris de qualitat (columnes 7 i 8).
    awk -F'\t' 'NR==1 || ($2>=0.9 && $5>=1000 && $7>=1 && $8>=5)' "temp_vs2_${SAMPLE}/final-viral-score.tsv" | cut -f1 > "temp_vs2_${SAMPLE}/vs2.ids"
    conda deactivate

    conda activate dvf
    log_msg " "
    log_msg "########################################"
    log_msg "# [$SAMPLE] Executant DeepVirFinder... #"
    log_msg "########################################"
    log_msg " "
    # -l estableix la longitud mínima dels contigs a analitzar.
    python "$DVF_SCRIPT" -i "$CONTIGS" -l 1000 -c "$THREADS" -o "temp_dvf_${SAMPLE}" >> "$LOG_FILE" 2>&1
    DVF_PRED=$(find "temp_dvf_${SAMPLE}" -name "*_dvfpred.txt")
    # Filtrem per aquells contigs que tenen una probabilitat de ser virals (columna 2) de 0.9 o superior, i un valor de p (columna 3) de 0.01 o inferior.
    awk -F'\t' 'NR==1 || ($2>=0.9 && $3<=0.01)' "$DVF_PRED" | cut -f1 > "temp_dvf_${SAMPLE}/dvf.ids"
    conda deactivate

    # 5. UNIÓ I FILTRATGE
    log_msg " "
    log_msg "##########################################"
    log_msg "# [$SAMPLE] Unim contigs de vs2 i dvf... #"
    log_msg "##########################################"
    log_msg " "
    mkdir -p "temp_union_${SAMPLE}"
    cat "temp_vs2_${SAMPLE}/vs2.ids" "temp_dvf_${SAMPLE}/dvf.ids" | sort | uniq > "temp_union_${SAMPLE}/viral_union.ids"
    # -r és per a fer una cerca recursiva, -n per a mostrar només els noms dels contigs que coincideixen amb els IDs de la llista.
    # -f indica el fitxer amb els IDs a buscar, i el darrer argument és el fitxer de contigs on es farà la cerca.
    seqkit grep -r -n -f "temp_union_${SAMPLE}/viral_union.ids" "$CONTIGS" > "temp_union_${SAMPLE}/viral_final.fasta"

    # 6. CHECKV (Qualitat)
    conda activate checkv
    log_msg " "
    log_msg "################################################"
    log_msg "# [$SAMPLE] Avaluant la qualitat amb CheckV... #"
    log_msg "################################################"
    log_msg " "
    export CHECKVDB="$CHECKV_DB"
    # end_to_end és el mode de CheckV que realitza una anàlisi completa de les seqüències d'entrada, incloent la predicció de gens, l'estimació de la completitud i la classificació de la qualitat dels contigs virals.
    # -t indica el nombre de fils a utilitzar, i -d especifica la base de dades de CheckV a utilitzar per a l'anàlisi.
    checkv end_to_end "temp_union_${SAMPLE}/viral_final.fasta" "checkv_res/$SAMPLE" -t "$THREADS" -d "$CHECKV_DB" >> "$LOG_FILE" 2>&1
    # Triem la columna 8 per la nostra versio de checkV, que és on es classifiquen els contigs segons la seva qualitat (Complete, High-quality, Medium-quality, Low-quality, etc.). Cal revisar-ho segons la versió de CheckV que s'estigui utilitzant, ja que les columnes poden variar.
    awk '$8=="Complete" || $8=="High-quality" || $8=="Medium-quality"{print $1}' \
        "checkv_res/$SAMPLE/quality_summary.tsv" > "checkv_res/$SAMPLE/good.ids"
    
    seqkit grep -r -n -f "checkv_res/$SAMPLE/good.ids" "temp_union_${SAMPLE}/viral_final.fasta" > "temp_union_${SAMPLE}/viral_good.fasta"
    conda deactivate

    # 7. GENOMAD (Taxonomia)
    conda activate genomad
    log_msg " "
    log_msg "################################################"
    log_msg "# [$SAMPLE] Assignant taxonomia amb Genomad... #"
    log_msg "################################################"
    log_msg " "
    genomad end-to-end --cleanup --min-score 0.7 \
        "temp_union_${SAMPLE}/viral_good.fasta" "genomad_res/$SAMPLE/" "$GENOMAD_DB" >> "$LOG_FILE" 2>&1
    conda deactivate

    # 8. MAPEIG I ABUNDÀNCIES
    mkdir -p "Mapped_res/$SAMPLE"
    conda activate bowtie2
    log_msg " "
    log_msg "##################################################"
    log_msg "# [$SAMPLE] Calculant abundàncies amb Bowtie2... #"
    log_msg "##################################################"
    log_msg " "
    # Li donem a Bowtie2 el fitxer de contigs virals de bona qualitat com a referència per construir l'índex, i després mapejem les lectures parells de kneaddata contra aquest índex.
    bowtie2-build "temp_union_${SAMPLE}/viral_good.fasta" "temp_union_${SAMPLE}/idx" >> "$LOG_FILE" 2>&1
    # -x indica l'índex de referència, -1 i -2 són els fitxers de lectures parells a mapejar.
    # -p el nombre de fils a utilitzar. El resultat del mapeig es converteix a format BAM (molt menys pesat que SAM) i s'ordena
    bowtie2 -x "temp_union_${SAMPLE}/idx" -1 "$K1" -2 "$K2" -p "$THREADS" 2>> "$LOG_FILE" | \
        samtools view -bS - | samtools sort -o "Mapped_res/$SAMPLE/${SAMPLE}_sorted.bam"
    samtools index "Mapped_res/$SAMPLE/${SAMPLE}_sorted.bam"
    # Indexa el fitxer BAM resultant per a una consulta més ràpida, i finalment s'obtenen les estadístiques d'abundància dels contigs virals mapejats.
    samtools idxstats "Mapped_res/$SAMPLE/${SAMPLE}_sorted.bam" > "Mapped_res/$SAMPLE/${SAMPLE}_counts.txt"
    conda deactivate

    conda activate prodigal
    log_msg " "
    log_msg "###############################################################"
    log_msg "# [$SAMPLE] Predicció de gens amb Prodigal i featureCounts... #"
    log_msg "###############################################################"
    log_msg " "
    # -f gff indica que volem la sortida en format GFF, i -p meta és el mode de Prodigal optimitzat per a dades metagenòmiques, que no assumeix que les seqüències d'entrada siguin genomes complets.
    prodigal -i "temp_union_${SAMPLE}/viral_good.fasta" -o "Mapped_res/$SAMPLE/${SAMPLE}_viral_genes.gff" -f gff -p meta >> "$LOG_FILE" 2>&1
    # -p indica que les lectures són parells, -a és el fitxer d'entrada amb les anotacions de gens en format GFF, -t és el tipus d'element a comptar (CDS per a gens codificants), 
    # -g és el camp del GFF que s'utilitzarà com a identificador de gen. 
    featureCounts -p -a "Mapped_res/$SAMPLE/${SAMPLE}_viral_genes.gff" -t CDS -g ID \
        -o "Mapped_res/$SAMPLE/gene_counts_matrix.txt" "Mapped_res/$SAMPLE/${SAMPLE}_sorted.bam" >> "$LOG_FILE" 2>&1
    conda deactivate

    # 9. NETEJA D'ARXIUS INTERMEDIS ---
    log_msg " "
    log_msg "#######################################################"
    log_msg "# [$SAMPLE] Netejant arxius temporals de la mostra... #"
    log_msg "#######################################################"
    log_msg " "
    rm -rf "temp_trim_${SAMPLE}" "temp_knead_${SAMPLE}" "temp_spades_${SAMPLE}" "temp_vs2_${SAMPLE}" "temp_dvf_${SAMPLE}" "temp_union_${SAMPLE}"
    
    log_msg "FINALITZAT: $SAMPLE"
    log_msg "========================================================================"
done

log_msg " "
log_msg "#################################################################"
log_msg "# Tot fet! Revisa el fitxer $LOG_FILE per a detalls del procés. #"
log_msg "#################################################################"
