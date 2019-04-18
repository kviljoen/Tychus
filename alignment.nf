#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017 Chris Dean

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:

 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.

 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

//General configuation variables moved to nextflow.config

// Display help menu
if(params.help) {
	log.info ''
	log.info 'Tychus - Alignment Pipeline'
	log.info ''
	log.info 'Usage: '
	log.info '    nextflow alignment.nf -profile alignment [options]'
	log.info ''
	log.info 'General Options: '
	log.info '    --read_pairs      DIR		Directory of paired FASTQ files'
	log.info '    --genome          FILE		Path to the FASTA formatted reference database'
	log.info '    --amr_db          FILE		Path to the (single line) FASTA formatted resistance database'
	log.info '    --vf_db           FILE		Path to the (single line) FASTA formatted virulence database'
	log.info '    --plasmid_db      FILE		Path to the (single line) FASTA formatted plasmid database'
	log.info '    --draft           FILE		Path to the FASTA formatted draft databases'
	log.info '    --threads         INT		Number of threads to use for each process'
	log.info '    --alignment_out_dir         DIR		Directory to write output files to'
	log.info ''
	log.info 'Trimmomatic Options: '
	log.info '    --leading         INT		Remove leading low quality or N bases'
	log.info '    --trailing        INT		Remove trailing low quality or N bases'
	log.info '    --slidingwindow   INT		Scan read with a sliding window'
	log.info '    --minlen          INT		Drop reads below INT bases long'
	log.info '    --adapters        FILE		FASTA formatted adapter sequences'
	log.info '    --saveTrimmed 	BOOL		Save trimmed FastQ file intermediates'
	log.info ''
	log.info 'kSNP Options: '
	log.info '    --ML              BOOL		Estimate maximum likelihood tree'
	log.info '    --NJ              BOOL		Estimate neighbor joining tree'
	log.info '    --min_frac        DECIMAL		Minimum fraction of genomes with locus'
	log.info ''
	log.info 'Figtree Options: '
	log.info '    --JPEG            BOOL		Convert newick tree to annotated JPEG'
	log.info '    --PDF             BOOL		Convert newick tree to annotated PDF'
	log.info '    --PNG             BOOL		Convert newick tree to annotated PNG'
	log.info '    --SVG             BOOL		Convert newick tree to annotated SVG'
	log.info ''
	return
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

// Returns a tuple of read pairs in the form
// [dataset_id, forward.fq, reverse.fq] where
// the dataset_id is the shared prefix from
// the two paired FASTQ files.
Channel
        .fromFilePairs(params.read_pairs, flat: true)
	.ifEmpty { exit 1, "Read pairs could not be found: ${params.read_pairs}" }
        .set { trimmomatic_read_pairs }

// Define multiple-use value channel to guide users to intermediate files in not saved to publishdir
Channel.fromPath("$baseDir/assets/where_are_my_files.txt")
       .into{ch_where_trimmomatic}
       

// Validate user-inputs
if( params.genome ) {
        genome = file(params.genome)
        if( !genome.exists() ) exit 1, "Genome file could not be found: ${params.genome}"
}

if( params.amr_db ) {
        amr_db = file(params.amr_db)
        if( !amr_db.exists() ) exit 1, "AMR file could not be found: ${params.amr_db}"
}

if( params.vf_db ) {
        vf_db = file(params.vf_db)
        if( !vf_db.exists() ) exit 1, "Virulence file could not be found: ${params.vf_db}"
}

if( params.plasmid_db ) {
        plasmid_db = file(params.plasmid_db)
        if( !plasmid_db.exists() ) exit 1, "Plasmid file could not be found: ${params.plasmid_db}"
}

//if( params.annot_db) {
//	annot_db = file(params.annot_db)
//	if( !annot_db.exists() ) exit 1, "Annotation file could not be found: ${params.annot_db}"
//}

if( params.draft ) {
        draft_path = params.draft.substring(0, params.draft.lastIndexOf("/"))
	draft_genomes = Channel.fromPath(params.draft).toSortedList()
        if( !draft_genomes.size() == 0 ) exit 1, "Draft genome file(s) could not be found: ${params.draft}"
}

/*
 * Build reference genome index with Bowtie2
 */
process BuildGenomeIndex {
	tag { "${genome.baseName}" }
	cache 'deep'
	
	input:
	file genome

	output:
	file 'genome.index*' into genome_index

	"""
	bowtie2-build $genome genome.index --threads ${task.cpus}
	"""
}

/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 */
process RunQC {
        publishDir "${params.alignment_out_dir}/Preprocessing", mode: 'copy',
            saveAs: {filename ->
                if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
                else if (params.saveTrimmed && filename != "where_are_my_files.txt") "trimmed_fastq/$filename"
                else null
            }
	tag { dataset_id }
	cache 'deep'
	
        input:
        set dataset_id, file(forward), file(reverse) from trimmomatic_read_pairs
	file wherearemyfiles from ch_where_trimmomatic.collect()

        output:
        set dataset_id, file("${dataset_id}_1P.fastq"), file("${dataset_id}_2P.fastq") into (amr_read_pairs, plasmid_read_pairs, vf_read_pairs, genome_read_pairs)
	file "where_are_my_files.txt"

        """
        ${JAVA}/java "-Xmx3G" -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar PE -threads ${task.cpus} $forward $reverse -baseout ${dataset_id} ILLUMINACLIP:Trimmomatic-0.36/adapters/${params.adapters}:2:30:10:3:TRUE LEADING:${params.leading} TRAILING:${params.trailing} SLIDINGWINDOW:${params.slidingwindow} MINLEN:${params.minlen}
        mv ${dataset_id}_1P ${dataset_id}_1P.fastq
        mv ${dataset_id}_2P ${dataset_id}_2P.fastq
        """
}

if( params.amr_db ) {
	/*
	 * Build resistance database index with Bowtie2
	 */
	process BuildAMRIndex {
		tag { "${amr_db.baseName}" }
		cache 'deep'
		
		input:
        	file amr_db

        	output:
        	file 'amr.index*' into amr_index

        	"""
        	bowtie2-build $amr_db amr.index --threads ${task.cpus}
		"""
	}

	/*
         * Align reads to resistance database with Bowtie2
         */
	process AMRAlignment {
        	publishDir "${params.alignment_out_dir}/Alignment", mode: "move", pattern: "*.bam"
		tag { dataset_id }
		cache 'deep'
		
        	input:
        	set dataset_id, file(forward), file(reverse) from amr_read_pairs
        	file index from amr_index

        	output:
        	set dataset_id, file("${dataset_id}_amr_alignment.sam") into amr_sam_files
        	set dataset_id, file("${dataset_id}_amr_alignment.bam") into amr_bam_files

        	"""
        	bowtie2 -p ${task.cpus} -x amr.index -1 $forward -2 $reverse -S ${dataset_id}_amr_alignment.sam
        	samtools view -bS ${dataset_id}_amr_alignment.sam | samtools sort -@ ${task.cpus} -o ${dataset_id}_amr_alignment.bam
        	"""
	}

	process AMRResistome {
        	publishDir "${params.alignment_out_dir}/Resistome", mode: "move"
        	tag { dataset_id }
		cache 'deep'
		
        	input:
        	file amr_db
        	set dataset_id, file(amr_sam) from amr_sam_files

        	output:
        	set dataset_id, file("${dataset_id}_amr_gene_resistome.tsv") into amr_gene_level

        	"""
		csa -ref_fp ${amr_db} -sam_fp ${amr_sam} -min 5 -max 100 -skip 5 -t 0 -samples 1 -out_fp "${dataset_id}_amr_gene_resistome.tsv"
        	"""
	}
}

if( params.vf_db ) {
	/*
         * Build virulence database index with Bowtie2
         */
	process BuildVFIndex {
		tag { "${vf_db.baseName}" }
		cache 'deep' 
		
		input:
        	file vf_db

        	output:
        	file 'vf.index*' into vf_index

        	"""
        	bowtie2-build $vf_db vf.index --threads ${task.cpus}
		"""
	}
	/*
         * Align reads to virulence factor database with Bowtie2
         */
	process VFAlignment {
        	publishDir "${params.alignment_out_dir}/Alignment", mode: "move", pattern: "*.bam"
		tag { dataset_id }
		cache 'deep'
		
        	input:
        	set dataset_id, file(forward), file(reverse) from vf_read_pairs
        	file index from vf_index

        	output:
        	set dataset_id, file("${dataset_id}_vf_alignment.sam") into vf_sam_files
        	set dataset_id, file("${dataset_id}_vf_alignment.bam") into vf_bam_files

        	"""
        	bowtie2 -p ${task.cpus} -x vf.index -1 $forward -2 $reverse -S ${dataset_id}_vf_alignment.sam
        	samtools view -bS ${dataset_id}_vf_alignment.sam | samtools sort -@ ${task.cpus} -o ${dataset_id}_vf_alignment.bam
        	"""
	}

	process VFResistome {
        	publishDir "${params.alignment_out_dir}/Resistome", mode: "move"
        	tag { dataset_id }
		cache 'deep'
		
        	input:
        	file vf_db
        	set dataset_id, file(vf_sam) from vf_sam_files

        	output:
        	set dataset_id, file("${dataset_id}_vf_gene_resistome.tsv") into vf_gene_level

        	"""
        	csa -ref_fp ${vf_db} -sam_fp ${vf_sam} -min 5 -max 100 -skip 5 -t 0 -samples 1 -out_fp "${dataset_id}_vf_gene_resistome.tsv"
        	"""
	}
}

if( params.plasmid_db ) { //KL: downloaded prebuilt from plsdb
	/*
         * Build plasmid index with Bowtie2
         */
	process BuildPlasmidIndex {
		tag { "${plasmid_db.baseName}" }
		cache 'deep'
		
		input:
       		file plasmid_db

        	output:
        	file 'plasmid.index*' into plasmid_index

        	"""
        	bowtie2-build $plasmid_db plasmid.index --threads ${task.cpus}
		"""
	}
	/*
         * Align reads to plasmid database with Bowtie2
         */
	process PlasmidAlignment {
        	publishDir "${params.alignment_out_dir}/Alignment", mode: "move", pattern: "*.bam"
        	tag { dataset_id }
		cache 'deep'
		
        	input:
        	set dataset_id, file(forward), file(reverse) from plasmid_read_pairs
        	file index from plasmid_index

        	output:
        	set dataset_id, file("${dataset_id}_plasmid_alignment.sam") into plasmid_sam_files
        	set dataset_id, file("${dataset_id}_plasmid_alignment.bam") into plasmid_bam_files

        	"""
        	bowtie2 -p ${task.cpus} -x plasmid.index -1 $forward -2 $reverse -S ${dataset_id}_plasmid_alignment.sam
        	samtools view -bS ${dataset_id}_plasmid_alignment.sam | samtools sort -@ ${task.cpus} -o ${dataset_id}_plasmid_alignment.bam
        	"""
	}

	process PlasmidResistome {
        	publishDir "${params.alignment_out_dir}/Resistome", mode: "move"
        	tag { dataset_id }
		cache 'deep'
		
        	input:
        	file plasmid_db
        	set dataset_id, file(plasmid_sam) from plasmid_sam_files

        	output:
        	set dataset_id, file("${dataset_id}_plasmid_gene_resistome.tsv") into plasmid_gene_level

        	"""
        	csa -ref_fp ${plasmid_db} -sam_fp ${plasmid_sam} -min 5 -max 100 -skip 5 -t 0 -samples 1 -out_fp "${dataset_id}_plasmid_gene_resistome.tsv"
        	"""
	}
}

/*
 * Align reads to reference genome with Bowtie2
 */
process GenomeAlignment {
	publishDir "${params.alignment_out_dir}/Alignment", mode: "copy", pattern: "*.bam"
	tag { dataset_id }
	cache 'deep'
	
	input:
	set dataset_id, file(forward), file(reverse) from genome_read_pairs
	file index from genome_index

	output:
	set dataset_id, file("${dataset_id}_genome_alignment.bam") into genome_bam_files
        set dataset_id, file("${dataset_id}_genome_alignment.bai") into genome_index_files

	"""
	bowtie2 -p ${task.cpus} -x genome.index -1 $forward -2 $reverse -S ${dataset_id}_genome_alignment.sam
	samtools view -bS ${dataset_id}_genome_alignment.sam | samtools sort -@ ${task.cpus} -o ${dataset_id}_genome_alignment.bam
        samtools index ${dataset_id}_genome_alignment.bam ${dataset_id}_genome_alignment.bai
	"""
}

/*
 * Call SNPs with Freebayes and integrate them into reference genome with BCFtools
 */
process BuildConsensusSequence {
	tag { dataset_id }
	publishDir "${params.alignment_out_dir}/Consensus", mode: "copy"
	cache 'deep'
	
	input:
	set dataset_id, file(bam) from genome_bam_files
	set dataset_id, file(bai) from genome_index_files
	file genome

	output:
	file("${dataset_id}_consensus.fa") into consensus_files
	file("${dataset_id}_in_list.txt") into ksnp3_configuration

	"""
	freebayes -p 1 -f ${genome} $bam | bgzip -c > ${dataset_id}_genome_variants.vcf.gz
	tabix ${dataset_id}_genome_variants.vcf.gz
	cat $genome | bcftools consensus ${dataset_id}_genome_variants.vcf.gz > ${dataset_id}_consensus.fa
	echo -e "$params.alignment_out_dir/Consensus/${dataset_id}_consensus.fa\t$dataset_id" >> ${dataset_id}_in_list.txt
	"""
}

if( params.draft ) {
	/*
	 * Create configuration file for kSNP3 using the draft assemblies and user-input reference genome
	 */
	process kSNPDraftAndGenomeConfiguration {
		echo true

		input:
                file draft from draft_genomes

                output:
                file("genome_paths.txt") into genome_config

                shell:
                '''
                #!/bin/sh
                echo "!{genome}\t!{genome.baseName}" > genome_paths.txt
                for d in !{draft};
                do
                        echo "!{draft_path}/${d}\t${d%.*}" >> genome_paths.txt
                done
                '''
	}
}

else {
	/*
	 * Create configuration file for kSNP3 using the user-input reference genome
	 */
	process kSNPGenomeConfiguration {
		echo true

		//storeDir 'temporary_files'
		
		//ref_genome = file(params.genome)
		
		input:
		file genome

		output:
		file("genome_paths.txt") into genome_config
		file("$genome") into kchooser_genome

		shell:
		'''
		#!/bin/sh
		base=`echo !{genome} | cut -f1 -d '.'`
		fp=`readlink !{genome}`
		echo "${fp}\t${base}" > genome_paths.txt
		'''
	}
}

/*
 * Build phylogenies with kSNP3
 */
process BuildPhylogenies {
	publishDir "${params.alignment_out_dir}/SNPsAndPhylogenies", mode: "copy"
	tag { "ConfigurationFiles" }

	input:
	file kchooser_config from genome_config
	file ksnp3_config from ksnp3_configuration.toList()

	output:
	file("Trees/*.tre") into phylogenetic_trees
	file("SNPs/*") into polymorphisms

	shell:
	'''
	#!/bin/sh
	MakeFasta !{kchooser_config} MF.fasta > /dev/null
	Kchooser MF.fasta > /dev/null
	optimum_k=`grep "The optimum" Kchooser.report | tr -dc '0-9'`
	cat !{kchooser_config} > in_list
	cat !{ksnp3_config} >> in_list
	if [ !{params.ML} && !{params.NJ} ]
	then
		kSNP3 -in in_list -outdir kSNP3_results -k "${optimum_k}" -ML -NJ -core -min_frac !{params.min_frac} >> /dev/null
	elif [ !{params.NJ} ]
	then
		kSNP3 -in in_list -outdir kSNP3_results -k "${optimum_k}" -NJ -core -min_frac !{params.min_frac} >> /dev/null
	else
		kSNP3 -in in_list -outdir kSNP3_results -k "${optimum_k}" -ML -core -min_frac !{params.min_frac} >> /dev/null
	fi
	mkdir Trees
	mkdir SNPs
	mv kSNP3_results/*.tre Trees
	mv kSNP3_results/* SNPs
	rm -rf !{params.work_dir}
	'''
}

// Ignore files with nothing in them. This will
// occur when there aren't enough samples for kSNP
// to build the appropriate trees.
phylogenetic_trees.flatten()
	.filter { file -> !file.isEmpty() }
	.set { trees }

/*
 * Build phylogenetic trees with Figtree
 */
process ConvertNewickToPDF {
	publishDir "${params.alignment_out_dir}/SNPsAndPhylogenies/TreeImages", mode: "move"
	cache 'deep'
	
	input:
	file tree from trees

	output:
	file "*"

	script:
	base = tree.baseName

	shell:
	"""
	#!/bin/sh
        if [ ${params.PDF} ]
        then
               	java -jar ${FIGTREE}/figtree.jar -graphic PDF ${tree} ${base}.pdf
        elif [ ${params.SVG} ]
        then
               	java -jar ${FIGTREE}/figtree.jar -graphic SVG ${tree} ${base}.svg
        elif [ ${params.JPEG} ]
        then
               	java -jar ${FIGTREE}/figtree.jar -graphic JPEG ${tree} ${base}.jpg
        else
               	java -jar ${FIGTREE}/figtree.jar -graphic PNG ${tree} ${base}.png
       	fi
	chmod o+w ${base}*
	"""
}


/*
 *
 * Completion e-mail notification
 *
 */
workflow.onComplete {
  
    def subject = "[uct-yamp] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[uct-yamp] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
   
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    //if(workflow.container) email_fields['summary']['Docker image'] = workflow.container

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[uct-yamp] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[uct-yamp] Sent summary e-mail to $params.email (mail)"
        }
    }
}
