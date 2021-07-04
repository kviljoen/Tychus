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

// General configuration variables moved to nextflow.config


// Display help menu
if(params.help) {
	log.info ''
	log.info 'Tychus - Assembly Pipeline'
	log.info ''
	log.info 'Usage: '
	log.info '    nextflow assembly.nf -profile assembly [options]'
	log.info ''
	log.info 'General Options: '
	log.info '    --read_pairs      DIR		Directory of paired FASTQ files'
	log.info '    --threads         INT		Number of threads to use for each process'
	log.info '    --assembly_out_dir          DIR		Directory to write output files to'
	log.info ''
	log.info 'Trimmomatic Options: '
	log.info '    --leading         INT		Remove leading low quality or N bases'
	log.info '    --trailing        INT		Remove trailing low quality or N bases'
	log.info '    --slidingwindow   INT		Scan read with a sliding window'
	log.info '    --minlen          INT		Drop reads below INT bases long'
	log.info '    --adapters        STR		FASTA formatted adapter file'
	log.info ''
	log.info 'Prokka Options: '
	log.info '    --genus           STR		Target genus'
	log.info '    --species         STR		Target species'
	log.info ''
	return
}

params.name = false
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


// Header log info
log.info "==================================="
log.info " uct-tychus "
log.info "==================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.read_pairs
summary['OS']		= System.getProperty("os.name")
summary['OS.arch']	= System.getProperty("os.arch") 
summary['OS.version']	= System.getProperty("os.version")
summary['javaversion'] = System.getProperty("java.version") //Java Runtime Environment version
summary['javaVMname'] = System.getProperty("java.vm.name") //Java Virtual Machine implementation name
summary['javaVMVersion'] = System.getProperty("java.vm.version") //Java Virtual Machine implementation version
//Gets starting time		
sysdate = new java.util.Date() 
summary['User']		= System.getProperty("user.name") //User's account name
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.assembly_out_dir
summary['Working dir']    = workflow.workDir
summary['Container']      = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.genus) {
	summary['Genus'] = params.genus
}
if(params.species) {
	summary['Species'] = params.species
}
if(params.email) {
    summary['E-mail Address'] = params.email
}
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="



/*
 * Remove adapter sequences and low quality base pairs with Trimmomatic
 */
process RunQC {
	publishDir "${params.assembly_out_dir}/PreProcessing", mode: "copy"
	tag { dataset_id }
	cache 'deep'
	
        input:
        set dataset_id, file(forward), file(reverse) from trimmomatic_read_pairs

        output:
        set dataset_id, file("${dataset_id}_1P.fastq"), file("${dataset_id}_2P.fastq") into (abyss_read_pairs, velvet_read_pairs, spades_read_pairs, idba_read_pairs, kmer_genie_read_pairs)

        """
        java -jar ${TRIMMOMATIC}/trimmomatic-0.36.jar PE -threads ${task.cpus} $forward $reverse -baseout ${dataset_id} ILLUMINACLIP:Trimmomatic-0.36/adapters/${params.adapters}:2:30:10:3:TRUE LEADING:${params.leading} TRAILING:${params.trailing} SLIDINGWINDOW:${params.slidingwindow} MINLEN:${params.minlen}
        mv ${dataset_id}_1P ${dataset_id}_1P.fastq
        mv ${dataset_id}_2P ${dataset_id}_2P.fastq
        """
}

/*
 * Choose a best kmer to build the underlying De Bruijn graph for Abyss and Velvet with KmerGenie
 */
process IdentifyBestKmer {
	tag { dataset_id }
	cache 'deep'
	
	input:
	set dataset_id, file(forward), file(reverse) from kmer_genie_read_pairs

	output:
	file("${dataset_id}_best-k.txt") into (best_abyss_kmer_results, best_velvet_kmer_results)
	set dataset_id, file("${dataset_id}_forward_kg.fq"), file("${dataset_id}_reverse_kg.fq") into (abyss_kg_pairs, velvet_kg_pairs)

	"""
	echo $forward > ${dataset_id}_read_pair_list.txt
	echo $reverse >> ${dataset_id}_read_pair_list.txt
	kmergenie "${dataset_id}_read_pair_list.txt" -t ${task.cpus} | tail -n 1 | awk '{print \$3}' > ${dataset_id}_best-k.txt
	cp $forward "${dataset_id}_forward_kg.fq"
	cp $reverse "${dataset_id}_reverse_kg.fq"
	"""
}

/*
 * Build assembly with Abyss
 */
process BuildAbyssAssembly {
	publishDir "${params.assembly_out_dir}/AbyssContigs", mode: "copy"
	tag { dataset_id }
	cache 'deep'
	
        input:
        set dataset_id, file(forward), file(reverse) from abyss_kg_pairs
	val best from best_abyss_kmer_results

	output:
        set dataset_id, file("${dataset_id}_abyss-contigs.fa") into (abyss_assembly_results, abyss_assembly_quast_contigs)

	shell:
	'''
	#!/bin/sh
	best_kmer=`cat !{best}`
	abyss-pe k=$best_kmer name=abyss j=!{task.cpus} in='!{forward} !{reverse}'
        mv abyss-contigs.fa !{dataset_id}_abyss-contigs.fa
	'''
}

/*
 * Build assembly with Velvet
 */
process BuildVelvetAssembly {
	publishDir "${params.assembly_out_dir}/VelvetContigs", mode: "copy"
	tag { dataset_id }
	cache 'deep'
	
	input:
	set dataset_id, file(forward), file(reverse) from velvet_kg_pairs
	val best from best_velvet_kmer_results

	output:
	set dataset_id, file("${dataset_id}_velvet-contigs.fa") into (velvet_assembly_results, velvet_assembly_quast_contigs)

	
	shell:
	'''
	#!/bin/sh
	best_kmer=`cat !{best}`
	velveth auto $best_kmer -separate -fastq -shortPaired !{forward} !{reverse}
	velvetg auto -exp_cov auto -cov_cutoff auto
	mv auto/contigs.fa !{dataset_id}_velvet-contigs.fa
	'''
}

/*
 * Build assembly with SPAdes
 */
process BuildSpadesAssembly {
	publishDir "${params.assembly_out_dir}/SPadesContigs", mode: "copy"
	tag { dataset_id }
	cache 'deep'

	input:
	set dataset_id, file(forward), file(reverse) from spades_read_pairs

	output:
	set dataset_id, file("${dataset_id}_spades-contigs.fa") into (spades_assembly_results, spades_assembly_quast_contigs)

	"""
	spades.py --pe1-1 ${forward} --pe1-2 ${reverse} -t ${task.cpus} -o spades_output
	mv spades_output/contigs.fasta ${dataset_id}_spades-contigs.fa
	"""
}

/*
 * Build assembly with IDBA-UD
 */
process BuildIDBAAssembly {
	publishDir "${params.assembly_out_dir}/IDBAContigs", mode: "copy"
	tag { dataset_id }
	cache 'deep'
	
	input:
	set dataset_id, file(forward), file(reverse) from idba_read_pairs

	output:
	set dataset_id, file("${dataset_id}_idba-contigs.fa") into (idba_assembly_results, idba_assembly_quast_contigs)

	"""
	fq2fa --merge --filter ${forward} ${reverse} ${dataset_id}_idba-paired-contigs.fa
	idba_ud -r ${dataset_id}_idba-paired-contigs.fa --num_threads ${task.cpus} -o ${dataset_id}_idba_output
	mv ${dataset_id}_idba_output/contig.fa ${dataset_id}_idba-contigs.fa	
	"""
}


// What's this do? Good question.
// I needed a way to group contigs produced from each assembler
// based on the reads that produced those contigs. If you don't do this
// you cannot guarantee that the assemblies passed to CISA arrived in the
// correct order.
// This function concatenates all of the tuples produced from each assembly
// (i.e., [dataset_id, dataset_id_[assembler]-contigs.fa]) and groups them
// into a single tuple based on the dataset_id. The result is the following:
// [dataset_id, abyss.fa, idba.fa, spades.fa, velvet.fa]. The order in which
// these contigs appear in the tuple is irrelevant as none of the downstream
// processes require me to know it.
abyss_assembly_results.concat(
		velvet_assembly_results,
		spades_assembly_results,
		idba_assembly_results
	)
	.groupTuple(sort: true, size: 4)
	.set { grouped_assembly_contigs }
	

/*
 * Integrate contigs produced from each assembler with CISA
 */
process IntegrateContigs {
	publishDir "${params.assembly_out_dir}/IntegratedContigs", mode: "copy"
	tag { dataset_id }
	cache 'deep'
	
	input:
	set dataset_id, file(contigs) from grouped_assembly_contigs

	output:
	set dataset_id, file("${dataset_id}_master_contigs.fa") into master_contigs
	set dataset_id, file("${dataset_id}_master_integrated_contigs.fa") into (cisa_integrated_contigs, cisa_integrated_quast_contigs)

	shell:
	'''
	#!/bin/sh
	echo count=4 > Merge.config
	echo data=!{contigs[0]},title=Contigs0 >> Merge.config
	echo data=!{contigs[1]},title=Contigs1 >> Merge.config
	echo data=!{contigs[2]},title=Contigs2 >> Merge.config
	echo data=!{contigs[3]},title=Contigs3 >> Merge.config
	echo Master_file=!{dataset_id}_master_contigs.fa >> Merge.config
	Merge.py Merge.config
	echo genome=`grep 'Whole Genome' 'Merge_info' | cut -d ':' -f2 | sort -rn | head -n 1 | tr -d [:space:]` > CISA.config
	echo infile=!{dataset_id}_master_contigs.fa >> CISA.config
	echo outfile=!{dataset_id}_master_integrated_contigs.fa >> CISA.config
	echo nucmer=`which nucmer` >> CISA.config
	echo R2_Gap=0.95 >> CISA.config
	echo CISA=${CISA} >> CISA.config
	echo makeblastdb=`which makeblastdb` >> CISA.config
	echo blastn=`which blastn` >> CISA.config
	CISA.py CISA.config
	'''
}

/*
 * Annotate the CISA integrated contigs with Prokka
 */
process AnnotateContigs {
	publishDir "${params.assembly_out_dir}/AnnotatedContigs", mode: "copy"
	tag { dataset_id }
	cache 'deep'

	input:
	set dataset_id, file(cisa_contigs) from cisa_integrated_contigs
	
	output:
	file("${dataset_id}.*") into prokka_annotations

	"""
	if [ -z ${params.species} ] && [ -z ${params.genus} ]
	then
		prokka ${cisa_contigs} --prefix ${dataset_id} --cpus ${task.cpus} --outdir annotations
	else
		prokka ${cisa_contigs} --genus ${params.genus} --species ${params.species} --centre tychus --prefix ${dataset_id} --cpus ${task.cpus} --outdir annotations

	fi
	mv annotations/* .
	"""
}

abyss_assembly_quast_contigs.concat(
                velvet_assembly_quast_contigs,
                spades_assembly_quast_contigs,
                idba_assembly_quast_contigs,
		cisa_integrated_quast_contigs
        )
        .groupTuple(sort: true, size: 5)
        .set { grouped_assembly_quast_contigs }


/*
 * Evaluate ALL assemblies with QUAST
 */
process EvaluateAssemblies {
	publishDir "${params.assembly_out_dir}/AssemblyReport", mode: "move"
	tag { dataset_id }
	cache 'deep'

	input:
	set dataset_id, file(quast_contigs) from grouped_assembly_quast_contigs

	output:
	file("${dataset_id}_*") into quast_evaluation
	
	shell:
	'''
	#!/bin/sh
	quast.py !{quast_contigs[0]} !{quast_contigs[1]} !{quast_contigs[2]} !{quast_contigs[3]} !{quast_contigs[4]} --space-efficient --threads !{params.threads} -o output
        mkdir quast_output
        find output/ -maxdepth 2 -type f | xargs mv -t quast_output
        cd quast_output
        ls * | xargs -I {} mv {} !{dataset_id}_{}
        mv * ../
	'''	
}

// Display information about the completed run
// See https://www.nextflow.io/docs/latest/metadata.html for more
// information about available onComplete options
workflow.onComplete {
	log.info "Nextflow Version:	$workflow.nextflow.version"
  	log.info "Command Line:		$workflow.commandLine"
	log.info "Container:		$workflow.container"
	log.info "Duration:		$workflow.duration"
	log.info "Output Directory:	$params.assembly_out_dir"
}
