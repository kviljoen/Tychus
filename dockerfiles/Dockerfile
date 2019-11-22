#########################################################################################################
# Dockerfile to run the Tychus assembly AND alignment pipeline
#########################################################################################################
FROM ubuntu:16.04
MAINTAINER Chris Dean <cdean11@rams.colostate.edu>


#########################################################################################################
# Install core packages
#########################################################################################################
RUN apt-get update && apt-get install -y \
        aufs-tools \
        automake \
        build-essential \
cmake \
        wget \
        git \
        libbz2-dev \
        unzip \
        python \
	tar \
	zlib1g-dev \
  libncurses5-dev
 

#########################################################################################################
# Install Java, Python, and R
#########################################################################################################
RUN apt-get install -y \
	openjdk-8-jre \
	openjdk-8-jdk \
	python \
	r-base-dev

#########################################################################################################
# Install KmerGenie
#########################################################################################################
RUN cd /opt && \
	wget http://kmergenie.bx.psu.edu/kmergenie-1.7016.tar.gz && \
	tar -xvf kmergenie-1.7016.tar.gz && \
	cd kmergenie-1.7016 && \
	make && \
	make install && \
	cd / && \
	rm /opt/kmergenie-1.7016.tar.gz

#########################################################################################################
# Install Trimmomatic
#########################################################################################################
RUN cd /opt && \
        wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
        unzip Trimmomatic-0.36.zip && \
        cd Trimmomatic-0.36 && \
        chmod u+x trimmomatic-0.36.jar && \
        cd / && \
        rm /opt/Trimmomatic-0.36.zip

#########################################################################################################
# Install Prokka and dependencies
#########################################################################################################
RUN cd /opt && \
  	apt-get install -y libdatetime-perl libxml-simple-perl libdigest-md5-perl bioperl && \
  	git clone git://github.com/tseemann/prokka.git && \
  	prokka/bin/prokka --setupdb && \
  	ln -s /opt/prokka/bin/* /usr/bin/ && \
  	cd /

#########################################################################################################
# Install Abyss and dependencies
#########################################################################################################
RUN cd /opt && \
	apt-get install -y libsparsehash-dev libboost1.58-dev && \
	git clone https://github.com/bcgsc/abyss.git && \
	cd abyss && \
	./autogen.sh && \
	./configure --enable-maxk=160 && \
	make && \
	make install && \
	cd / && \
	rm -rf /opt/abyss

#########################################################################################################
# Install IDBA
#########################################################################################################
RUN cd /opt && \
	git clone https://github.com/cdeanj/idba.git && \
	cd idba && \
	./build.sh && \
	ln -s /opt/idba/bin/idba_ud /usr/bin && \
	ln -s /opt/idba/bin/fq2fa /usr/bin && \
	cd /


#########################################################################################################
# Install SPAdes and dependencies
#########################################################################################################
RUN cd /opt && \
	apt-get install -y cmake zlib1g-dev libbz2-dev && \
	wget http://cab.spbu.ru/files/release3.9.1/SPAdes-3.9.1.tar.gz && \
    	tar -xzf SPAdes-3.9.1.tar.gz && \
    	cd SPAdes-3.9.1 && \
	PREFIX=/usr ./spades_compile.sh && \
	cd / && \
	rm -rf /opt/SPAdes-3.9.1.tar.gz /opt/SPAdes-3.9.1

#########################################################################################################
# Install Velvet
#########################################################################################################
RUN cd /opt && \
	git clone https://github.com/dzerbino/velvet.git && \
	cd velvet && \
	make MAXKMERLENGTH=100 && \
	ln -s /opt/velvet/velvetg /opt/velvet/velveth /usr/bin && \
	cd /


#########################################################################################################
# Install CISA and dependencies
#########################################################################################################
RUN cd /opt && \
	apt-get install -y mummer ncbi-blast+ && \
	wget http://sb.nhri.org.tw/CISA/upload/en/2012/9/CISA_20120904-14032356.tar && \
	tar -xvf CISA_20120904-14032356.tar && \
	chmod 755 -R CISA1.2 && \
	ln -s /opt/CISA1.2/*.py /usr/bin && \
	cd / && \
	rm /opt/CISA_20120904-14032356.tar


#########################################################################################################
# Install QUAST and dependencies
#########################################################################################################
RUN cd /opt && \
	apt-get install -y pkg-config libfreetype6-dev libpng-dev python-matplotlib python-setuptools && \
	rm -rf /var/lib/apt/lists/* && \
	wget https://downloads.sourceforge.net/project/quast/quast-4.4.tar.gz && \
	tar -xzf quast-4.4.tar.gz && \
	cd quast-4.4 && \
	./setup.py install && \
	cd / && \
	rm -rf /opt/quast-4.4 /opt/quast-4.4.tar.gz

#########################################################################################################
# Install Velvet Optimiser
#########################################################################################################
RUN cd /opt && \
	git clone https://github.com/tseemann/VelvetOptimiser.git && \
	ln -s /opt/VelvetOptimiser/VelvetOptimiser.pl /usr/bin && \
	ln -s /opt/VelvetOptimiser/VelvetOpt /usr/bin && \
	cd /

#########################################################################################################
# Install Tbl2san
#########################################################################################################
RUN wget -N ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux64.tbl2asn.gz && \
	gunzip linux64.tbl2asn.gz && \
	chmod +x linux64.tbl2asn && \
	mv linux64.tbl2asn /opt/prokka/binaries/linux/tbl2asn && \
	cd /
	
	


#################################################################################################
# Dockerfile to run the Nextflow alignment pipeline
#################################################################################################

#################################################################################################
# Install Bowtie2
#################################################################################################
RUN git clone https://github.com/BenLangmead/bowtie2.git && \
        cd bowtie2 && \
        make NO_TBB=1 && \
        make install && \
	cd && \
	rm -rf /bowtie2

#################################################################################################
# Install HTSLib
#################################################################################################
RUN git clone https://github.com/samtools/htslib.git && \
	cd htslib && \
	autoheader && \
	autoconf && \
	./configure && \
	make && \
	make install && \
	cd / && \
	rm -rf htslib

#################################################################################################
# Install Samtools
#################################################################################################
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
	tar -xvf samtools-1.3.1.tar.bz2 && \
	cd samtools-1.3.1 && \
	make && \
	make install && \
	cd / && \
	rm -rf samtools-1.3.1 && \
	rm samtools-1.3.1.tar.bz2

#################################################################################################
# Install BCFtools
#################################################################################################
RUN wget https://github.com/samtools/bcftools/releases/download/1.3.1/bcftools-1.3.1.tar.bz2 && \
	tar -xvf bcftools-1.3.1.tar.bz2 && \
	cd bcftools-1.3.1 && \
	make && \
	make install && \
	cd / && \
	rm -rf bcftools-1.3.1 && \
	rm bcftools-1.3.1.tar.bz2

#################################################################################################
# Install Freebayes
#################################################################################################
RUN git clone --recursive git://github.com/ekg/freebayes.git && \
	cd freebayes && \
	make && \
	make install && \
	cd / && \
	rm -rf freebayes

################################################################################################
# Install Jellyfish
################################################################################################
RUN cd /opt && \
	wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.6/jellyfish-2.2.6.tar.gz && \
        tar -xvf jellyfish-2.2.6.tar.gz && \
        cd jellyfish-2.2.6 && \
        ./configure && \
        make -j 4 && \
	ln -s /opt/jellyfish-2.2.6/bin/jellyfish /usr/bin && \
	cd / && \
	rm /opt/jellyfish-2.2.6.tar.gz

################################################################################################
# Install Java
################################################################################################
RUN apt-get install -y \
	openjdk-8-jre \
        openjdk-8-jdk

################################################################################################
# Install Figtree
################################################################################################
RUN cd /opt && \
	git clone https://github.com/cdeanj/figtree.git && \
	chmod u+x figtree/lib/figtree.jar && \
	cd /


################################################################################################
# Install kSNP3 and dependencies
################################################################################################
RUN cd /opt && \
	apt-get install -y tcsh && \
	rm -rf /var/lib/apt/lists/* && \
	cd /opt && \
	git clone https://github.com/cdeanj/kSNP3.git && \
	mv kSNP3 /usr/local && \
	cd /

################################################################################################
# Install ResistomeAnalyzer
################################################################################################
RUN cd /opt && \
	git clone https://github.com/cdeanj/resistomeanalyzer.git && \
	cd resistomeanalyzer && \
	make && \
	ln -s /opt/resistomeanalyzer/resistome /usr/bin && \
	cd /

################################################################################################
# Install CoverageSampler
################################################################################################
RUN cd /opt && \
	git clone https://github.com/cdeanj/coverage_sampler.git && \
	cd coverage_sampler && \
	make && \
	ln -s /opt/coverage_sampler/csa /usr/bin && \
	cd /

RUN mkdir /SCRATCH
WORKDIR /SCRATCH
CMD ["bash"]	
