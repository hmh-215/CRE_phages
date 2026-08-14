<h1>Genomics analysis of Carbapenem-resistant <i>Enterobacteriaceae</i> (CRE) and their targeting phages</h1>

<h2>Background</h2>
	<p class="subtitle">Carbapenem-resistant <i>Enterobacteriaceae</i> (CRE) are among the list of the WHO 2024 critical priority group, posing as among the most threatening to human health drug-resistant bacteria. On the other hand, the need for development and application of strategies for either replacing or complementary with antibiotics for the treatment of drug-resistant bacterial infection had become crucial in the current time. Phage therapy is recognized as one of the potential alternative therapeutic strategies to combat antimicrobial resistance (AMR).</p>
  <p class="subtitle">The purpose of this study is to identify the potential mechanism of phage-susceptible and phage-resistance in bacteria, as well as identification of genomic features by which phages can target CRE.</p>
<h2 id="conda-envs">Conda environments</h2>
	<p>See <code>conda_envs.txt</code> for the authoritative, versioned list. Summary of environments referenced by these scripts:</p>
	<table>
		<tr><th>Environment</th><th>Key tools</th></tr>
		<tr><td><code>preprocessing</code></td><td>Trimmomatic, FastQC, MultiQC</td></tr>
		<tr><td><code>mapping</code></td><td>BWA-MEM, SAMtools, Picard</td></tr>
		<tr><td><code>assembly</code></td><td>SPAdes, QUAST, seqkit, Unicycler</td></tr>
		<tr><td><code>BPannotation</code></td><td>Bakta, Prokka, Pharokka, PHANOTATE, eggNOG-mapper, ABRicate</td></tr>
		<tr><td><code>BPstructure</code></td><td>CheckM, CheckV</td></tr>
		<tr><td><code>BPtyping</code></td><td>mlst, Kleborate</td></tr>
		<tr><td><code>recombination</code></td><td>ISEScan, IntegronFinder, EMBOSS, PhiSpy</td></tr>
		<tr><td><code>ncbi</code></td><td>BLAST, NCBI datasets CLI, AMRFinderPlus</td></tr>
		<tr><td><code>rgi_env</code></td><td>RGI + CARD database</td></tr>
		<tr><td><code>plasmid</code></td><td>Platon, PlasmidFinder, MOB-suite</td></tr>
		<tr><td><code>phageterm_env</code></td><td>PhageTerm</td></tr>
		<tr><td><code>gipsy_env</code></td><td>Gipsy2</td></tr>
	</table>

<h2 id="repo-structure">Pipelines</h2>
	<table>
		<tr><th>File</th><th>Purpose</th></tr>
		<tr><td><code>phages_assembly_mapping_annotation.bash</code></td><td>Host-read removal, assembly, terminus prediction, annotation</td></tr>
    <tr><td><code>phages_autoBLAST.bash</code></td><td>Connection to NCBI API and automatically performing BLASTN for all phage contigs</td></tr>
    <tr><td><code>nithesis_annotation.bash</code></td><td>Mini script focusing on terminus prediction, annotation of a <i>nithesis</i>-like phage isolated</td></tr>
    <tr><td><code>Ecoli_assembly_annotation.bash</code><br><code>Ecoli_assembly_annotation_2.bash</code><td>De novo assembly and functional annotation of <i>E. coli</i> samples</td></tr>
    <tr><td><code>Kleb_assembly_annotation.bash</code><br><code>Kleb_assembly_annotation_2.bash</code></td><td>De novo assembly and functional annotation of <i>K. pneumoniae</i> samples</td></tr>
	<tr><td><code>Kleb_AMRprofiling_reassembly.bash</code></td><td>Further analysis on <i>K. pneumoniae</i> samples for AMR profiling</td></tr>
  </table>

<div class="card">
		<div class="card-title">
			<h3><code>phages_assembly_mapping_annotation.bash</code></h3>
			<span class="tag">Bacteriophage</span>
		</div>
		<div class="meta-row">
			<span><strong>Input:</strong> Illumina FASTQs, <em>E. coli</em> E72 reference assembly</span>
			<span><strong>Samples:</strong> 4 phage-bacteria mixtures (A21&ndash;A24)</span>
		</div>
		<p>Maps reads to the bacterial host to isolate unmapped (phage) reads (BWA-MEM, samtools); assembly and post-assembly QC (metaSPAdes, CheckV); termini/orientation prediction (PhageTerm); functional annotations (Pharokka, PHANOTATE, eggNOG-mapper).</p>
		<pre><code>bash ./phages_assembly_mapping_annotation.bash</code></pre>
</div>

<div class="card">
		<div class="card-title">
			<h3><code>phages_autoBLAST.bash</code></h3>
			<span class="tag">Bacteriophage</span>
		</div>
		<div class="meta-row">
			<span><strong>Input:</strong> Phages assemblies in FASTAs</span>
			<span><strong>Samples:</strong> 4 phage-bacteria mixtures (A21&ndash;A24)</span>
		</div>
		<p>Connection to NCBI API for BLASTN service with automatic contig submissions, merge all contigs per sample reports into .tsv format</p>
		<pre><code>bash ./phages_assembly_mapping_annotation.bash</code></pre>
</div>

<div class="card">
		<div class="card-title">
			<h3><code>nithesis_annotation.bash</code></h3>
			<span class="tag">Bacteriophage</span>
		</div>
		<div class="meta-row">
			<span><strong>Input:</strong> Annotated phage samples</span>
			<span><strong>Samples:</strong> Extracted for a <i>nithesis</i>-like phage genome as FASTAs</span>
		</div>
		<p>Re-checking completeness (CheckV); extraction of assembled NODEs with >90% completeness; termini/orientation prediction (PhageTerm); functional annotations (Pharokka, PHANOTATE, eggNOG-mapper).</p>
		<pre><code>bash ./phages_assembly_mapping_annotation.bash</code></pre>
</div>

<div class="card">
		<div class="card-title">
			<h3><code>Ecoli_assembly_annotation.bash</code> and <code>Ecoli_assembly_annotation_2.bash</code></h3>
			<span class="tag">Bacterial &middot; assembly</span>
  </div>
  <div class="meta-row">
			<span><strong>Input:</strong> paired-end Illumina FASTQs (<code>WS2762512A*_R1/R2.fastq.gz</code>)</span>
			<span><strong>Samples:</strong> 11 <i>E. coli</i> clinical isolates</span>
		</div>
		<p>Raw reads QC and trimming (FastQC/Trimmomatic); assembly and post-assembly QC (SPAdes, QUAST/CheckM); sequence typing (MLST); functional annotation (Bakta, AMRFinderPlus, eggNOG-mapper).</p>
		<pre><code>bash ./Ecoli_assembly_annotation*.bash</code></pre>
</div>

<div class="card">
		<div class="card-title">
			<h3><code>Kleb_assembly_annotation.bash</code> and <code>Kleb_assembly_annotation_2.bash</code></h3>
			<span class="tag">Bacterial &middot; assembly</span>
  </div>
  <div class="meta-row">
			<span><strong>Input:</strong> paired-end Illumina FASTQs (<code>WS2762512A*_R1/R2.fastq.gz</code>)</span>
			<span><strong>Samples:</strong> 9 <i>K. pneumoniae</i> clinical isolates</span>
		</div>
		<p>Raw reads QC and trimming (FastQC/Trimmomatic); assembly and post-assembly QC (SPAdes, QUAST/CheckM); sequence typing (MLST); functional annotation (Bakta, AMRFinderPlus, eggNOG-mapper).</p>
		<pre><code>bash ./Kleb_assembly_annotation*.bash</code></pre>
</div>

<div class="card">
		<div class="card-title">
			<h3><code>Kleb_AMRprofiling_reassembly.bash</code></h3>
			<span class="tag">Bacterial &middot; AMRprofilling &middot;Re-assembly</span>
  </div>
  <div class="meta-row">
			<span><strong>Input:</strong> paired-end Illumina FASTQs and filtered assemblies (<code>WS2762512A*_R1/R2.fastq.gz</code>)</span>
			<span><strong>Samples:</strong> 9 <i>K. pneumoniae</i> clinical isolates</span>
		</div>
		<p>AMR profiling (RGI, ABRicate); MGEs detection (ISEscan, intergron_finder); Re-annotation (Prokka); pathogenicity island and virulence detection (Gipsy2, PhiSpy, Kleborate); Plasmid assignment (Platon, plasmidfinder, MOBsuite); Re-mapping and chimera detection (BWA-MEM, samtools); Re-assembly for circularization detection and QC (Unicycler, QUAST, CheckM); </p>
		<pre><code>bash ./Kleb_AMRprofiling_reassembly.bash</code></pre>
</div>


<h2>Publications</h2>
<ul>
	<li>Huynh, M. H., Tran, D. Q., Pham, T. T. H., Le, T. T. H., Nguyen, C. L., Tran, T. T. T., & Nguyen, C. H. (Year). Isolation and genomic characterizations of nithesis-like bacteriophage targeting multidrug-resistant Escherichia coli from sewage water. (In press)</li>
</ul>

<h2>Notes</h2>
		Internal lab pipelines &middot; run on <code>/storage/student9/</code>
