rule target:
	input:
		uv_cc = expand("shuffledbyzip/uv_zones/{run_name}.{dataset}.zone{zone}.zonegwas.shuffledbyzip.run{run}.assoc.txt",
			run_name = "190402_RAN",
			zone = ["1", "2", "3", "5", "6", "7", "8", "9"],
			dataset = ["purebred"],
			run = list(range(1,2))),
		mv_cc = expand("shuffledbyzip/mv_zones/{run_name}.{dataset}.{mv}.zonegwas.shuffledbyzip.run{run}.assoc.txt",
			run_name = "190402_RAN",
			mv =["mv_all", "mv_3_5_7_8_9"],
			dataset = ["purebred"],
			run = list(range(1,2))),
		env = expand("shuffledbyzip/envgwas/{run_name}.{dataset}.{variable}.envgwas.shuffledbyzip.run{run}.assoc.txt",
			run_name = "190402_RAN",
			variable =["temp", "precip", "elev", "mv"],
			dataset = ["purebred"],
			run = list(range(1,2)))
		# env = expand("output/{run_name}.{dataset}.chr{chr}.gpsm.assoc.txt",
		# 	run_name = "190402_RAN",
		# 	variable =["temp", "precip", "elev"],
		# 	dataset = ["all", "purebred"],
		# 	chr = list(range(1,30)))
	shell:
		"rm .snakemake/*_tracking/*"

# rule grm:
# 	input:
# 		genotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes{run_name}.dose.mgf.gz",
# 		phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes{run_name}_casecontrol_zones.{dataset}.txt"
# 	params:
# 		oprefix = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/190326_SIM/phenotypes/{run_name}"
# 	benchmark:
# 		"benchmarks/grm/{run_name}.benchmark.txt"
# 	log:
# 		"logs/grm/{run_name}.log"
# 	output:
# 		grm = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/190326_SIM/phenotypes/output/{run_name}.sXX.txt"
# 	shell:
# 		"(gemma -g {input.genotypes} -p {input.phenotypes} -gk 2 -o {params.oprefix}) > {log}"

mv_zone_dict = {"mv_all":"1 2 3 5 6 7 8 9", "mv_3_5_7_8_9":"3 5 7 8 9"}
def mvchooser(WC):
	col_num = mv_zone_dict[WC.mv]
	return col_num

# rule mv_case_control_gwas:
# 	input:
# 		genotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/single_chrom/{run_name}.chr{chr}.dose.replaced.mgf.gz",
# 		phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_casecontrol_shuffled_zones.{dataset}.txt",
# 		#phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_casecontrol_zones.{dataset}.txt",
# 		grm = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/190402_RAN/grm/output/{run_name}.sXX.txt"
# 		#grm = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_casecontrol_zones.txt"
# 	params:
# 		cols = mvchooser,
# 		oprefix = "{run_name}.{dataset}.{mv}.chr{chr}.shuffled.casecontrol"
# 		#oprefix = "{run_name}.{dataset}.{mv}.chr{chr}.casecontrol"
# 	benchmark:
# 		"benchmarks/case_control_gwas/{run_name}.{mv}.chr{chr}.benchmark.txt"
# 	log:
# 		"logs/case_control_gwas/{run_name}.{mv}.chr{chr}.log"
# 	output:
# 		#"output/{run_name}.{dataset}.{mv}.chr{chr}.casecontrol.assoc.txt"
# 		"output/{run_name}.{dataset}.{mv}.chr{chr}.shuffled.casecontrol.assoc.txt"
# 	shell:
# 		"(gemma -g {input.genotypes} -p {input.phenotypes} -k {input.grm} -lmm 4 -n {params.cols} -o {params.oprefix}) > {log}"


env_dict = {"temp":"1", "precip":"2", "elev":"3", "mv":"1 2 3"}
def varchooser(WC):
	col_num = env_dict[WC.variable]
	return col_num

# rule env_gwas:
# 	input:
# 		genotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/single_chrom/{run_name}.chr{chr}.dose.replaced.mgf.gz",
# 		#phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_env_vars.{dataset}.txt",
# 		phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_shuffled_env_vars.{dataset}.txt",
# 		grm = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/190402_RAN/grm/output/{run_name}.sXX.txt"
# 	params:
# 		variable = varchooser,
# 		oprefix = "{run_name}.{dataset}.{variable}.chr{chr}.shuffled.envgwas"
# 		#oprefix = "{run_name}.{dataset}.{variable}.chr{chr}.envgwas"
# 	benchmark:
# 		"benchmarks/env_gwas/{run_name}.{dataset}.{variable}.chr{chr}.benchmark.txt"
# 	log:
# 		"logs/env_gwas/{run_name}.{dataset}.{variable}.chr{chr}.log"
# 	output:
# 		"output/{run_name}.{dataset}.{variable}.chr{chr}.shuffled.envgwas.assoc.txt"
# 		#"output/{run_name}.{dataset}.{variable}.chr{chr}.envgwas.assoc.txt"
# 	shell:
# 		"(gemma -g {input.genotypes} -p {input.phenotypes} -k {input.grm} -lmm 4 -n {params.variable} -o {params.oprefix}) > {log}"

rule env_gwas_shuffledbyzip:
	input:
		genotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/{run_name}.dose2.mgf.gz",
		#phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_env_vars.{dataset}.txt",
		phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_env_vars_shuffledbyzip.run{run}.{dataset}.txt",
		grm = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/190402_RAN/grm/output/{run_name}.sXX.txt"
	params:
		variable = varchooser,
		oprefix = "{run_name}.{dataset}.{variable}.envgwas.shuffledbyzip.run{run}",
		outdir = "shuffledbyzip/envgwas"
		#oprefix = "{run_name}.{dataset}.{variable}.chr{chr}.envgwas"
	benchmark:
		"benchmarks/env_gwas_shuffledbyzip/{run_name}.{dataset}.{variable}.run{run}.benchmark.txt"
	log:
		"logs/env_gwas_shuffledbyzip/{run_name}.{dataset}.{variable}.run{run}.log"
	output:
		"shuffledbyzip/envgwas/{run_name}.{dataset}.{variable}.envgwas.shuffledbyzip.run{run}.assoc.txt"
		#"output/{run_name}.{dataset}.{variable}.chr{chr}.envgwas.assoc.txt"
	shell:
		"(gemma -g {input.genotypes} -p {input.phenotypes} -k {input.grm} -lmm 4 -n {params.variable} -outdir {params.outdir} -o {params.oprefix}) > {log}"

rule mv_zone_gwas_shuffledbyzip:
	input:
		genotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/{run_name}.dose2.mgf.gz",
		#phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_env_vars.{dataset}.txt",
		phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_zones_shuffledbyzip.run{run}.{dataset}.txt",
		grm = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/190402_RAN/grm/output/{run_name}.sXX.txt"
	params:
		variable = mvchooser,
		oprefix = "{run_name}.{dataset}.{mv}.zonegwas.shuffledbyzip.run{run}",
		outdir = "shuffledbyzip/mv_zones"
		#oprefix = "{run_name}.{dataset}.{variable}.chr{chr}.envgwas"
	benchmark:
		"benchmarks/mv_zone_gwas_shuffledbyzip/{run_name}.{dataset}.{mv}.run{run}.benchmark.txt"
	log:
		"logs/mv_zone_gwas_shuffledbyzip/{run_name}.{dataset}.{mv}.run{run}.log"
	output:
		"shuffledbyzip/mv_zones/{run_name}.{dataset}.{mv}.zonegwas.shuffledbyzip.run{run}.assoc.txt"
		#"output/{run_name}.{dataset}.{variable}.chr{chr}.envgwas.assoc.txt"
	shell:
		"(gemma -g {input.genotypes} -p {input.phenotypes} -k {input.grm} -lmm 4 -n {params.variable} -outdir {params.outdir} -o {params.oprefix}) > {log}"

rule univ_zone_gwas_shuffledbyzip:
	input:
		genotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/{run_name}.dose2.mgf.gz",
		#phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_env_vars.{dataset}.txt",
		phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_zones_shuffledbyzip.run{run}.{dataset}.txt",
		grm = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/190402_RAN/grm/output/{run_name}.sXX.txt"
	params:
		variable = "{zone}",
		oprefix = "{run_name}.{dataset}.zone{zone}.zonegwas.shuffledbyzip.run{run}",
		outidr = "shuffledbyzip/uv_zones"
		#oprefix = "{run_name}.{dataset}.{variable}.chr{chr}.envgwas"
	benchmark:
		"benchmarks/univ_zone_gwas_shuffledbyzip/{run_name}.{dataset}.{zone}.run{run}.benchmark.txt"
	log:
		"logs/univ_zone_gwas_shuffledbyzip/{run_name}.{dataset}.{zone}.run{run}.log"
	output:
		"shuffledbyzip/uv_zones/{run_name}.{dataset}.zone{zone}.zonegwas.shuffledbyzip.run{run}.assoc.txt"
		#"output/{run_name}.{dataset}.{variable}.chr{chr}.envgwas.assoc.txt"
	shell:
		"(gemma -g {input.genotypes} -p {input.phenotypes} -k {input.grm} -lmm 4 -n {params.variable} -outdir {params.outidr} -o {params.oprefix}) > {log}"


# rule gpsm:
# 	input:
# 		genotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/single_chrom/{run_name}.chr{chr}.dose.mgf.gz",
# 		phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_gpsm.{dataset}.txt",
# 		grm = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/190326_SIM/output/{run_name}.sXX.txt"
# 	params:
# 		variable = "1",
# 		oprefix = "{run_name}.{dataset}.chr{chr}.gpsm"
# 	benchmark:
# 		"benchmarks/gpsm/{run_name}.{dataset}.chr{chr}.benchmark.txt"
# 	log:
# 		"logs/gpsm/{run_name}.{dataset}.chr{chr}.log"
# 	output:
# 		"output/{run_name}.{dataset}.chr{chr}.gpsm.assoc.txt"
# 	shell:
# 		"(gemma -g {input.genotypes} -p {input.phenotypes} -k {input.grm} -lmm 4 -n {params.variable} -o {params.oprefix}) > {log}"
