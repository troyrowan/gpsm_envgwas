rule target:
	input:
		# cc = expand("output/{run_name}.{dataset}.{mv}.chr{chr}.shuffled.casecontrol.assoc.txt",
		# 	run_name = "190402_RAN",
		# 	mv = ["mv_3_5_7_8_9"],
		# 	dataset = ["all", "purebred"],
		# 	chr = list(range(1,30))),
		# env = expand("output/{run_name}.{dataset}.{variable}.chr{chr}.shuffled.envgwas.assoc.txt",
		# 	run_name = "190402_RAN",
		# 	variable =["mv"],
		# 	dataset = ["all", "purebred"],
		# 	chr = list(range(1,30)))
		env = expand("gxe/output/{run_name}.{dataset}.chr{chr}.{variable}.gpsm.assoc.txt",
			run_name = "190402_RAN",
			variable =["temp", "precip", "elev"],
			dataset = ["purebred"],
			chr = list(range(1,30)))
	shell:
		"rm .snakemake/*_tracking/*"


env_dict = {"temp":"1", "precip":"2", "elev":"3", "mv":"1 2 3"}
def varchooser(WC):
	col_num = env_dict[WC.variable]
	return col_num


rule gpsm_gxe:
	input:
		genotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/single_chrom/{run_name}.chr{chr}.dose.replaced.mgf.gz",
		phenotypes = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/{run_name}/phenotypes/{run_name}_gpsm.{dataset}.txt",
		grm = "/CIFS/MUG01_N/deckerje/tnr343/local_adaptation_genotypes/red_angus/190402_RAN/output/{run_name}.sXX.txt"
	params:
		variable = env_dict,
		oprefix = "gxe/{run_name}.{dataset}.chr{chr}.{variable}.gpsm"
	benchmark:
		"benchmarks/gxe_gpsm/{run_name}.{dataset}.chr{chr}.{variable}.benchmark.txt"
	log:
		"logs/gxe_gpsm/{run_name}.{dataset}.chr{chr}.{variable}.log"
	output:
		"gxe/output/{run_name}.{dataset}.chr{chr}.{variable}.gpsm.assoc.txt"
	shell:
		"(gemma -g {input.genotypes} -p {input.phenotypes} -k {input.grm} -lmm 4 -n {params.variable} -o {params.oprefix}) > {log}"
